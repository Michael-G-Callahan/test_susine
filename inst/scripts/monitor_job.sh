#!/bin/bash
# monitor_job.sh — Monitor progress of a SuSiNE SLURM array job.
#
# Usage:
#   bash monitor_job.sh <job_name> <parent_job_id> [output_root]
#
# Examples:
#   bash monitor_job.sh pilot_phase_a 49012345
#   bash monitor_job.sh pilot_phase_a 49012345 /storage/work/mgc5166/SuSiNE\ -\ 2.0/test_susine/output

set -euo pipefail

# --- arguments ---
JOB_NAME="${1:?Usage: monitor_job.sh <job_name> <parent_job_id> [output_root]}"
PARENT_ID="${2:?Usage: monitor_job.sh <job_name> <parent_job_id> [output_root]}"
OUTPUT_ROOT="${3:-}"

# --- auto-detect output_root if not provided ---
if [ -z "$OUTPUT_ROOT" ]; then
  for candidate in \
    "/storage/work/mgc5166/SuSiNE - 2.0/test_susine/output" \
    "./output" \
    "../output"; do
    if [ -d "$candidate/slurm_output" ]; then
      OUTPUT_ROOT="$candidate"
      break
    fi
  done
  if [ -z "$OUTPUT_ROOT" ]; then
    echo "ERROR: Could not auto-detect output_root. Provide it as the 3rd argument."
    exit 1
  fi
fi

STAGING_DIR="${OUTPUT_ROOT}/slurm_output/${JOB_NAME}/${PARENT_ID}"
PRINTS_DIR="${OUTPUT_ROOT}/slurm_prints/${JOB_NAME}/${PARENT_ID}"
CONFIG_PATH="${OUTPUT_ROOT}/temp/${JOB_NAME}/job_config.json"

# --- validate paths ---
if [ ! -d "$STAGING_DIR" ]; then
  echo "ERROR: Staging directory not found: $STAGING_DIR"
  echo "  Has the job started writing output yet?"
  exit 1
fi

# --- helper: extract a single integer from grep output (first match only) ---
extract_int() {
  head -1 | tr -dc '0-9'
}

# --- get expected totals from config JSON ---
BUNDLES_PER_TASK=0
N_TASKS=0
if [ -f "$CONFIG_PATH" ]; then
  BUNDLES_PER_TASK=$(grep -oP '"bundles_per_task"\s*:\s*\K[0-9]+' "$CONFIG_PATH" 2>/dev/null | extract_int || echo 0)
  N_TASKS=$(grep -oP '"n_tasks"\s*:\s*\K[0-9]+' "$CONFIG_PATH" 2>/dev/null | extract_int || echo 0)
fi

# Fall back to SLURM if config didn't have n_tasks
if [ "$N_TASKS" = "0" ] || [ -z "$N_TASKS" ]; then
  N_TASKS=$(scontrol show job "$PARENT_ID" 2>/dev/null | grep -oP 'ArrayTaskId=\K[0-9]+-[0-9]+' | awk -F- '{print $2}' || echo 0)
fi
N_TASKS=${N_TASKS:-0}
BUNDLES_PER_TASK=${BUNDLES_PER_TASK:-0}

# Total expected bundles
TOTAL_BUNDLES=0
if [ "$BUNDLES_PER_TASK" -gt 0 ] && [ "$N_TASKS" -gt 0 ]; then
  TOTAL_BUNDLES=$((BUNDLES_PER_TASK * N_TASKS))
fi

# --- count task directories ---
N_TASK_DIRS=$(find "$STAGING_DIR" -maxdepth 1 -type d -name "task-*" | wc -l | tr -d ' ')

# --- count completed flushes ---
FLUSH_COUNT=$(find "$STAGING_DIR" -name "flush-*_model_metrics.csv" 2>/dev/null | wc -l | tr -d ' ')

# --- count completed tasks ---
COMPLETED_TASKS=0
if [ -d "$PRINTS_DIR" ]; then
  COMPLETED_TASKS=$(grep -rl "Completed task" "$PRINTS_DIR"/ 2>/dev/null | wc -l | tr -d ' ')
fi

# --- get job start time from SLURM ---
ELAPSED=""
ELAPSED_SEC=0
START_TIME=$(scontrol show job "$PARENT_ID" 2>/dev/null | grep -oP 'StartTime=\K[^ ]+' | head -1 || echo "")
if [ -n "$START_TIME" ] && [ "$START_TIME" != "Unknown" ]; then
  START_EPOCH=$(date -d "$START_TIME" +%s 2>/dev/null || echo "")
  if [ -n "$START_EPOCH" ]; then
    NOW_EPOCH=$(date +%s)
    ELAPSED_SEC=$((NOW_EPOCH - START_EPOCH))
    ELAPSED="$(printf '%dh %02dm %02ds' $((ELAPSED_SEC/3600)) $(((ELAPSED_SEC%3600)/60)) $((ELAPSED_SEC%60)))"
  fi
fi

# --- estimate time remaining ---
ETA=""
if [ "$TOTAL_BUNDLES" -gt 0 ] && [ "$FLUSH_COUNT" -gt 0 ] && [ "$ELAPSED_SEC" -gt 0 ]; then
  REMAINING=$((TOTAL_BUNDLES - FLUSH_COUNT))
  SEC_PER_BUNDLE=$((ELAPSED_SEC / FLUSH_COUNT))
  ETA_SEC=$((REMAINING * SEC_PER_BUNDLE))
  ETA="$(printf '~%dh %02dm' $((ETA_SEC/3600)) $(((ETA_SEC%3600)/60)))"
fi

# --- compute percentages ---
BUNDLE_PCT=""; TASK_PCT=""; TASKDIR_PCT=""
if [ "$TOTAL_BUNDLES" -gt 0 ]; then
  BUNDLE_PCT=$((FLUSH_COUNT * 100 / TOTAL_BUNDLES))
fi
if [ "$N_TASKS" -gt 0 ]; then
  TASK_PCT=$((COMPLETED_TASKS * 100 / N_TASKS))
  TASKDIR_PCT=$((N_TASK_DIRS * 100 / N_TASKS))
fi

# --- print summary ---
echo "=============================================="
echo "  SuSiNE Job Monitor: ${JOB_NAME} / ${PARENT_ID}"
echo "=============================================="
echo ""
if [ "$TOTAL_BUNDLES" -gt 0 ]; then
  echo "  Bundles completed:  ${FLUSH_COUNT} / ${TOTAL_BUNDLES} (${BUNDLE_PCT}%)"
  BAR_LEN=30
  FILLED=$((BUNDLE_PCT * BAR_LEN / 100))
  EMPTY=$((BAR_LEN - FILLED))
  BAR=$(printf '%0.s#' $(seq 1 $FILLED 2>/dev/null) || true)
  BAR="${BAR}$(printf '%0.s-' $(seq 1 $EMPTY 2>/dev/null) || true)"
  echo "  Progress:           [${BAR}] ${BUNDLE_PCT}%"
else
  echo "  Bundles completed:  ${FLUSH_COUNT} (total unknown)"
fi
if [ "$N_TASKS" -gt 0 ]; then
  echo "  Tasks completed:    ${COMPLETED_TASKS} / ${N_TASKS} (${TASK_PCT}%)"
  echo "  Task dirs created:  ${N_TASK_DIRS} / ${N_TASKS} (${TASKDIR_PCT}%)"
else
  echo "  Tasks completed:    ${COMPLETED_TASKS}"
  echo "  Task dirs created:  ${N_TASK_DIRS}"
fi
echo ""
if [ -n "$ELAPSED" ]; then
  echo "  Elapsed:            ${ELAPSED}"
fi
if [ -n "$ETA" ]; then
  echo "  Est. remaining:     ${ETA}"
fi
echo ""

# --- per-task flush histogram ---
echo "  Per-task flush distribution:"
echo "  ----------------------------"
if [ "$N_TASK_DIRS" -gt 0 ]; then
  find "$STAGING_DIR" -maxdepth 1 -type d -name "task-*" -print0 | while IFS= read -r -d '' tdir; do
    find "$tdir" -name "flush-*_model_metrics.csv" 2>/dev/null | wc -l | tr -d ' '
  done | sort -n | uniq -c | while read -r count flushes; do
    printf "    %3d tasks with %d flushes\n" "$count" "$flushes"
  done
else
  echo "    (no task directories found yet)"
fi

echo ""
echo "  Paths:"
echo "    Staging: ${STAGING_DIR}"
echo "    Prints:  ${PRINTS_DIR}"
echo "=============================================="

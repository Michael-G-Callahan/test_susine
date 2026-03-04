#!/bin/bash
# monitor_job.sh — Monitor progress of a SuSiNE SLURM array job.
#
# Usage:
#   bash monitor_job.sh <job_name> <parent_job_id> [output_root]
#
# Examples:
#   bash monitor_job.sh phase_a 49012345
#   bash monitor_job.sh phase_a 49012345 /storage/work/mgc5166/SuSiNE\ -\ 2.0/test_susine/output
#
# What it does:
#   - Counts flush files (one per completed bundle) to measure progress
#   - Shows per-task breakdown
#   - Estimates time remaining from elapsed time and completion rate
#
# Requirements:
#   - Run from the cluster where output files are accessible
#   - The job_config.json must exist in the temp directory

set -euo pipefail

# --- arguments ---
JOB_NAME="${1:?Usage: monitor_job.sh <job_name> <parent_job_id> [output_root]}"
PARENT_ID="${2:?Usage: monitor_job.sh <job_name> <parent_job_id> [output_root]}"
OUTPUT_ROOT="${3:-}"

# --- auto-detect output_root from job_config if not provided ---
if [ -z "$OUTPUT_ROOT" ]; then
  # Try common locations
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

# --- get expected totals ---
# Try to read bundles_per_task and n_tasks from the config JSON.
# Fall back to counting task directories if config isn't available.
BUNDLES_PER_TASK=""
N_TASKS=""
if [ -f "$CONFIG_PATH" ]; then
  # Extract bundles_per_task from JSON (simple grep — no jq dependency)
  BUNDLES_PER_TASK=$(grep -oP '"bundles_per_task"\s*:\s*\K[0-9]+' "$CONFIG_PATH" 2>/dev/null || true)
  N_TASKS=$(grep -oP '"n_tasks"\s*:\s*\K[0-9]+' "$CONFIG_PATH" 2>/dev/null || true)
fi

# Count task directories that exist
TASK_DIRS=$(find "$STAGING_DIR" -maxdepth 1 -type d -name "task-*" | sort)
N_TASK_DIRS=$(echo "$TASK_DIRS" | grep -c . 2>/dev/null || echo 0)

# If we couldn't get from config, estimate from SLURM
if [ -z "$N_TASKS" ] || [ "$N_TASKS" = "0" ]; then
  # Try to get array size from scontrol
  N_TASKS=$(scontrol show job "$PARENT_ID" 2>/dev/null | grep -oP 'ArrayTaskId=\K[0-9]+-[0-9]+' | awk -F- '{print $2}' || echo "")
  if [ -z "$N_TASKS" ]; then
    N_TASKS="$N_TASK_DIRS"
  fi
fi

if [ -z "$BUNDLES_PER_TASK" ] || [ "$BUNDLES_PER_TASK" = "0" ]; then
  BUNDLES_PER_TASK="?"
fi

TOTAL_BUNDLES=""
if [ "$BUNDLES_PER_TASK" != "?" ] && [ -n "$N_TASKS" ] && [ "$N_TASKS" != "0" ]; then
  TOTAL_BUNDLES=$((BUNDLES_PER_TASK * N_TASKS))
fi

# --- count completed flushes ---
FLUSH_COUNT=$(find "$STAGING_DIR" -name "flush-*_model_metrics.csv" 2>/dev/null | wc -l)

# --- count completed tasks (have "Completed task" in .out file) ---
COMPLETED_TASKS=0
if [ -d "$PRINTS_DIR" ]; then
  COMPLETED_TASKS=$(grep -l "Completed task" "$PRINTS_DIR"/*.out 2>/dev/null | wc -l || echo 0)
fi

# --- get job start time from SLURM ---
START_TIME=$(scontrol show job "$PARENT_ID" 2>/dev/null | grep -oP 'StartTime=\K[^ ]+' || echo "")
ELAPSED=""
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
if [ -n "$TOTAL_BUNDLES" ] && [ "$FLUSH_COUNT" -gt 0 ] && [ -n "$ELAPSED" ]; then
  REMAINING_BUNDLES=$((TOTAL_BUNDLES - FLUSH_COUNT))
  SEC_PER_BUNDLE=$((ELAPSED_SEC / FLUSH_COUNT))
  ETA_SEC=$((REMAINING_BUNDLES * SEC_PER_BUNDLE))
  ETA="$(printf '~%dh %02dm' $((ETA_SEC/3600)) $(((ETA_SEC%3600)/60)))"
fi

# --- print summary ---
echo "=============================================="
echo "  SuSiNE Job Monitor: ${JOB_NAME} / ${PARENT_ID}"
echo "=============================================="
echo ""
echo "  Bundles completed:  ${FLUSH_COUNT}${TOTAL_BUNDLES:+ / ${TOTAL_BUNDLES}}"
if [ -n "$TOTAL_BUNDLES" ] && [ "$TOTAL_BUNDLES" -gt 0 ]; then
  PCT=$((FLUSH_COUNT * 100 / TOTAL_BUNDLES))
  # Simple progress bar
  BAR_LEN=30
  FILLED=$((PCT * BAR_LEN / 100))
  EMPTY=$((BAR_LEN - FILLED))
  BAR=$(printf '%0.s#' $(seq 1 $FILLED 2>/dev/null) || true)
  BAR="${BAR}$(printf '%0.s-' $(seq 1 $EMPTY 2>/dev/null) || true)"
  echo "  Progress:           [${BAR}] ${PCT}%"
fi
echo "  Tasks completed:    ${COMPLETED_TASKS} / ${N_TASKS:-?}"
echo "  Task dirs created:  ${N_TASK_DIRS}"
echo ""
if [ -n "$ELAPSED" ]; then
  echo "  Elapsed:            ${ELAPSED}"
fi
if [ -n "$ETA" ]; then
  echo "  Est. remaining:     ${ETA}"
fi
echo ""

# --- per-task detail ---
echo "  Per-task flush counts:"
echo "  ----------------------"
if [ -n "$TASK_DIRS" ] && [ "$N_TASK_DIRS" -gt 0 ]; then
  echo "$TASK_DIRS" | while read -r tdir; do
    tname=$(basename "$tdir")
    tflush=$(find "$tdir" -name "flush-*_model_metrics.csv" 2>/dev/null | wc -l)
    # Check if task is done
    tid=$(echo "$tname" | sed 's/task-0*//')
    done_marker=""
    if [ -d "$PRINTS_DIR" ] && [ -f "$PRINTS_DIR/${tid}.out" ]; then
      if grep -q "Completed task" "$PRINTS_DIR/${tid}.out" 2>/dev/null; then
        done_marker=" [DONE]"
      fi
    fi
    printf "    %-12s %d flushes%s\n" "$tname" "$tflush" "$done_marker"
  done
else
  echo "    (no task directories found yet)"
fi

echo ""
echo "  Paths:"
echo "    Staging: ${STAGING_DIR}"
echo "    Prints:  ${PRINTS_DIR}"
echo "=============================================="

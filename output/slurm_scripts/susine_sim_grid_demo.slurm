#!/bin/bash
#SBATCH --job-name=susine_sim_grid_demo
#SBATCH --array=1-1
#SBATCH --time=02:59:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=mgc5166@psu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
# Log to /dev/null; we'll redirect ourselves below.
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

JOB_ROOT="C:/Users/mgcal/OneDrive/Documents/School/Research/Genetics/cTWAS Generalization/Code/test_susine/output"
CONFIG_PATH="C:/Users/mgcal/OneDrive/Documents/School/Research/Genetics/cTWAS Generalization/Code/test_susine/output/temp/susine_sim_grid_demo/job_config.json"
RUN_TASK_SCRIPT="C:/Users/mgcal/OneDrive/Documents/School/Research/Genetics/cTWAS Generalization/Code/test_susine/inst/scripts/run_task.R"
SLURM_PRINTS_BASE="C:/Users/mgcal/OneDrive/Documents/School/Research/Genetics/cTWAS Generalization/Code/test_susine/output/slurm_prints"
SLURM_OUTPUT_BASE="C:/Users/mgcal/OneDrive/Documents/School/Research/Genetics/cTWAS Generalization/Code/test_susine/output/slurm_output"
TEMP_DIR="C:/Users/mgcal/OneDrive/Documents/School/Research/Genetics/cTWAS Generalization/Code/test_susine/output/temp/susine_sim_grid_demo"
RUN_HISTORY_BASE="C:/Users/mgcal/OneDrive/Documents/School/Research/Genetics/cTWAS Generalization/Code/test_susine/output/run_history"

# --- identifiers (stable across Slurm versions) ---
JOBNAME="${SLURM_JOB_NAME}"
PARENT_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"   # parent array ID if array, else job id
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
PAD=4
SHARD_SIZE_OUTPUT=1000
SHARD_SIZE_PRINTS=1000

# --- compute shard dir for PRINTS ---
if [ "${SHARD_SIZE_PRINTS}" -gt 0 ]; then
  TI="${TASK_ID}"; if [ "${TI}" -le 0 ]; then TI=1; fi
  SHARD_IDX_PRINTS=$(( (TI - 1) / SHARD_SIZE_PRINTS ))
  SHARD_DIR_PRINTS="$(printf "shard-%03d" "${SHARD_IDX_PRINTS}")"
  PRINTS_DIR="${SLURM_PRINTS_BASE}/${JOBNAME}/${PARENT_ID}/${SHARD_DIR_PRINTS}"
else
  PRINTS_DIR="${SLURM_PRINTS_BASE}/${JOBNAME}/${PARENT_ID}"
fi
mkdir -p "${PRINTS_DIR}"

# Redirect logs (after dirs exist) — one file per task
exec >"${PRINTS_DIR}/${TASK_ID}.out" 2>"${PRINTS_DIR}/${TASK_ID}.err"
echo "[$(date -Is)] Starting task ${TASK_ID} for job ${JOBNAME} (parent ${PARENT_ID})"

# --- export job info for R to compute run-specific output dirs ---
export SUSINE_JOB_NAME="${JOBNAME}"
export SUSINE_PARENT_ID="${PARENT_ID}"
# --- site/module setup ---
module load r

export R_LIBS_USER="/storage/home/mgc5166/R/x86_64-pc-linux-gnu-library/4.3"


# --- task 1 pre-creates all shard directories and copies run_history ---
if [ "${TASK_ID}" = "1" ]; then
  FINAL_HISTORY_DIR="${RUN_HISTORY_BASE}/${JOBNAME}/${PARENT_ID}"
  mkdir -p "${FINAL_HISTORY_DIR}"
  cp "${TEMP_DIR}"/* "${FINAL_HISTORY_DIR}/"
  echo "[$(date -Is)] Task 1: copied run_history from temp to ${FINAL_HISTORY_DIR}"
  
  # Pre-create all shard directories to avoid race conditions
  SLURM_OUTPUT_JOB_DIR="${SLURM_OUTPUT_BASE}/${JOBNAME}/${PARENT_ID}"
  if [ "${SHARD_SIZE_OUTPUT}" -gt 0 ]; then
    # Read max run_id from run_table.csv (header + data)
    MAX_RUN_ID=$(tail -1 "${FINAL_HISTORY_DIR}/run_table.csv" | cut -d, -f1)
    MAX_SHARD=$(( (MAX_RUN_ID - 1) / SHARD_SIZE_OUTPUT ))
    for ((i=0; i<=MAX_SHARD; i++)); do
      SHARD_DIR="$(printf "shard-%03d" "$i")"
      mkdir -p "${SLURM_OUTPUT_JOB_DIR}/${SHARD_DIR}"
    done
    echo "[$(date -Is)] Task 1: pre-created shards 0-${MAX_SHARD} for ${MAX_RUN_ID} runs"
  fi
fi

# --- run task ---
Rscript "$RUN_TASK_SCRIPT" \
  --job-name "$JOBNAME" \
  --task-id "$TASK_ID" \
  --job-root "$JOB_ROOT" \
  --config-path "$CONFIG_PATH"

echo "[$(date -Is)] Completed task ${TASK_ID}"

#!/bin/bash
#SBATCH --job-name=effect_matching_refit_5arm
#SBATCH --time=05:59:59
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=mgc5166@psu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --account=statsresearch_sc_default
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/storage/work/mgc5166/SuSiNE - 2.0/test_susine}"
SLURM_PRINTS_BASE="${REPO_ROOT}/output/slurm_prints"

JOBNAME="${SLURM_JOB_NAME}"
PARENT_ID="${SLURM_JOB_ID}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
PRINTS_DIR="${SLURM_PRINTS_BASE}/${JOBNAME}/${PARENT_ID}"
mkdir -p "${PRINTS_DIR}"

exec >"${PRINTS_DIR}/${TASK_ID}.out" 2>"${PRINTS_DIR}/${TASK_ID}.err"
echo "[$(date -Is)] Starting ${JOBNAME} (${PARENT_ID})"

module load r

export R_LIBS_USER="${R_LIBS_USER:-/storage/home/mgc5166/R/x86_64-pc-linux-gnu-library/4.3}"
export SUSINE_DEV="1"

FRESH="${EFFECT_MATCHING_FRESH:-true}"
CHECKPOINT_EVERY="${EFFECT_MATCHING_CHECKPOINT_EVERY:-5}"
MAX_DATASETS="${EFFECT_MATCHING_MAX_DATASETS:-}"

CMD=(Rscript "${REPO_ROOT}/inst/scripts/render_effect_matching_refit_5arm.R"
  "--repo-root=${REPO_ROOT}"
  "--fresh=${FRESH}"
  "--checkpoint-every=${CHECKPOINT_EVERY}")

if [ -n "${MAX_DATASETS}" ]; then
  CMD+=("--max-datasets=${MAX_DATASETS}")
fi

echo "[$(date -Is)] Command: ${CMD[*]}"
"${CMD[@]}"

echo "[$(date -Is)] Completed ${JOBNAME}"

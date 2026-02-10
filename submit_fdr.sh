#!/usr/bin/env bash
#SBATCH --job-name=betAS_fdr
#SBATCH --array=1-4
#SBATCH --no-requeue
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/fdr_%A_%a.out
#SBATCH --error=logs/fdr_%A_%a.err

# ============================================================================
# submit_fdr.sh
# Submit betAS FDR calculations as a SLURM job array (4 comparisons)
#
# Array tasks:
#   1 = Myoblast: Taxol vs DMSO
#   2 = Myotube:  Taxol vs DMSO
#   3 = Differentiation (DMSO)
#   4 = Differentiation (Taxol)
#
# USAGE:
#   sbatch submit_fdr.sh
#
# CONTAINER: andresgordoortiz/splicing_analysis_r_crg:v1.5
# ============================================================================

set -e
set -u
set -o pipefail

# --- Configuration ---
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
CONTAINER="docker://andresgordoortiz/splicing_analysis_r_crg:v1.5"
SIF_CACHE="${WORKDIR}/singularity_cache"

# --- Create output directories ---
mkdir -p "${WORKDIR}/logs"
mkdir -p "${WORKDIR}/results"
mkdir -p "${SIF_CACHE}"

# --- Load modules ---
module load Singularity 2>/dev/null || module load singularity 2>/dev/null || true

# --- Info ---
echo "============================================================"
echo "betAS FDR Calculation — Job Array"
echo "============================================================"
echo "Job ID:        ${SLURM_JOB_ID:-local}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "Working dir:   ${WORKDIR}"
echo "Container:     ${CONTAINER}"
echo "Date:          $(date)"
echo "Node:          $(hostname)"
echo "============================================================"

# --- Pull/cache the container if not already present ---
export SINGULARITY_CACHEDIR="${SIF_CACHE}"

# --- Run the R script inside the container ---
singularity exec \
    --bind "${WORKDIR}:${WORKDIR}" \
    "${CONTAINER}" \
    Rscript "${WORKDIR}/01_fdr_calculation.R" "${SLURM_ARRAY_TASK_ID}"

echo ""
echo "============================================================"
echo "Task ${SLURM_ARRAY_TASK_ID} finished at $(date)"
echo "============================================================"

#!/bin/bash
#SBATCH --job-name=gse154386_sweep
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$HOME/CS5}"
CONFIG_FILE="${REPO_DIR}/gse154386_cogaps_sweep_config.sh"

if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "ERROR: shared config file not found: ${CONFIG_FILE}" >&2
  exit 2
fi

source "${CONFIG_FILE}"

JOBS_TSV="${JOBS_TSV_OVERRIDE:-${OUTDIR}/jobs.tsv}"
COGAPS_INPUT_H5AD="${OUTDIR}/cache/gse154386_experimental_discovery_genes_x_cells.h5ad"
DISCOVERY_H5AD="${OUTDIR}/cache/gse154386_experimental_discovery_cells_x_genes.h5ad"

echo "=== GSE154386 sweep array task ==="
date
hostname
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-}"
echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-}"
echo "SUBMIT_DIR=${SUBMIT_DIR}"
echo "OUTDIR=${OUTDIR}"
echo "JOBS_TSV=${JOBS_TSV}"
echo "=================================="

cd "${SUBMIT_DIR}"

export PYTHONPATH="${PYTHONPATH:-}"

if [[ ! -f "${JOBS_TSV}" ]]; then
  echo "ERROR: jobs.tsv not found: ${JOBS_TSV}" >&2
  exit 2
fi

if [[ ! -f "${COGAPS_INPUT_H5AD}" ]]; then
  echo "ERROR: CoGAPS input cache not found: ${COGAPS_INPUT_H5AD}" >&2
  exit 2
fi

if [[ ! -f "${DISCOVERY_H5AD}" ]]; then
  echo "ERROR: discovery cache not found: ${DISCOVERY_H5AD}" >&2
  exit 2
fi

if [[ ! -f "${CONDA_SH}" ]]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

if [[ -d "${PYCOGAPS_SOURCE}" ]]; then
  export PYTHONPATH="${PYCOGAPS_SOURCE}:${PYTHONPATH:-}"
fi

NTHREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS="${NTHREADS}"
export OPENBLAS_NUM_THREADS="${NTHREADS}"
export MKL_NUM_THREADS="${NTHREADS}"
export VECLIB_MAXIMUM_THREADS="${NTHREADS}"
export NUMEXPR_NUM_THREADS="${NTHREADS}"

LINE_NO=$((SLURM_ARRAY_TASK_ID + 1))
JOB_LINE="$(sed -n "${LINE_NO}p" "${JOBS_TSV}")"
if [[ -z "${JOB_LINE}" ]]; then
  echo "ERROR: no job line found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} (line ${LINE_NO})" >&2
  exit 2
fi

IFS=$'\t' read -r K SEED N_ITER <<< "${JOB_LINE}"
echo "Resolved job: K=${K} seed=${SEED} n_iter=${N_ITER}"

python3 gse154386_run_one_singleprocess.py \
  --outdir "${OUTDIR}" \
  --cogaps-input-h5ad "${COGAPS_INPUT_H5AD}" \
  --discovery-h5ad "${DISCOVERY_H5AD}" \
  --k "${K}" \
  --seed "${SEED}" \
  --n-iter "${N_ITER}" \
  --cogaps-threads "${NTHREADS}" \
  --blas-threads "${NTHREADS}" \
  --top-genes "${TOP_GENES}"

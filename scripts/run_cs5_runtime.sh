#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

IMAGE="${IMAGE:-othomas2/cs5-cogaps-runtime:0.1.0}"
PLATFORM="${PLATFORM:-linux/amd64}"
PORT="${PORT:-8787}"
PASSWORD="${PASSWORD:-Password12}"

CS5_SWEEP_RESULTS="${CS5_SWEEP_RESULTS:-/Users/othomas/Desktop/CS5_sweep_results}"

CS5_PRIMARY_K="${CS5_PRIMARY_K:-10}"
CS5_PRIMARY_SEED="${CS5_PRIMARY_SEED:-2}"
CS5_PRIMARY_N_ITER="${CS5_PRIMARY_N_ITER:-2000}"
CS5_PRIMARY_STEM="${CS5_PRIMARY_STEM:-cogaps_K${CS5_PRIMARY_K}_seed${CS5_PRIMARY_SEED}_iter${CS5_PRIMARY_N_ITER}}"

CS5_COMPARATOR_K="${CS5_COMPARATOR_K:-5}"
CS5_COMPARATOR_SEED="${CS5_COMPARATOR_SEED:-2}"
CS5_COMPARATOR_N_ITER="${CS5_COMPARATOR_N_ITER:-8000}"
CS5_COMPARATOR_STEM="${CS5_COMPARATOR_STEM:-cogaps_K${CS5_COMPARATOR_K}_seed${CS5_COMPARATOR_SEED}_iter${CS5_COMPARATOR_N_ITER}}"

mkdir -p \
  "${REPO_ROOT}/data/external" \
  "${REPO_ROOT}/data/processed/selected_model_k10/r" \
  "${REPO_ROOT}/data/processed/selected_model_k10/python" \
  "${REPO_ROOT}/data/processed/comparator_model_k5/r" \
  "${REPO_ROOT}/data/processed/comparator_model_k5/python" \
  "${REPO_ROOT}/data/processed/directionality" \
  "${REPO_ROOT}/data/processed/k_selection" \
  "${REPO_ROOT}/data/processed/interpretation" \
  "${REPO_ROOT}/data/processed/r_python_comparison"

echo "Launching ${IMAGE}"
echo "RStudio: http://localhost:${PORT}"
echo "User: rstudio"
echo "Password: ${PASSWORD}"
echo
echo "Mounts:"
echo "  case-study repo: ${REPO_ROOT}"
echo "  sweep results:   ${CS5_SWEEP_RESULTS}"
echo

docker run --platform "${PLATFORM}" -it --rm \
  -p "${PORT}:8787" \
  -e PASSWORD="${PASSWORD}" \
  -e CS5_PRIMARY_K="${CS5_PRIMARY_K}" \
  -e CS5_PRIMARY_SEED="${CS5_PRIMARY_SEED}" \
  -e CS5_PRIMARY_N_ITER="${CS5_PRIMARY_N_ITER}" \
  -e CS5_PRIMARY_STEM="${CS5_PRIMARY_STEM}" \
  -e CS5_COMPARATOR_K="${CS5_COMPARATOR_K}" \
  -e CS5_COMPARATOR_SEED="${CS5_COMPARATOR_SEED}" \
  -e CS5_COMPARATOR_N_ITER="${CS5_COMPARATOR_N_ITER}" \
  -e CS5_COMPARATOR_STEM="${CS5_COMPARATOR_STEM}" \
  -e CS5_CHOSEN_K="${CS5_PRIMARY_K}" \
  -e CS5_CHOSEN_SEED="${CS5_PRIMARY_SEED}" \
  -e CS5_CHOSEN_N_ITER="${CS5_PRIMARY_N_ITER}" \
  -e CS5_SWEEP_OUTDIR=/home/rstudio/project/data/external/cs5_sweep_results \
  -e CS5_PRIMARY_R_RESULTS_DIR=/home/rstudio/project/data/processed/selected_model_k10/r \
  -e CS5_PRIMARY_PYTHON_RESULTS_DIR=/home/rstudio/project/data/processed/selected_model_k10/python \
  -e CS5_COMPARATOR_R_RESULTS_DIR=/home/rstudio/project/data/processed/comparator_model_k5/r \
  -e CS5_COMPARATOR_PYTHON_RESULTS_DIR=/home/rstudio/project/data/processed/comparator_model_k5/python \
  -e CS5_DIRECTIONALITY_DIR=/home/rstudio/project/data/processed/directionality \
  -e CS5_K_SELECTION_DIR=/home/rstudio/project/data/processed/k_selection \
  -e CS5_INTERPRETATION_DIR=/home/rstudio/project/data/processed/interpretation \
  -e CS5_R_PYTHON_COMPARISON_DIR=/home/rstudio/project/data/processed/r_python_comparison \
  -v "${REPO_ROOT}:/home/rstudio/project:rw" \
  -v "${CS5_SWEEP_RESULTS}:/home/rstudio/project/data/external/cs5_sweep_results:ro" \
  "${IMAGE}"

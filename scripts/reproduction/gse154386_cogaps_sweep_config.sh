#!/usr/bin/env bash

# Shared configuration for the Case Study 5 CoGAPS Slurm sweep.
# Edit values here once; the sbatch wrappers source this file.

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Cluster / environment defaults for the same Rhino setup as the reference sweep.
CONDA_SH="${CONDA_SH:-/app/software/Miniforge3/24.1.2-0/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-oshane-jlab}"
PYCOGAPS_SOURCE="${PYCOGAPS_SOURCE:-${HOME}/src/pycogaps}"

# Sweep location defaults.
# SUBMIT_DIR can safely be fixed now because you said the run target is /home/othomas/CS5.
SUBMIT_DIR="${SUBMIT_DIR:-/home/othomas/CS5}"
OUTDIR="${OUTDIR:-GSE154386/cogaps_sweep_singleprocess_hpc}"

# Sweep grid defaults.
# These are sweep design choices, not values present in the original single-run script.
# The current default keeps the original dense low-K sweep and the approved high-K expansion.
K_GRID="${K_GRID:-5,6,7,8,9,10,11,12,16,20,24,28,32,36,40}"
SEEDS="${SEEDS:-1,2,3,4,5}"
ITERS="${ITERS:-2000,4000,6000,8000}"

# These come directly from the local non-distributed Case Study 5 script.
DISCOVERY_MAX_CELLS_PER_SAMPLE="${DISCOVERY_MAX_CELLS_PER_SAMPLE:-200}"
DISCOVERY_SEED="${DISCOVERY_SEED:-42}"
TOP_GENES="${TOP_GENES:-50}"

# Aggregation defaults introduced for the sweep.
PLATEAU_TOL="${PLATEAU_TOL:-0.03}"
MIN_SUCCESSFUL_RUNS="${MIN_SUCCESSFUL_RUNS:-3}"

# Cluster throughput default copied from the reference sweep example.
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-15}"

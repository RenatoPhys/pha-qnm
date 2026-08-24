#!/bin/sh
#SBATCH --job-name=pha-qnm
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --output=results/manifests/slurm-%j.log

set -eu
: "${PHA_QNM_BIN:?Set PHA_QNM_BIN to the hpc-slepc executable}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
srun "$PHA_QNM_BIN" qnm scan --config configs/runs/full_scan.yaml --resume


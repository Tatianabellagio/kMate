#!/bin/bash
#SBATCH --job-name=g9win_diag
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=0:40:00
#SBATCH --requeue
#SBATCH --output=logs/g9win_diag_%A_%a.out
#SBATCH --error=logs/g9win_diag_%A_%a.err
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
$PY analysis/grenenet_gea/gen9_window/founder_vs_evolved_dynld.py "Chr${SLURM_ARRAY_TASK_ID}"

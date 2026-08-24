#!/bin/bash
#SBATCH --job-name=driftnull_lin
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --output=logs/driftnull_%j.out
#SBATCH --error=logs/driftnull_%j.out
# WF-drift-through-projection FWER null for the sampling-variance block LD-LMM (REML-refit per sim).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
SITE=${SITE:-4} NSIM=${NSIM:-2000}
echo "[$(date)] driftnull SITE=$SITE NSIM=$NSIM on $(hostname)"
SITE=$SITE NSIM=$NSIM $PY analysis/grenenet_gea/block_ld_lmm_driftnull_linear.py
echo "[$(date)] done"

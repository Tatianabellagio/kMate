#!/bin/bash
#SBATCH --job-name=gen1_allsites
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=1:00:00
#SBATCH --array=0-29
#SBATCH --output=logs/gen1_allsites_%A_%a.out
#SBATCH --error=logs/gen1_allsites_%A_%a.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
SITES=(1 2 4 5 6 9 10 11 12 13 23 24 25 26 27 28 32 37 42 43 45 46 48 49 52 53 54 55 57 60)
SITE=${SITES[$SLURM_ARRAY_TASK_ID]}
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] gen1 SITE=$SITE on $(hostname)"
SITE=$SITE $PY analysis/grenenet_gea/block_ld_lmm_gen1.py
echo "[$(date)] done SITE=$SITE"

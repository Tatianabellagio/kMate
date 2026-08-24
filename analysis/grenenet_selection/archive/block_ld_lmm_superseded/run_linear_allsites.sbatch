#!/bin/bash
#SBATCH --job-name=lin_allsites
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=1:00:00
#SBATCH --array=0-16
#SBATCH --output=logs/lin_allsites_%A_%a.out
#SBATCH --error=logs/lin_allsites_%A_%a.out
# Per-site LINEAR + sampling-floor block LD-LMM across all 17 usable GrENE-net sites.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
SITES=(1 4 5 9 11 13 27 28 32 42 45 46 52 53 54 55 57)
SITE=${SITES[$SLURM_ARRAY_TASK_ID]}
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] linear_sampvar SITE=$SITE on $(hostname)"
SITE=$SITE $PY analysis/grenenet_gea/block_ld_lmm_linear_sampvar.py 2>&1 | grep -E "LINEAR|hits|ONE-SIDED" | head -4
echo "[$(date)] done SITE=$SITE"

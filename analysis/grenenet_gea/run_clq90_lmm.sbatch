#!/bin/bash
#SBATCH --job-name=clq90_lmm
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=1:30:00
#SBATCH --output=logs/clq90_lmm_%A_%a.out
#SBATCH --error=logs/clq90_lmm_%A_%a.out
# Calibrated temporal LD-LMM (drift-null / linked-background) on the clq0.9 hapfreq
# matrix, per site, then the size-matched SV-enrichment vs the FDR-selected blocks.
# Only sites with >=3 persistent plots (LMM needs replicate plots for its per-plot SE).
#   sbatch --array=0-16 run_clq90_lmm.sbatch
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
SITES=(27 1 42 46 55 52 57 53 54 9 11 5 45 28 4 32 13)
SITE=${SITES[$SLURM_ARRAY_TASK_ID]}
KPY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
export HF_DIR=results/grenenet_gea/hapfreq_clq90 MEMB_TAG=clq90 SITE=$SITE
echo "[$(date)] host=$(hostname) site=$SITE  clq90 LD-LMM"
$KPY -u analysis/grenenet_gea/block_ld_lmm_temporal.py
echo "[$(date)] site=$SITE  SV-enrichment vs LMM-selected blocks"
$KPY -u analysis/grenenet_gea/sv_enrichment_lmm.py --site "$SITE" --qcut 0.05 --n-perm 5000
echo "[$(date)] DONE site=$SITE"

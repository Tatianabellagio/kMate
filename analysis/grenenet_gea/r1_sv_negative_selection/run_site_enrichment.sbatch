#!/bin/bash
#SBATCH --job-name=sv_enrich_site
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/sv_enrich_%A_%a.out
#SBATCH --error=logs/sv_enrich_%A_%a.out
# Per-site temporal SV-enrichment on clq0.9 haploblocks (GLOBAL-mode AF).
# Array index -> site. Each task computes per-variant s (in-memory) + the clq0.9
# block enrichment (size-matched null on n_nonsv). Small outputs only (json/parquet/npz).
#   sbatch --array=0-16 run_site_enrichment.sbatch
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
# all 31 sites with >=1 post-founding generation (variable-length pooled trajectories).
# 19 have gen1-3, 8 have 2 gens, 4 are single-gen Δp-only (23 26 37 60).
SITES=(1 2 4 5 6 9 10 11 12 13 23 24 25 26 27 28 32 33 37 42 43 45 46 48 49 52 53 54 55 57 60)
SITE=${SITES[$SLURM_ARRAY_TASK_ID]}
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] host=$(hostname) site=$SITE"
MIN_MAC=${MIN_MAC:-12}  # founder-panel MAC floor; 12=common(MAF5%, current choice), 2=no-singletons(~MAF1%)
$PY -u analysis/grenenet_gea/r1_sv_negative_selection/site_sv_enrichment.py --site "$SITE" --top-frac 0.01 --n-perm 5000 --min-mac "$MIN_MAC"
echo "[$(date)] DONE site=$SITE"

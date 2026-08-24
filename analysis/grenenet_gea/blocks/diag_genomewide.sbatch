#!/bin/bash
#SBATCH --job-name=diag_gw
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/diag_gw_%A_%a.out
#SBATCH --error=logs/diag_gw_%A_%a.out
# Fill genome-wide PRE-FILTER diagnostics on the OLD 0.5 map (blocks_recompute, MINCF=0.5):
#  - missing-data sensitivity (Part E) for Chr2-5
#  - panel-support tag for Chr5 (Part D gap; chr1-4 already on disk)
# These justify the panel-support filter decision, so they run on the unfiltered map.

set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
CHRS=(Chr2 Chr3 Chr4 Chr5)
CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python

echo "[$(date)] $CHR : missing-data sensitivity (ARI imputed vs masked)"
$PY -u analysis/grenenet_gea/blocks/block_missing_sensitivity.py --chrom "$CHR"

if [ "$CHR" = "Chr5" ]; then
  echo "[$(date)] Chr5 : panel-support tag (Part D gap)"
  $PY -u analysis/grenenet_gea/blocks/block_panel_support_tag.py --chrom Chr5
fi
echo "[$(date)] $CHR done"

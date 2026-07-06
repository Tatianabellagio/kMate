#!/bin/bash
#SBATCH --job-name=dynld
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=4:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/dynld_%A_%a.out
#SBATCH --error=logs/dynld_%A_%a.out
# Dynamic LD-guided block growth (post-processing on CLQ0.9 blocks; HapFM untouched).
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
CHRS=(Chr1 Chr2 Chr3 Chr4 Chr5); CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] $CHR dynamic LD merge"
$PY -u analysis/grenenet_gea/dynamic_ld_blocks.py --chrom "$CHR" --K 500 \
  2>&1 | grep -vE "QUALITY CHECK|window_size|variants kept"
echo "[$(date)] $CHR done"

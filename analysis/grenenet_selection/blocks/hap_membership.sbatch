#!/bin/bash
#SBATCH --job-name=hapmemb
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=4:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/hapmemb_%A_%a.out
#SBATCH --error=logs/hapmemb_%A_%a.out
# STEP 1: founder->haplotype membership (projection matrix M_b) on the dynld-K500 units.
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/grenenet_selection/blocks_mcf90/hap_membership
CHRS=(Chr1 Chr2 Chr3 Chr4 Chr5); CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] $CHR haplotype membership build"
$PY -u analysis/grenenet_selection/blocks/build_hap_membership.py --chrom "$CHR" --K 500 \
  2>&1 | grep -vE "variants kept|window_size|QUALITY"
echo "[$(date)] $CHR done"

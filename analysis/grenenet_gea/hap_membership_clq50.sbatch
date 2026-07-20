#!/bin/bash
#SBATCH --job-name=hapmemb_clq50
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=4:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/hapmemb_clq50_%A_%a.out
#SBATCH --error=logs/hapmemb_clq50_%A_%a.out
# founder->haplotype membership on the r2>=0.5 (clq0.5) LD blocks (coarser than clq0.9).
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/grenenet_gea/blocks_mcf90/hap_membership
CHRS=(Chr1 Chr2 Chr3 Chr4 Chr5); CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}
CHRLC=$(echo "$CHR" | tr 'A-Z' 'a-z')
UNITS="analysis/grenenet_gea/blocks_mcf90/${CHRLC}_clq0.5_blocks_clq0.5.tsv"
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] $CHR clq0.5 haplotype membership build from $UNITS"
$PY -u analysis/grenenet_gea/build_hap_membership.py --chrom "$CHR" --units "$UNITS" --tag clq50 \
  2>&1 | grep -vE "variants kept|window_size|QUALITY"
echo "[$(date)] $CHR done"

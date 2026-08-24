#!/bin/bash
#SBATCH --job-name=hapmemb_clq90
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=4:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/hapmemb_clq90_%A_%a.out
#SBATCH --error=logs/hapmemb_clq90_%A_%a.out
# founder->haplotype membership on the r2>=0.9 (clq0.9) LD blocks instead of dynld-K500.
# Rationale: the founder GWAS genotype is panel haplotype membership and the trait is
# genome-wide ecotype h -- neither estimates a block frequency from pool-seq, so the
# k-mer-coverage motivation for the coarse dynld-K500 units does not apply; finer clq0.9
# units give better localization. Identical HapFM-xmeans clustering rule, just finer blocks.
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/grenenet_selection/blocks_mcf90/hap_membership
CHRS=(Chr1 Chr2 Chr3 Chr4 Chr5); CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}
CHRLC=$(echo "$CHR" | tr 'A-Z' 'a-z')
UNITS="analysis/grenenet_selection/blocks_mcf90/${CHRLC}_clq0.9_blocks_clq0.9.tsv"
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] $CHR clq0.9 haplotype membership build from $UNITS"
$PY -u analysis/grenenet_selection/blocks/build_hap_membership.py --chrom "$CHR" --units "$UNITS" --tag clq90 \
  2>&1 | grep -vE "variants kept|window_size|QUALITY"
echo "[$(date)] $CHR done"

#!/bin/bash
#SBATCH --job-name=hapmemb_clq90nosv
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=4:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/hapmemb_clq90nosv_%A_%a.out
#SBATCH --error=logs/hapmemb_clq90nosv_%A_%a.out
# SELF-TAGGING CONTROL: identical clq0.9 haplotype membership build, but --exclude-sv drops SV
# columns (|alt_len-ref_len|>50) from the founder clustering matrix so the clq90 haplotypes are
# defined by SNP+indel only. The SVs later tested for tagging these haplotypes then no longer
# helped DEFINE them -> removes the circularity in the surviving z~2-3 ecotype-axis signal.
# Output: chr{N}_hapmemb_clq90nosv.npz (originals chr{N}_hapmemb_clq90.npz untouched).
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/grenenet_selection/blocks_mcf90/hap_membership
CHRS=(Chr1 Chr2 Chr3 Chr4 Chr5); CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}
CHRLC=$(echo "$CHR" | tr 'A-Z' 'a-z')
UNITS="analysis/grenenet_selection/blocks_mcf90/${CHRLC}_clq0.9_blocks_clq0.9.tsv"
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] $CHR clq0.9 SNP+indel-only haplotype membership build from $UNITS"
$PY -u analysis/grenenet_selection/blocks/build_hap_membership.py --chrom "$CHR" --units "$UNITS" \
  --tag clq90nosv --exclude-sv \
  2>&1 | grep -vE "variants kept|window_size|QUALITY"
echo "[$(date)] $CHR done"

#!/bin/bash
#SBATCH --job-name=overlap_t2_gtchk
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/tier2_gtcheck_%j.out
#SBATCH --error=logs/tier2_gtcheck_%j.err
mkdir -p logs
set -euo pipefail

# Tier 2 (substitute for k-mer Jaccard): pairwise SNP-genotype identity
# between the 53 extras and 82 GrENE-Net cactus, via bcftools gtcheck on
# biallelic SNPs from pang_135 raw VCF.
#
# Quasi-duplicates of any 82-cactus accession contribute no new genotyping
# information to the 153 PG founders. This is the same metric that caught
# 5772/6150 originally (99.81% SNP identity).

source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh
conda activate sequencing_pipeline

WORK=/global/scratch/users/tbellg/kmate/panel_overlap_135_vs_82
RAW=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.raw.vcf.gz

echo "[$(date +%H:%M:%S)] building biallelic-SNP-only subset for gtcheck"
SNPS=$WORK/data/pang135_biallelic_snps.vcf.gz
if [ ! -s "$SNPS" ]; then
  bcftools view --threads 4 -v snps -m2 -M2 -Oz -o "$SNPS" "$RAW"
  tabix -p vcf -f "$SNPS"
fi
echo "[$(date +%H:%M:%S)] biallelic SNPs file: $(stat -c%s $SNPS) bytes"
echo "[$(date +%H:%M:%S)] n records: $(bcftools index -n $SNPS)"

echo "[$(date +%H:%M:%S)] running gtcheck on 53 extras vs 82 cactus..."

# Build a 53 x 82 pairs file: for each (extras, cactus) emit one line.
PAIRS=$WORK/data/gtcheck_pairs.tsv
> $PAIRS
while read e; do
  while read c; do
    printf "%s\t%s\n" "$e" "$c" >> $PAIRS
  done < $WORK/data/grenenet_82_in_pang135.txt
done < $WORK/data/extras_53.txt
echo "[$(date +%H:%M:%S)] pairs file: $(wc -l < $PAIRS) pairs"

OUT=$WORK/results/tier2_gtcheck.tsv
bcftools gtcheck --use GT --no-HWE-prob -P $PAIRS -O t -o $OUT $SNPS
echo "[$(date +%H:%M:%S)] gtcheck done, output: $OUT"
echo "[$(date +%H:%M:%S)] head:"
head -20 $OUT

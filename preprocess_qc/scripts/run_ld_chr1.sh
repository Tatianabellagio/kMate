#!/bin/bash
#SBATCH --job-name=ld_chr1
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/preprocess_qc/logs/ld_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/preprocess_qc/logs/ld_%j.err

# =============================================================================
# Compute pairwise LD on Chr1 for both panels:
#   1. GrENE-Net 231-founder SNP-only (greneNet_final_v1.1.recode.vcf.gz, chrom='1')
#   2. Production merged 231-founder (founders_231_chr.vcf.gz, chrom='Chr1')
#       — separately for SNP anchors and SV anchors
#
# Outputs land in preprocess_qc/output/ld/
# =============================================================================
set -eo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
OUT=$BASE/preprocess_qc/output/ld
mkdir -p $OUT

GN_VCF=$BASE/data/vcf/greneNet_final_v1.1.recode.vcf.gz
MERGED_VCF=$BASE/pangenie_genotyping/data/merged/founders_231_chr.vcf.gz

echo "[$(date)] LD on GrENE-Net SNPs (chr1, anchor=SNP)"
$PY $BASE/preprocess_qc/scripts/compute_ld.py \
    --vcf $GN_VCF --chrom 1 --out $OUT/grenenet_snp_chr1.tsv \
    --n-anchors 5000 --max-dist 1000000 --min-maf 0.05 \
    --anchor-classes SNP

echo "[$(date)] LD on production merged VCF (chr1, anchor=SNP)"
$PY $BASE/preprocess_qc/scripts/compute_ld.py \
    --vcf $MERGED_VCF --chrom Chr1 --out $OUT/merged_snp_chr1.tsv \
    --n-anchors 5000 --max-dist 1000000 --min-maf 0.05 \
    --anchor-classes SNP

echo "[$(date)] LD on production merged VCF (chr1, anchor=SVs only — SV-x-anything)"
$PY $BASE/preprocess_qc/scripts/compute_ld.py \
    --vcf $MERGED_VCF --chrom Chr1 --out $OUT/merged_sv_chr1.tsv \
    --n-anchors 5000 --max-dist 1000000 --min-maf 0.05 \
    --anchor-classes small_sv,medium_sv,large_sv

echo "[$(date)] DONE"
ls -lh $OUT/

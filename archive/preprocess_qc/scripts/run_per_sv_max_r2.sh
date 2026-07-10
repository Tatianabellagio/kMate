#!/bin/bash
#SBATCH --job-name=per_sv_r2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --array=0-9%8
#SBATCH --requeue
#SBATCH --output=logs/per_sv_r2_%A_%a.out
#SBATCH --error=logs/per_sv_r2_%A_%a.err

# Per-SV max-r² to nearest SNP (±50 kb), all 231 founders.
# SVs are sourced from the production VCF founders_231_chr.vcf.gz (80 cactus +
# 151 PanGenie, no Beagle imputation), with multi-allelics decomposed inline.
#
#   tasks 0..4 → within   Chr{1..5}   (SNPs from production VCF)
#   tasks 5..9 → cross    Chr{1..5}   (SNPs from xwu's GrENE-Net SNP-only)

set -eo pipefail
BASE=/global/scratch/users/tbellg/hapfire_sv
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
OUT=$BASE/preprocess_qc/output/ld
mkdir -p $OUT $BASE/preprocess_qc/logs

PROD_VCF=$BASE/pangenie_genotyping/data/merged/founders_231_chr.vcf.gz
GN_VCF=$BASE/archive/hapfire_projection/data/vcf/greneNet_final_v1.1.recode.vcf.gz

T=$SLURM_ARRAY_TASK_ID
if [[ $T -lt 5 ]]; then
    MODE=within
    CHR_NUM=$((T + 1))
    EXTRA=""
else
    MODE=cross
    CHR_NUM=$((T - 4))
    EXTRA="--snp-vcf $GN_VCF --snp-chrom ${CHR_NUM}"
fi

PROD_CHR="Chr${CHR_NUM}"
OUTFILE=$OUT/per_sv_max_r2_${MODE}_chr${CHR_NUM}.tsv
echo "[$(date)] task=$T mode=$MODE chr=$PROD_CHR"
echo "  out: $OUTFILE"

$PY $BASE/preprocess_qc/scripts/compute_per_sv_max_r2.py \
    --prod-vcf "$PROD_VCF" \
    --mode "$MODE" --prod-chrom "$PROD_CHR" \
    --out "$OUTFILE" \
    --window-bp 50000 --min-maf 0.05 --min-call-frac 0.7 \
    $EXTRA

gzip -f "$OUTFILE"
echo "[$(date)] done -> ${OUTFILE}.gz"

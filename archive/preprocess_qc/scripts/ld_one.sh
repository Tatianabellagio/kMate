#!/bin/bash
#SBATCH --job-name=ld_arr
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/ld_%A_%a.out
#SBATCH --error=logs/ld_%A_%a.err

# =============================================================================
# ld_one.sh — array task computes LD on one chromosome × one panel × one
# anchor-class set. Array index encodes (chrom, panel) via index % 15:
#   0..4   : GrENE-Net SNP anchors, chr 1..5
#   5..9   : production SNP anchors, chr 1..5
#   10..14 : production SV anchors, chr 1..5
# =============================================================================
mkdir -p logs
set -eo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
OUT=$BASE/preprocess_qc/output/ld
mkdir -p $OUT

GN_VCF=$BASE/data/vcf/greneNet_final_v1.1.recode.vcf.gz
MERGED_VCF=$BASE/pangenie_genotyping/data/merged/founders_231_chr.vcf.gz

IDX=${SLURM_ARRAY_TASK_ID:?must be set 0..14}
GROUP=$((IDX / 5))   # 0=grenenet,SNP / 1=merged,SNP / 2=merged,SV
CHR=$((IDX % 5 + 1)) # 1..5

case $GROUP in
    0) VCF=$GN_VCF;     CHROM=$CHR;       ANCHOR=SNP;                          OUTNAME=grenenet_snp ;;
    1) VCF=$MERGED_VCF; CHROM=Chr$CHR;    ANCHOR=SNP;                          OUTNAME=merged_snp ;;
    2) VCF=$MERGED_VCF; CHROM=Chr$CHR;    ANCHOR=small_sv,medium_sv,large_sv;  OUTNAME=merged_sv ;;
    *) echo "ERROR: bad group $GROUP" >&2; exit 1 ;;
esac

OUT_TSV=$OUT/${OUTNAME}_chr${CHR}.tsv
if [ -s "$OUT_TSV" ]; then
    echo "[$(date)] $OUT_TSV exists, skipping"; exit 0
fi

echo "[$(date)] LD: vcf=$VCF chrom=$CHROM anchor=$ANCHOR -> $OUT_TSV"
$PY $BASE/preprocess_qc/scripts/compute_ld.py \
    --vcf $VCF --chrom $CHROM --out $OUT_TSV \
    --n-anchors 5000 --max-dist 1000000 --min-maf 0.05 \
    --anchor-classes $ANCHOR

echo "[$(date)] DONE  $(wc -l < $OUT_TSV) lines"

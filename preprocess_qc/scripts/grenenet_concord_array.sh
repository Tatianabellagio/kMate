#!/bin/bash
#SBATCH --job-name=gn_concord
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/preprocess_qc/logs/gn_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/preprocess_qc/logs/gn_%A_%a.err

# =============================================================================
# grenenet_concord_array.sh — array runner for GrENE-Net SNP concordance.
# Compares PanGenie genotype VCFs against the independent GrENE-Net 231-founder
# SNP catalog. PANEL=main|loo selects the manifest + genotyped dir.
# =============================================================================
set -eo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PANG=$BASE/pangenie_genotyping
QC=$BASE/preprocess_qc

PANEL=${PANEL:-main}
if [ "$PANEL" = "main" ]; then
    MANIFEST=$PANG/data/ena_manifest.tsv
    GT_DIR=$PANG/data/genotyped
    ECO_COL=2
elif [ "$PANEL" = "loo" ]; then
    MANIFEST=$PANG/data/loo_ena_manifest.tsv
    GT_DIR=$PANG/data/loo_genotyped
    ECO_COL=1
else
    echo "ERROR: PANEL must be main or loo" >&2; exit 1
fi

OUT_DIR=$QC/output/grenenet_concordance
mkdir -p $OUT_DIR

TRUTH_VCF=$BASE/data/vcf/greneNet_final_v1.1.recode.vcf.gz

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $MANIFEST)
[ -n "$LINE" ] || { echo "ERROR: no row at $IDX" >&2; exit 1; }
ECOTYPE=$(echo "$LINE" | cut -f$ECO_COL)
PG_VCF=$GT_DIR/${ECOTYPE}_genotyping.vcf.gz

if [ ! -f "$PG_VCF" ]; then
    echo "ERROR: $PG_VCF not found" >&2; exit 1
fi
if [ -f "$OUT_DIR/${PANEL}_${ECOTYPE}_summary.tsv" ]; then
    echo "[$(date)] $ECOTYPE: already done, skipping"; exit 0
fi

PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
$PY $QC/scripts/grenenet_snp_concordance.py \
    --pangenie-vcf $PG_VCF \
    --truth-vcf    $TRUTH_VCF \
    --sample       $ECOTYPE \
    --panel        $PANEL \
    --out          $OUT_DIR/${PANEL}_${ECOTYPE}

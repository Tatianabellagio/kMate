#!/bin/bash
#SBATCH --job-name=qc_gt
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --requeue
#SBATCH --output=logs/qc_gt_%A_%a.out
#SBATCH --error=logs/qc_gt_%A_%a.err

# =============================================================================
# qc_genotyped_array.sh — array runner for per-sample genotyped-VCF QC.
# Reads ena_manifest.tsv (or loo_ena_manifest.tsv via PANEL=loo) and runs
# qc_genotyped_one.py on the corresponding <eco>_genotyping.vcf.gz.
# =============================================================================
set -eo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv
PANGENIE=$BASE/pangenie_genotyping
QCBASE=$BASE/preprocess_qc

PANEL=${PANEL:-main}
if [ "$PANEL" = "main" ]; then
    MANIFEST=$PANGENIE/data/ena_manifest.tsv
    GT_DIR=$PANGENIE/data/genotyped
    ECO_COL=2
elif [ "$PANEL" = "loo" ]; then
    MANIFEST=$PANGENIE/data/loo_ena_manifest.tsv
    GT_DIR=$PANGENIE/data/loo_genotyped
    ECO_COL=1
else
    echo "ERROR: PANEL must be main or loo" >&2; exit 1
fi

OUT_DIR=$QCBASE/output/genotyped_qc
mkdir -p $OUT_DIR $QCBASE/logs

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $MANIFEST)
[ -n "$LINE" ] || { echo "ERROR: no row at $IDX" >&2; exit 1; }
ECOTYPE=$(echo "$LINE" | cut -f$ECO_COL)

VCF=$GT_DIR/${ECOTYPE}_genotyping.vcf.gz
OUT=$OUT_DIR/${PANEL}_${ECOTYPE}.tsv
if [ ! -f "$VCF" ]; then echo "ERROR: $VCF not found" >&2; exit 1; fi

PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
$PY $QCBASE/scripts/qc_genotyped_one.py \
    --vcf $VCF \
    --sample $ECOTYPE \
    --out $OUT

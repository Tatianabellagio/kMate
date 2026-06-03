#!/bin/bash
#SBATCH --job-name=kmate_phaseC_cat
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:30:00
#SBATCH --requeue
#SBATCH --output=logs/kmate_phaseC_%A_%a.out
#SBATCH --error=logs/kmate_phaseC_%A_%a.err

# ============================================================================
# TWO-PHASE COHORT RUNNER — PHASE C: concat per-chrom TSVs + drop the DB.
#
# ⚠️  OPTIONAL / NOT YET RUN AT SCALE — single-phase run_site_array_perchrom.sh
# is the validated production path. See grenenet/README.md.
#
# One array task per sample. Concatenates ${SAMPLE}_${CHR}.tsv (all CHROMS) into
# the genome-wide ${SAMPLE}.tsv (header kept once), then removes ${SAMPLE}.jf to
# reclaim scratch. Idempotent: skips if the final TSV already exists.
#
# Usage:
#   sbatch --array=1-N%C --dependency=afterok:<phaseB_jobid> \
#     --export=ALL,MANIFEST=...,DB_DIR=...,OUT_DIR=... \
#     grenenet/run_phaseC_concat.sh
#
# Required env: MANIFEST, DB_DIR, OUT_DIR
# Optional env: CHROMS ("Chr1 Chr2 Chr3 Chr4 Chr5"), KEEP_DB=1 (don't delete .jf)
# ============================================================================
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs

: ${MANIFEST:?Set MANIFEST}
: ${DB_DIR:?Set DB_DIR}
: ${OUT_DIR:?Set OUT_DIR}
CHROMS=${CHROMS:-"Chr1 Chr2 Chr3 Chr4 Chr5"}

OFFSET=${OFFSET:-0}
LINE=$((SLURM_ARRAY_TASK_ID + OFFSET + 1))
ROW=$(awk -F'\t' -v n=$LINE 'NR==n' "$MANIFEST")
SAMPLE=$(echo "$ROW" | cut -f1)

FINAL_OUT=${OUT_DIR}/${SAMPLE}.tsv
if [ -s "$FINAL_OUT" ]; then
    echo "[$(date)] $SAMPLE already concatenated — (re)moving DB only"
    [ -n "${KEEP_DB:-}" ] || rm -f "${DB_DIR}/${SAMPLE}.jf"
    exit 0
fi

TMP_OUT=${FINAL_OUT}.tmp
: > "$TMP_OUT"
first=1
for CHR in $CHROMS; do
    f=${OUT_DIR}/${SAMPLE}_${CHR}.tsv
    [ -s "$f" ] || { echo "[$(date)] ERROR: missing per-chrom $f"; rm -f "$TMP_OUT"; exit 1; }
    if [ $first -eq 1 ]; then cat "$f" >> "$TMP_OUT"; first=0
    else tail -n +2 "$f" >> "$TMP_OUT"; fi
done
mv "$TMP_OUT" "$FINAL_OUT"
[ -n "${KEEP_DB:-}" ] || rm -f "${DB_DIR}/${SAMPLE}.jf"

echo "[$(date)] DONE $SAMPLE -> $FINAL_OUT"
ls -lh "$FINAL_OUT"; wc -l "$FINAL_OUT"

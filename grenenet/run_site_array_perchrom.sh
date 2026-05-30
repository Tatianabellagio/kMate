#!/bin/bash
#SBATCH --job-name=kmate_grenenet
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --output=logs/kmate_grenenet_%A_%a.out
#SBATCH --error=logs/kmate_grenenet_%A_%a.err

# GrENE-Net production scale-out — run kMate across many pool-seq samples in
# parallel (one SLURM array task per sample of the ~2,415 evolved GrENE-Net
# libraries). kMate the *algorithm* lives in src/; this is the GrENE-Net
# *application*. The 64 GB allocation (vs 200 GB) fits memex nodes, unlocking
# ~3x more cluster concurrency.
#
# STATUS (2026-05-30): updated to the current production recipe but NOT yet run
# — executing the full GrENE-Net cohort is a future session's task. Before
# running, confirm the matrices below (see notes) and that Chr2-5 are built.
#
# Usage:
#   sbatch --array=1-N%C --export=MANIFEST=...,OUT_DIR=... grenenet/run_site_array_perchrom.sh
#
# Required env:
#   MANIFEST     — TSV with header + columns: sample_id, reads_path[, reads_path2]
#   OUT_DIR      — output dir (relative to project root or absolute)
# Optional env (defaults = current production recipe):
#   BLOCK_MODE   — "window" (evolved samples, default) | "global" (F0 / SEEDMIX pools)
#   WINDOW_BP    — window-mode width, default 10000 (the production window recipe)
#   KMER_PA_PREFIX, VAR_PA, VAR_CALLED, VAR_META — panel matrices (see notes below)
#   CHROMS       — quoted space-separated, default "Chr1" (only Chr1 is built)

set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs

: ${MANIFEST:?Set MANIFEST to a TSV path with sample_id, reads_path[, reads_path2]}
: ${OUT_DIR:?Set OUT_DIR to a results subdir}
BLOCK_MODE=${BLOCK_MODE:-window}
WINDOW_BP=${WINDOW_BP:-10000}

# Production matrices (2026-05-30):
#   kmer_pa  — the production filter is filt2inv (filt2 + invariant-column cut;
#              ALGORITHM.md §2.1). That matrix is not yet rebuilt on disk, so the
#              default points at the existing _filt2; switch to
#              kmer_pa_231_v3qc_v3_filt2inv once the rebuild is done.
#   var_pa / var_called / var_meta — arch3 decomposition (Chr1 built; Chr2-5 pending).
KMER_PA_PREFIX=${KMER_PA_PREFIX:-data/kmer_pa_231_v3qc_v3_filt2/kmer_pa}
VAR_PA=${VAR_PA:-panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz}
VAR_CALLED=${VAR_CALLED:-panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz}
VAR_META=${VAR_META:-panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz}
CHROMS=${CHROMS:-"Chr1"}

mkdir -p $OUT_DIR

LINE=$((SLURM_ARRAY_TASK_ID + 1))   # +1 to skip header
ROW=$(awk -F'\t' -v n=$LINE 'NR==n' $MANIFEST)
SAMPLE=$(echo "$ROW" | cut -f1)
R1=$(echo "$ROW" | cut -f2)
R2=$(echo "$ROW" | cut -f3)

OUT=${OUT_DIR}/${SAMPLE}.tsv
if [ -s "$OUT" ]; then
    echo "[$(date)] $SAMPLE already done — skipping"
    exit 0
fi

echo "[$(date)] task=${SLURM_ARRAY_TASK_ID}  sample=${SAMPLE}  block_mode=${BLOCK_MODE}"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  kmer_pa_prefix: $KMER_PA_PREFIX  var_pa: $VAR_PA"

# Build reads args (handle single FASTQ or paired)
if [ -n "${R2:-}" ] && [ -s "$R2" ]; then
    READS_ARGS="--reads $R1 $R2"
else
    READS_ARGS="--reads $R1"
fi

WINDOW_ARGS=""
if [ "$BLOCK_MODE" = "window" ]; then
    WINDOW_ARGS="--window-bp $WINDOW_BP"
fi

/usr/bin/time -v /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --kmer-pa-prefix $KMER_PA_PREFIX \
    --var-pa $VAR_PA \
    --var-called $VAR_CALLED \
    --var-meta $VAR_META \
    $READS_ARGS \
    --sample $SAMPLE \
    --out $OUT \
    --threads 8 \
    --block-mode $BLOCK_MODE \
    --kmer-weight inv_mb \
    --chroms $CHROMS \
    $WINDOW_ARGS

echo "[$(date)] DONE — $OUT"
ls -lh $OUT

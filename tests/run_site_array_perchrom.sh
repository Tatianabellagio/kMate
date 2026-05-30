#!/bin/bash
#SBATCH --job-name=cem_perchrom
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --output=logs/cem_perchrom_%A_%a.out
#SBATCH --error=logs/cem_perchrom_%A_%a.err

# Production cactus_em runner — per-chromosome driver, window-mode EM.
# Replaces run_site04_array.sh. Uses 64 GB allocation (vs 200 GB) so it fits
# on memex nodes, unlocking ~3x more cluster concurrency.
#
# Usage:
#   sbatch --array=1-N%C --export=MANIFEST=...,OUT_DIR=...,BLOCK_MODE=window run_site_array_perchrom.sh
#
# Required env:
#   MANIFEST     — TSV with header + columns: sample_id, reads_path[, reads_path2]
#   OUT_DIR      — output dir (relative to project root or absolute)
# Optional env:
#   BLOCK_MODE   — "window" (evolved samples, default) | "global" (F0 pools)
#   WINDOW_BP    — for window mode, default 200000
#   CN_PREFIX    — default data/kmer_pa_231_v2/kmer_pa
#   CN_VAR       — default data/var_pa_231_v2.var_pa.npz
#   CN_VAR_META  — default data/var_pa_231_v2.meta.npz
#   CHROMS       — quoted space-separated, default "Chr1 Chr2 Chr3 Chr4 Chr5"

set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p tests/logs

: ${MANIFEST:?Set MANIFEST to a TSV path with sample_id, reads_path[, reads_path2]}
: ${OUT_DIR:?Set OUT_DIR to a results subdir}
BLOCK_MODE=${BLOCK_MODE:-window}
WINDOW_BP=${WINDOW_BP:-200000}
CN_PREFIX=${CN_PREFIX:-data/kmer_pa_231_v2/kmer_pa}
CN_VAR=${CN_VAR:-data/var_pa_231_v2.var_pa.npz}
CN_VAR_META=${CN_VAR_META:-data/var_pa_231_v2.meta.npz}
CHROMS=${CHROMS:-"Chr1 Chr2 Chr3 Chr4 Chr5"}

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
echo "  kmer_pa_prefix: $CN_PREFIX  var_pa: $CN_VAR"

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
    --kmer-pa-prefix $CN_PREFIX \
    --var-pa $CN_VAR \
    --var-meta $CN_VAR_META \
    $READS_ARGS \
    --sample $SAMPLE \
    --out $OUT \
    --threads 8 \
    --block-mode $BLOCK_MODE \
    --chroms $CHROMS \
    $WINDOW_ARGS

echo "[$(date)] DONE — $OUT"
ls -lh $OUT

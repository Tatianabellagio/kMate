#!/bin/bash
#SBATCH --job-name=kmate_phaseA_db
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --requeue
#SBATCH --output=logs/kmate_phaseA_%A_%a.out
#SBATCH --error=logs/kmate_phaseA_%A_%a.err

# ============================================================================
# TWO-PHASE COHORT RUNNER — PHASE A: build one k-mer DB per sample.
#
# ⚠️  OPTIONAL / NOT YET RUN AT SCALE — single-phase run_site_array_perchrom.sh
# is the validated production path. See grenenet/README.md and the note in
# run_phaseB_em_per_chrom.sh before choosing two-phase.
#
# This is the "count once" step factored OUT of the per-chrom EM so it can be
# right-sized independently: jellyfish count is CPU-bound and scales with
# threads (8), but needs little RAM (~the 3G hash). Phase B (the per-chrom EM)
# is low-CPU / high-RAM and packs many more tasks per node — keeping the two
# phases separate lets each run at its efficient footprint.
#
# One array task per sample. Builds ${DB_DIR}/${SAMPLE}.jf (canonical, k=31),
# skipped if already present (preemption/requeue-safe). Phase C removes each DB
# after that sample's genome-wide TSV is concatenated.
#
# Usage:
#   sbatch --array=1-N%C \
#     --export=ALL,MANIFEST=...,DB_DIR=... \
#     grenenet/run_phaseA_build_dbs.sh
#
# Required env:
#   MANIFEST  — TSV: header + columns sample_id, reads_path[, reads_path2]
#   DB_DIR    — where to write ${SAMPLE}.jf (needs ~5 GB/sample free)
# Optional env:
#   K         — k-mer length (default 31)
#   HASH_SIZE — jellyfish initial hash (default 3G)
#   THREADS   — jellyfish threads (default 8)
# ============================================================================
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs

PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
JELLYFISH=$(dirname "$PY")/jellyfish
SAMTOOLS=$(dirname "$PY")/samtools

: ${MANIFEST:?Set MANIFEST}
: ${DB_DIR:?Set DB_DIR (where ${SAMPLE}.jf is written)}
K=${K:-31}
HASH_SIZE=${HASH_SIZE:-3G}
THREADS=${THREADS:-8}
mkdir -p "$DB_DIR"

# OFFSET lets a wave address manifest rows past the MaxArraySize=1001 cap:
# row = OFFSET + TASK_ID (+1 for the header). Default 0.
OFFSET=${OFFSET:-0}
LINE=$((SLURM_ARRAY_TASK_ID + OFFSET + 1))   # +1 to skip header
ROW=$(awk -F'\t' -v n=$LINE 'NR==n' "$MANIFEST")
SAMPLE=$(echo "$ROW" | cut -f1)
R1=$(echo "$ROW" | cut -f2)
R2=$(echo "$ROW" | cut -f3)

JF_DB=${DB_DIR}/${SAMPLE}.jf
if [ -s "$JF_DB" ]; then
    echo "[$(date)] $SAMPLE DB already present ($JF_DB) — skipping"
    exit 0
fi

echo "[$(date)] PhaseA build DB  sample=$SAMPLE  -> $JF_DB"
echo "  R1: $R1"; echo "  R2: $R2"
if [ -n "${R2:-}" ] && [ -s "$R2" ]; then SRC="$R1 $R2"; else SRC="$R1"; fi

if [[ "$R1" == *.bam ]]; then
    $JELLYFISH count -m $K -s $HASH_SIZE -t $THREADS -C -o "$JF_DB" \
        <($SAMTOOLS fastq -F 0x900 "$R1" 2>/dev/null) \
      || { echo "[$(date)] ERROR: jellyfish count failed ($SAMPLE)"; rm -f "$JF_DB"; exit 1; }
else
    bash -c "set -o pipefail; zcat -f $SRC | $JELLYFISH count -m $K -s $HASH_SIZE -t $THREADS -C -o '$JF_DB' /dev/fd/0" \
      || { echo "[$(date)] ERROR: jellyfish count failed ($SAMPLE)"; rm -f "$JF_DB"; exit 1; }
fi
echo "[$(date)] DONE $SAMPLE DB: $(ls -lh "$JF_DB" | awk '{print $5}')"

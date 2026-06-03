#!/bin/bash
#SBATCH --job-name=kmate_grenenet
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --requeue
#SBATCH --output=logs/kmate_grenenet_%A_%a.out
#SBATCH --error=logs/kmate_grenenet_%A_%a.err

# GrENE-Net production scale-out — run kMate across many pool-seq samples in
# parallel (one SLURM array task per sample of the ~2,415 evolved GrENE-Net
# libraries). kMate the *algorithm* lives in src/; this is the GrENE-Net
# *application*. This is the VALIDATED production launcher (pilot: sites 4+54,
# 175 samples, 2026-06-01). See grenenet/README.md for the full scale-out guide.
#
# Each array task processes ONE sample across ALL chromosomes by calling the
# driver once per chrom (the driver dense-loads one chrom's K_pa at a time, and
# our arch3 V_pa lives in per-chrom files), writing ${SAMPLE}_${CHR}.tsv per
# chrom, then concatenating into a genome-wide ${SAMPLE}.tsv. The per-chrom TSVs
# and the final TSV are skipped if already present, so a preempted (lowprio +
# --requeue) task resumes at the chrom it left off.
#
# COUNT ONCE: k-mer counting (scanning the full read pool) dominated each chrom's
# wall and was otherwise repeated identically for all 5 chroms. So we build ONE
# canonical Jellyfish DB per sample up front and each per-chrom driver call
# QUERIES it (--kmer-db) instead of re-scanning the reads (~1.7-2x, byte-
# identical). The DB goes on fast node-local storage (JF_DIR, default auto ->
# /dev/shm) so concurrent tasks don't thrash the shared filesystem reading it.
# It is REBUILT cheaply (~40s) on requeue and removed on every exit (incl.
# preemption) by the trap — the per-chrom TSVs are what actually survive.
#
# Production recipe: env=kmate; K_pa=arch3 filt2inv; V_pa=arch3 per-chrom triplet;
# EM weight inv_mb; GLOBAL mode (the evolved-cohort choice; window only on request).
#
# Usage (chunk by OFFSET if N>1000; see grenenet/README.md):
#   sbatch --array=1-N%C --mem=32G \
#     --export=ALL,MANIFEST=...,OUT_DIR=...,BLOCK_MODE=global \
#     grenenet/run_site_array_perchrom.sh
#
# Required env:
#   MANIFEST     — TSV with header + columns: sample_id, reads_path[, reads_path2]
#   OUT_DIR      — output dir (relative to project root or absolute)
# Optional env (defaults = current production recipe):
#   BLOCK_MODE   — "global" (chrom-wide EM, DEFAULT) | "window" (per-10kb block-EM;
#                  only when explicitly requested)
#   JF_DIR       — where the transient k-mer DB lives (default auto: $SLURM_TMPDIR
#                  -> /dev/shm -> OUT_DIR). The key scale-out knob (see below).
#   OFFSET       — manifest-row offset for chunking past MaxArraySize=1001 (default 0)
#   WINDOW_BP    — window-mode width, default 10000 (only used when BLOCK_MODE=window)
#   KMER_WEIGHT  — "inv_mb" (production, default) | "uniform"
#   CHROMS       — quoted space-separated, default "Chr1 Chr2 Chr3 Chr4 Chr5"
#   KMER_PA_PREFIX — K_pa prefix; driver appends _<CHR>.{kmer_pa,meta}.npz
#   VAR_PA_DIR / VAR_PA_TAG — per-chrom V_pa location; paths are built as
#                  $VAR_PA_DIR/<chr>/${VAR_PA_TAG}_<chr>.{var_pa,var_called,meta}.npz

set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs

# --- whole-job runtime reporting (fires on ANY exit, incl. error/preemption) ---
_T0=$SECONDS
_report_runtime() {
    local rc=$? s=$((SECONDS - _T0))
    # Always free the k-mer DB — critical when JF_DIR=/dev/shm (RAM), else a
    # preempted/failed task leaks gigabytes of node RAM.
    rm -f "${JF_DB:-}"
    printf '[%s] RUNTIME %dh%02dm%02ds (%ds total wall)  sample=%s  exit=%d\n' \
        "$(date)" $((s/3600)) $(((s%3600)/60)) $((s%60)) "$s" "${SAMPLE:-?}" "$rc"
}
trap _report_runtime EXIT

PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python

: ${MANIFEST:?Set MANIFEST to a TSV path with sample_id, reads_path[, reads_path2]}
: ${OUT_DIR:?Set OUT_DIR to a results subdir}
BLOCK_MODE=${BLOCK_MODE:-global}
WINDOW_BP=${WINDOW_BP:-10000}
KMER_WEIGHT=${KMER_WEIGHT:-inv_mb}
CHROMS=${CHROMS:-"Chr1 Chr2 Chr3 Chr4 Chr5"}

# Production matrices (2026-05-31, arch3 single-source — docs/PIPELINE_STATE.md §0):
#   K_pa : in-house index + merged_231 VCF, filt2inv, N-on  (per-chrom under one prefix)
#   V_pa : arch3 decomposition, per-chrom triplet (var_pa / var_called / meta)
KMER_PA_PREFIX=${KMER_PA_PREFIX:-data/kmer_pa_231_arch3_filt2inv/kmer_pa}
VAR_PA_DIR=${VAR_PA_DIR:-panel/arch3}
VAR_PA_TAG=${VAR_PA_TAG:-var_pa_231_arch3}

# Where to put the transient per-sample k-mer DB. This is a SCALE-OUT tuning
# knob, not part of the kMate method: the DB is unique + cold per sample, so on a
# congested shared filesystem many concurrent tasks reading their 2-3 GB DBs is
# the dominant stall ("Vessel B"; measured up to 48x query slowdown at 75-way
# concurrency). Putting it on fast node-local storage removes it from the shared
# FS entirely (query tail 32 min -> 39 s in an A/B test).
#
# Auto-pick the fastest SAFE location (override by setting JF_DIR explicitly):
#   1. node-local SSD if the scheduler exposes one ($SLURM_TMPDIR)
#   2. else /dev/shm (RAM) when it has >4 GB free  [needs --mem to cover ~25 GB
#      peak + ~3 GB DB; the exit trap frees it even on preemption]
#   3. else OUT_DIR (shared scratch) — always works, just slower under load
if [ -z "${JF_DIR:-}" ]; then
    if [ -n "${SLURM_TMPDIR:-}" ] && [ -d "${SLURM_TMPDIR}" ]; then
        JF_DIR=$SLURM_TMPDIR
    elif [ -d /dev/shm ] && [ -w /dev/shm ] && \
         [ "$(df -Pk /dev/shm 2>/dev/null | awk 'NR==2{print $4+0}')" -gt 4194304 ]; then
        JF_DIR=/dev/shm
    else
        JF_DIR=$OUT_DIR
    fi
fi
echo "[$(date)] k-mer DB location (JF_DIR): $JF_DIR"

mkdir -p $OUT_DIR

# OFFSET lets a chunk address manifest rows beyond the cluster MaxArraySize
# (1001) cap: submit the cohort as chunks of <=1000 with OFFSET=0,1000,2000...
# and the manifest row = OFFSET + TASK_ID (+1 for the header). Default 0.
OFFSET=${OFFSET:-0}
LINE=$((SLURM_ARRAY_TASK_ID + OFFSET + 1))   # +1 to skip header
ROW_RAW=$(awk -F'\t' -v n=$LINE 'NR==n' $MANIFEST)
if [ -z "$ROW_RAW" ]; then
    echo "[$(date)] no manifest row at line $LINE (OFFSET=$OFFSET task=$SLURM_ARRAY_TASK_ID) — past end, nothing to do"
    exit 0
fi
ROW="$ROW_RAW"
SAMPLE=$(echo "$ROW" | cut -f1)
R1=$(echo "$ROW" | cut -f2)
R2=$(echo "$ROW" | cut -f3)

FINAL_OUT=${OUT_DIR}/${SAMPLE}.tsv
if [ -s "$FINAL_OUT" ]; then
    echo "[$(date)] $SAMPLE already complete ($FINAL_OUT) — skipping"
    exit 0
fi

echo "[$(date)] task=${SLURM_ARRAY_TASK_ID}  sample=${SAMPLE}  block_mode=${BLOCK_MODE}  weight=${KMER_WEIGHT}"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  K_pa prefix: $KMER_PA_PREFIX   V_pa: $VAR_PA_DIR/<chr>/${VAR_PA_TAG}_<chr>"
echo "  chroms: $CHROMS"

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

# ---- build the k-mer DB ONCE (count the read pool a single time) ----
# Each per-chrom driver call queries this DB (--kmer-db) instead of re-scanning
# the reads, removing the ~5x redundant counting. Reused within this task if
# present; rebuilt on requeue (it lives in JF_DIR, typically RAM, so it doesn't
# persist across nodes). Removed on every exit by the trap + after the concat.
JELLYFISH=$(dirname "$PY")/jellyfish
mkdir -p "$JF_DIR"
# Task-unique name so concurrent tasks sharing a node's /dev/shm never collide.
JF_DB=${JF_DIR}/${SAMPLE}.${SLURM_ARRAY_JOB_ID:-x}_${SLURM_ARRAY_TASK_ID:-0}.jf
if [ ! -s "$JF_DB" ]; then
    echo "[$(date)] building k-mer DB once -> $JF_DB"
    if [ -n "${R2:-}" ] && [ -s "$R2" ]; then SRC="$R1 $R2"; else SRC="$R1"; fi
    if [[ "$R1" == *.bam ]]; then
        $JELLYFISH count -m 31 -s 3G -t 8 -C -o "$JF_DB" \
            <(/global/home/users/tbellg/miniforge3/envs/kmate/bin/samtools fastq -F 0x900 "$R1" 2>/dev/null) \
          || { echo "[$(date)] ERROR: jellyfish count failed for $SAMPLE"; rm -f "$JF_DB"; exit 1; }
    else
        bash -c "set -o pipefail; zcat -f $SRC | $JELLYFISH count -m 31 -s 3G -t 8 -C -o '$JF_DB' /dev/fd/0" \
          || { echo "[$(date)] ERROR: jellyfish count failed for $SAMPLE"; rm -f "$JF_DB"; exit 1; }
    fi
    echo "[$(date)] k-mer DB built: $(ls -lh "$JF_DB" | awk '{print $5}')"
fi

# ---- per-chrom runs ----
PER_CHROM_OUTS=()
for CHR in $CHROMS; do
    chrlc=${CHR,,}
    VAR_PA=${VAR_PA_DIR}/${chrlc}/${VAR_PA_TAG}_${chrlc}.var_pa.npz
    VAR_CALLED=${VAR_PA_DIR}/${chrlc}/${VAR_PA_TAG}_${chrlc}.var_called.npz
    VAR_META=${VAR_PA_DIR}/${chrlc}/${VAR_PA_TAG}_${chrlc}.meta.npz
    OUT_CHR=${OUT_DIR}/${SAMPLE}_${CHR}.tsv

    for f in "$VAR_PA" "$VAR_CALLED" "$VAR_META" "${KMER_PA_PREFIX}_${CHR}.kmer_pa.npz"; do
        [ -s "$f" ] || { echo "[$(date)] ERROR: missing matrix $f — aborting $SAMPLE"; exit 1; }
    done

    PER_CHROM_OUTS+=("$OUT_CHR")
    if [ -s "$OUT_CHR" ]; then
        echo "[$(date)]   $CHR already done ($OUT_CHR) — skipping"
        continue
    fi

    echo "[$(date)]   --- $SAMPLE $CHR ---"
    /usr/bin/time -v $PY -u src/per_sample_per_chrom.py \
        --kmer-pa-prefix $KMER_PA_PREFIX \
        --var-pa $VAR_PA \
        --var-called $VAR_CALLED \
        --var-meta $VAR_META \
        $READS_ARGS \
        --kmer-db $JF_DB \
        --sample $SAMPLE \
        --out $OUT_CHR \
        --threads 8 \
        --block-mode $BLOCK_MODE \
        --kmer-weight $KMER_WEIGHT \
        --chroms $CHR \
        $WINDOW_ARGS \
      || { echo "[$(date)] ERROR: driver failed on $SAMPLE $CHR"; exit 1; }
done

# ---- concatenate per-chrom TSVs into the genome-wide per-sample TSV ----
# Keep header from the first file, append record bodies from each.
echo "[$(date)] concatenating ${#PER_CHROM_OUTS[@]} per-chrom TSVs -> $FINAL_OUT"
first=1
TMP_OUT=${FINAL_OUT}.tmp
: > "$TMP_OUT"
for f in "${PER_CHROM_OUTS[@]}"; do
    [ -s "$f" ] || { echo "[$(date)] ERROR: expected per-chrom output $f missing"; rm -f "$TMP_OUT"; exit 1; }
    if [ $first -eq 1 ]; then
        cat "$f" >> "$TMP_OUT"; first=0
    else
        tail -n +2 "$f" >> "$TMP_OUT"   # skip the 8-col header on subsequent chroms
    fi
done
mv "$TMP_OUT" "$FINAL_OUT"

# Genome-wide TSV is complete; the shared k-mer DB is no longer needed.
rm -f "$JF_DB"

echo "[$(date)] DONE — $FINAL_OUT"
ls -lh "$FINAL_OUT"
wc -l "$FINAL_OUT"

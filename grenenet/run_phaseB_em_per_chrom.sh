#!/bin/bash
#SBATCH --job-name=kmate_phaseB_em
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=logs/kmate_phaseB_%A_%a.out
#SBATCH --error=logs/kmate_phaseB_%A_%a.err

# ============================================================================
# TWO-PHASE COHORT RUNNER — PHASE B: one (sample, chrom) EM per array task.
#
# ⚠️  OPTIONAL / NOT YET RUN AT SCALE — single-phase run_site_array_perchrom.sh
# is the validated path. Phase B reads each sample's DB from DB_DIR (shared
# scratch, built by Phase A on a possibly different node), so it canNOT put the
# DB node-local and IS subject to the "Vessel-B" shared-FS stall the single-phase
# runner avoids. Prefer single-phase. See grenenet/README.md.
#
# The genuinely parallel unit: a sample's 5 chroms are independent once its DB
# exists (Phase A). So we expose ALL N_samples x 5 chrom solves as one flat
# array — maximum fan-out, perfectly packed. Each task QUERIES the prebuilt DB
# (no read re-scan) for its chrom, runs the EM, and writes ${SAMPLE}_${CHR}.tsv.
#
# Footprint is right-sized for this phase: the per-chrom work is mostly
# single-threaded (query + dense-load + write) with a BLAS-threaded EM, so 4
# cores suffice (the old 8-core task averaged only ~2.7 cores here). Measured
# peak RAM ~21 GB (dense kmer_pa on Chr1) -> 32 GB.
#
# ARRAY INDEXING (no 12k-row manifest needed): task T (1-based) maps to
#   sample_idx = (T-1) / 5   (0-based row into the sample manifest)
#   chrom_idx  = (T-1) % 5   (into CHROMS)
# So --array=1-(N*5). Mind the cluster MaxArraySize (see submit helper for
# chunking if N*5 exceeds it).
#
# Usage:
#   sbatch --array=1-$((N*5))%C \
#     --export=ALL,MANIFEST=...,DB_DIR=...,OUT_DIR=... \
#     grenenet/run_phaseB_em_per_chrom.sh
#
# Required env: MANIFEST, DB_DIR (Phase A output), OUT_DIR
# Optional env (defaults = production recipe; mirror run_site_array_perchrom.sh):
#   BLOCK_MODE (global, DEFAULT; window only if explicitly set) WINDOW_BP (10000) KMER_WEIGHT (inv_mb)
#   CHROMS ("Chr1 Chr2 Chr3 Chr4 Chr5")
#   KMER_PA_PREFIX VAR_PA_DIR VAR_PA_TAG
# ============================================================================
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs

PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python

: ${MANIFEST:?Set MANIFEST}
: ${DB_DIR:?Set DB_DIR (Phase A output)}
: ${OUT_DIR:?Set OUT_DIR}
BLOCK_MODE=${BLOCK_MODE:-global}
WINDOW_BP=${WINDOW_BP:-10000}
KMER_WEIGHT=${KMER_WEIGHT:-inv_mb}
CHROMS=${CHROMS:-"Chr1 Chr2 Chr3 Chr4 Chr5"}
KMER_PA_PREFIX=${KMER_PA_PREFIX:-data/kmer_pa_231_arch3_filt2inv/kmer_pa}
VAR_PA_DIR=${VAR_PA_DIR:-panel/arch3}
VAR_PA_TAG=${VAR_PA_TAG:-var_pa_231_arch3}
mkdir -p "$OUT_DIR"

# OFFSET lets a wave address (sample,chrom) units past the MaxArraySize=1001
# cap: global unit G = OFFSET + TASK_ID, then split into sample/chrom. Default 0.
OFFSET=${OFFSET:-0}
GUNIT=$(( SLURM_ARRAY_TASK_ID + OFFSET ))            # 1-based global (sample,chrom) unit
CHROMS_ARR=($CHROMS)
NCHR=${#CHROMS_ARR[@]}
SAMPLE_IDX=$(( (GUNIT - 1) / NCHR ))                 # 0-based sample row
CHR_IDX=$(( (GUNIT - 1) % NCHR ))
CHR=${CHROMS_ARR[$CHR_IDX]}

LINE=$(( SAMPLE_IDX + 2 ))                           # +2: 1-based + header
ROW=$(awk -F'\t' -v n=$LINE 'NR==n' "$MANIFEST")
SAMPLE=$(echo "$ROW" | cut -f1)
[ -n "$SAMPLE" ] || { echo "[$(date)] ERROR: no sample at manifest line $LINE (task $SLURM_ARRAY_TASK_ID)"; exit 1; }

OUT_CHR=${OUT_DIR}/${SAMPLE}_${CHR}.tsv
if [ -s "$OUT_CHR" ]; then
    echo "[$(date)] $SAMPLE $CHR already done ($OUT_CHR) — skipping"; exit 0
fi

JF_DB=${DB_DIR}/${SAMPLE}.jf
[ -s "$JF_DB" ] || { echo "[$(date)] ERROR: DB missing $JF_DB — run Phase A first"; exit 1; }

chrlc=${CHR,,}
VAR_PA=${VAR_PA_DIR}/${chrlc}/${VAR_PA_TAG}_${chrlc}.var_pa.npz
VAR_CALLED=${VAR_PA_DIR}/${chrlc}/${VAR_PA_TAG}_${chrlc}.var_called.npz
VAR_META=${VAR_PA_DIR}/${chrlc}/${VAR_PA_TAG}_${chrlc}.meta.npz
for f in "$VAR_PA" "$VAR_CALLED" "$VAR_META" "${KMER_PA_PREFIX}_${CHR}.kmer_pa.npz"; do
    [ -s "$f" ] || { echo "[$(date)] ERROR: missing matrix $f"; exit 1; }
done

WINDOW_ARGS=""
[ "$BLOCK_MODE" = "window" ] && WINDOW_ARGS="--window-bp $WINDOW_BP"

echo "[$(date)] PhaseB task=$SLURM_ARRAY_TASK_ID (gunit=$GUNIT)  sample=$SAMPLE  $CHR  mode=$BLOCK_MODE"
/usr/bin/time -v $PY -u src/per_sample_per_chrom.py \
    --kmer-pa-prefix $KMER_PA_PREFIX \
    --var-pa $VAR_PA --var-called $VAR_CALLED --var-meta $VAR_META \
    --reads /dev/null \
    --kmer-db $JF_DB \
    --sample $SAMPLE --out $OUT_CHR \
    --threads ${SLURM_CPUS_PER_TASK:-4} \
    --block-mode $BLOCK_MODE --kmer-weight $KMER_WEIGHT \
    --chroms $CHR $WINDOW_ARGS \
  || { echo "[$(date)] ERROR: driver failed $SAMPLE $CHR"; rm -f "$OUT_CHR"; exit 1; }
echo "[$(date)] DONE $SAMPLE $CHR -> $OUT_CHR"

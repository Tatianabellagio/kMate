#!/bin/bash
#SBATCH --job-name=moi_kmate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=benchmarks/ruth_outcross/logs/moi_%A_%a.out
#SBATCH --error=benchmarks/ruth_outcross/logs/moi_%A_%a.err

# ============================================================================
# MOI / outcrossing benchmark — run kMate (global) on individually-sequenced
# greenhouse plants from GrENE-net seeds, to recover the founder-mixture h and
# call inbred (h ~ one-hot, eff_n_founders ~ 1) vs outcrossed (h spread over
# >=2 founders, eff_n_founders >= 2). Compares kMate to hapFIRE's call.
#
# One array task per sample (1-10). Each task:
#   1. builds a node-local k-mer DB ONCE from the sample's deduped BAM
#      (jellyfish reads the BAM via `samtools fastq`), then
#   2. runs the per-chrom global EM for all 5 chroms in --h-only mode (skips the
#      AF projection + var_pa load — we only need the founder mixture h for the
#      MOI call), writing just ${SAMPLE}_${CHR}.h_per_chrom.npz.
#
# Usage:  sbatch --array=1-10 benchmarks/ruth_outcross/run_moi.sbatch
# ============================================================================
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
B=benchmarks/ruth_outcross
mkdir -p "$B/logs" "$B/out"

PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
JELLYFISH=$(dirname "$PY")/jellyfish
SAMTOOLS=$(dirname "$PY")/samtools

MANIFEST=$B/manifest.tsv
OUT_DIR=$B/out
KMER_PA_PREFIX=data/kmer_pa_231_arch3_filt2inv/kmer_pa
CHROMS="Chr1 Chr2 Chr3 Chr4 Chr5"
K=31; HASH_SIZE=3G; THREADS=${SLURM_CPUS_PER_TASK:-8}

LINE=$((SLURM_ARRAY_TASK_ID + 1))   # +1 skip header
ROW=$(awk -F'\t' -v n=$LINE 'NR==n' "$MANIFEST")
SAMPLE=$(echo "$ROW" | cut -f1)
BAM=$(echo "$ROW" | cut -f2)
LABEL=$(echo "$ROW" | cut -f3)
[ -n "$SAMPLE" ] || { echo "ERROR: no sample at line $LINE"; exit 1; }
[ -s "$BAM" ]    || { echo "ERROR: BAM missing $BAM"; exit 1; }

# node-local k-mer DB (RAM if available, else node scratch) — avoids the
# shared-FS query stall; removed at task end.
JF_BASE=${TMPDIR:-/tmp}; [ -d /dev/shm ] && JF_BASE=/dev/shm
JF_DB=$JF_BASE/${SAMPLE}.jf
trap 'rm -f "$JF_DB"' EXIT

echo "[$(date)] MOI sample=$SAMPLE label=$LABEL"
echo "  BAM: $BAM"
echo "  DB : $JF_DB"

if [ ! -s "$JF_DB" ]; then
    $JELLYFISH count -m $K -s $HASH_SIZE -t $THREADS -C -o "$JF_DB" \
        <($SAMTOOLS fastq -F 0x900 "$BAM" 2>/dev/null) \
      || { echo "ERROR: jellyfish count failed ($SAMPLE)"; exit 1; }
fi
echo "[$(date)] DB built: $(ls -lh "$JF_DB" | awk '{print $5}')"

for CHR in $CHROMS; do
    OUT_CHR=$OUT_DIR/${SAMPLE}_${CHR}.tsv
    if [ -s "${OUT_CHR%.tsv}.h_per_chrom.npz" ]; then
        echo "[$(date)] $SAMPLE $CHR already done — skipping"; continue
    fi
    echo "[$(date)] EM $SAMPLE $CHR (global, h-only)"
    $PY -u src/per_sample_per_chrom.py \
        --kmer-pa-prefix $KMER_PA_PREFIX \
        --reads /dev/null --kmer-db "$JF_DB" \
        --sample $SAMPLE --out $OUT_CHR \
        --threads $THREADS --block-mode global --kmer-weight inv_mb \
        --chroms $CHR --h-only \
      || { echo "ERROR: EM failed $SAMPLE $CHR"; exit 1; }
done
echo "[$(date)] DONE $SAMPLE"

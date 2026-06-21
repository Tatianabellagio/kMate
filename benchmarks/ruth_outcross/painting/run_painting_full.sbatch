#!/bin/bash
#SBATCH --job-name=paint_full
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --output=benchmarks/ruth_outcross/painting/logs/paintfull_%A_%a.out
#SBATCH --error=benchmarks/ruth_outcross/painting/logs/paintfull_%A_%a.err

# Full 5-chromosome chromosome painting (10 kb windows + HMM smoothing) on a
# spread of hapFIRE "extreme" individuals (inbred / F1 / recombinant / complex).
# One array task per sample: build the k-mer DB once from the BAM, then ONE
# window-mode call over all 5 chroms -> <SAMPLE>_full_win10k.h_blocks_per_chrom.npz
# (holds per-chrom h_blocks for every chromosome).
#
# Usage: sbatch --array=1-6 benchmarks/ruth_outcross/painting/run_painting_full.sbatch
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
B=benchmarks/ruth_outcross/painting
mkdir -p "$B/logs" "$B/out"
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
JELLYFISH=$(dirname "$PY")/jellyfish; SAMTOOLS=$(dirname "$PY")/samtools

MANIFEST=$B/extremes_samples.tsv
OUT=$B/out
KMER_PA_PREFIX=data/kmer_pa_231_arch3_filt2inv/kmer_pa
VAR_PA_DIR=panel/arch3; VAR_PA_TAG=var_pa_231_arch3
CHROMS="Chr1 Chr2 Chr3 Chr4 Chr5"
THREADS=${SLURM_CPUS_PER_TASK:-8}

LINE=$((SLURM_ARRAY_TASK_ID + 1))
ROW=$(awk -F'\t' -v n=$LINE 'NR==n' "$MANIFEST")
SAMPLE=$(echo "$ROW" | cut -f1); BAM=$(echo "$ROW" | cut -f2); CAT=$(echo "$ROW" | cut -f3)
[ -s "$BAM" ] || { echo "ERROR: BAM missing $BAM"; exit 1; }

OUT_TSV=$OUT/${SAMPLE}_full_win10k.tsv
if [ -s "${OUT_TSV%.tsv}.h_blocks_per_chrom.npz" ]; then echo "done, skip"; exit 0; fi

JF_BASE=${TMPDIR:-/tmp}; [ -d /dev/shm ] && JF_BASE=/dev/shm
JF_DB=$JF_BASE/${SAMPLE}.paintfull.jf
trap 'rm -f "$JF_DB"' EXIT

echo "[$(date)] PAINT-FULL $SAMPLE ($CAT) all chroms"
$JELLYFISH count -m 31 -s 3G -t $THREADS -C -o "$JF_DB" \
    <($SAMTOOLS fastq -F 0x900 "$BAM" 2>/dev/null) \
  || { echo "ERROR: jellyfish count failed"; exit 1; }
echo "[$(date)] DB built"

# build per-chrom var_pa arg lists is per-chrom; the driver loops chroms internally
# but var-pa is one file per chrom -> we pass Chr1's and the driver auto-resolves? No:
# the driver takes ONE var_pa for all chroms. var_pa here is per-chrom, so run per chrom.
for CHR in $CHROMS; do
    chrlc=${CHR,,}
    OUT_CHR=$OUT/${SAMPLE}_${CHR}_full_win10k.tsv
    [ -s "${OUT_CHR%.tsv}.h_blocks_per_chrom.npz" ] && { echo "$CHR done, skip"; continue; }
    VAR_PA=$VAR_PA_DIR/$chrlc/${VAR_PA_TAG}_${chrlc}.var_pa.npz
    VAR_CALLED=$VAR_PA_DIR/$chrlc/${VAR_PA_TAG}_${chrlc}.var_called.npz
    VAR_META=$VAR_PA_DIR/$chrlc/${VAR_PA_TAG}_${chrlc}.meta.npz
    echo "[$(date)] window 10kb $SAMPLE $CHR"
    $PY -u src/per_sample_per_chrom.py \
        --kmer-pa-prefix $KMER_PA_PREFIX \
        --var-pa $VAR_PA --var-called $VAR_CALLED --var-meta $VAR_META \
        --reads /dev/null --kmer-db "$JF_DB" \
        --sample $SAMPLE --out "$OUT_CHR" --threads $THREADS \
        --block-mode window --kmer-weight inv_mb --window-bp 10000 --chroms $CHR \
      || { echo "ERROR: window EM failed $CHR"; rm -f "$OUT_CHR"; exit 1; }
done
echo "[$(date)] DONE $SAMPLE"

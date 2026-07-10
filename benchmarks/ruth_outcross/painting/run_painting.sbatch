#!/bin/bash
#SBATCH --job-name=moi_paint
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=3:00:00
#SBATCH --requeue
#SBATCH --output=benchmarks/ruth_outcross/painting/logs/paint_%A_%a.out
#SBATCH --error=benchmarks/ruth_outcross/painting/logs/paint_%A_%a.err

# ============================================================================
# CHROMOSOME-PAINTING PILOT — per-window founder mixture (local ancestry) along
# the chromosome, to see WHERE recombination happens in each individual.
#
# kMate window mode solves a local h per window; for a single individual that is
# the local ancestry (which founder(s) dominate locally + dosage). Switches along
# the chromosome = recombination breakpoints. Inbred -> one founder flat.
#
# Per sample (array task) we run TWO window definitions on Chr1 to compare:
#   1. fixed 10 kb windows + HMM smoothing (the "star2" recipe)         -> _win10k
#   2. dynld_K500 LD-blocks (analysis/grenenet_gea/blocks_mcf90/)        -> _dynld
# Each writes <SAMPLE>_Chr1_<tag>.h_blocks_per_chrom.npz (the per-window h).
#
# Usage: sbatch --array=1-3 benchmarks/ruth_outcross/painting/run_painting.sbatch
# ============================================================================
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
B=benchmarks/ruth_outcross/painting
mkdir -p "$B/logs" "$B/out"

PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
JELLYFISH=$(dirname "$PY")/jellyfish
SAMTOOLS=$(dirname "$PY")/samtools

MANIFEST=$B/samples.tsv
OUT=$B/out
KMER_PA_PREFIX=data/kmer_pa_231_arch3_filt2inv/kmer_pa
VAR_PA_DIR=panel/arch3; VAR_PA_TAG=var_pa_231_arch3
CHR=Chr1; chrlc=chr1
BLOCKS=analysis/grenenet_gea/blocks_mcf90/${chrlc}_units_dynld_K500.tsv
THREADS=${SLURM_CPUS_PER_TASK:-8}

LINE=$((SLURM_ARRAY_TASK_ID + 1))
ROW=$(awk -F'\t' -v n=$LINE 'NR==n' "$MANIFEST")
SAMPLE=$(echo "$ROW" | cut -f1); BAM=$(echo "$ROW" | cut -f2); LABEL=$(echo "$ROW" | cut -f3)
[ -s "$BAM" ] || { echo "ERROR: BAM missing $BAM"; exit 1; }

JF_BASE=${TMPDIR:-/tmp}; [ -d /dev/shm ] && JF_BASE=/dev/shm
JF_DB=$JF_BASE/${SAMPLE}.paint.jf
trap 'rm -f "$JF_DB"' EXIT

echo "[$(date)] PAINT $SAMPLE ($LABEL) $CHR"
if [ ! -s "$JF_DB" ]; then
    $JELLYFISH count -m 31 -s 3G -t $THREADS -C -o "$JF_DB" \
        <($SAMTOOLS fastq -F 0x900 "$BAM" 2>/dev/null) \
      || { echo "ERROR: jellyfish count failed"; exit 1; }
fi
echo "[$(date)] DB built"

VAR_PA=$VAR_PA_DIR/$chrlc/${VAR_PA_TAG}_${chrlc}.var_pa.npz
VAR_CALLED=$VAR_PA_DIR/$chrlc/${VAR_PA_TAG}_${chrlc}.var_called.npz
VAR_META=$VAR_PA_DIR/$chrlc/${VAR_PA_TAG}_${chrlc}.meta.npz

run_one () {  # $1 = tag, rest = extra args
    local tag=$1; shift
    local o=$OUT/${SAMPLE}_${CHR}_${tag}.tsv
    [ -s "${o%.tsv}.h_blocks_per_chrom.npz" ] && { echo "[$(date)] $tag done, skip"; return; }
    echo "[$(date)] window mode ($tag) $SAMPLE $CHR"
    $PY -u src/per_sample_per_chrom.py \
        --kmer-pa-prefix $KMER_PA_PREFIX \
        --var-pa $VAR_PA --var-called $VAR_CALLED --var-meta $VAR_META \
        --reads /dev/null --kmer-db "$JF_DB" \
        --sample $SAMPLE --out "$o" --threads $THREADS \
        --block-mode window --kmer-weight inv_mb --chroms $CHR "$@" \
      || { echo "ERROR: window EM failed ($tag)"; rm -f "$o"; return 1; }
}

run_one win10k --window-bp 10000
run_one dynld  --blocks-tsv "$BLOCKS"
echo "[$(date)] DONE $SAMPLE"

#!/bin/bash
#SBATCH --job-name=kmer_pa_arch3
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --array=1-5
#SBATCH --output=logs/kmer_pa_arch3_%A_%a.out
#SBATCH --error=logs/kmer_pa_arch3_%A_%a.err
#
# PRODUCTION K_pa (kmer_pa) builder — arch3 single-source. One array task / chrom.
# Builds the founder×k-mer matrix from:
#   index : in-house ours_Chr{N} (default)  |  PanGenie pang_135 index (INDEX=pg, for the
#           Level-B equivalence check)
#   VCF   : panel/arch3/chr{N}/merged_231_chr{N}_final.vcf.gz  (THE production panel; §0)
#   REF   : TAIR10.chr.iupacN.fa
# Applies the production filter inline (--filter-production: keep 2<=ac<=F-1 = filt2inv)
# and --treat-missing-as-n (N-on). Output written already-filtered.
#
# Usage:
#   sbatch scripts/build_kmer_pa_arch3.sh           # INDEX=ours -> data/kmer_pa_231_arch3_filt2inv
#   INDEX=pg sbatch ... scripts/build_kmer_pa_arch3.sh   # -> data/kmer_pa_231_arch3_pgidx_filt2inv
# (env override: `sbatch --export=ALL,INDEX=pg scripts/build_kmer_pa_arch3.sh`)
#
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
$PY -c "import pysam,scipy,numpy" || { echo "ERROR: kmate env lacks pysam/scipy/numpy"; exit 1; }

INDEX=${INDEX:-ours}
N=${SLURM_ARRAY_TASK_ID:-1}
CHR=Chr${N}

case "$INDEX" in
  ours) KMERS=$BASE/panel/pangenie_index/pang_135_haploid/ours_${CHR}_kmers.tsv.gz
        OUT_DIR=$BASE/data/kmer_pa_231_arch3_filt2inv ;;
  pg)   KMERS=$BASE/panel/pangenie_genotyping/data/pang_135_pangenie_index_${CHR}_kmers.tsv.gz
        OUT_DIR=$BASE/data/kmer_pa_231_arch3_pgidx_filt2inv ;;
  *)    echo "ERROR: INDEX must be 'ours' or 'pg' (got '$INDEX')"; exit 1 ;;
esac

VCF=$BASE/panel/arch3/${CHR,,}/merged_231_${CHR,,}_final.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
mkdir -p "$OUT_DIR"
OUT_PREFIX=$OUT_DIR/kmer_pa_${CHR}

[ -s "$KMERS" ] || { echo "ERROR: missing index $KMERS"; exit 1; }
[ -s "$VCF" ]   || { echo "ERROR: missing arch3 VCF $VCF (run arch3 build first)"; exit 1; }
[ -s "$REF" ]   || { echo "ERROR: missing $REF"; exit 1; }
[ ! -s "${OUT_PREFIX}.kmer_pa.npz" ] || { echo "$CHR ($INDEX) exists, skipping"; exit 0; }

echo "[$(date)] $CHR K_pa  index=$INDEX  N-on + filt2inv"
echo "  KMERS=$KMERS"
echo "  VCF=$VCF"
echo "  OUT=$OUT_PREFIX"
$PY -u $BASE/src/build_kmer_pa.py \
    --kmers "$KMERS" --vcf "$VCF" --ref "$REF" \
    --chrom "$CHR" --out "$OUT_PREFIX" \
    --treat-missing-as-n \
    --filter-production --min-ac 2 --invariant-margin 1
echo "[$(date)] $CHR ($INDEX) DONE"
ls -lh ${OUT_PREFIX}.kmer_pa.npz ${OUT_PREFIX}.meta.npz

#!/bin/bash
#SBATCH --job-name=cnfull_idx
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=logs/cnfull_idx_%j.out
#SBATCH --error=logs/cnfull_idx_%j.err
mkdir -p logs
set -euo pipefail

# Build kmer_pa Chr1 from a given kmers.tsv.gz index, against the 231 haploid panel.
# Everything downstream of the index is held IDENTICAL to the production PG-index
# build (build_kmer_pa_v3qc_v3_chr1.sh): same VCF, same ref, --treat-missing-as-n.
# The ONLY thing that varies between runs is the k-mer index -> isolates the
# index builder's effect on kmer_pa (and hence on the EM h estimate).
#
# Usage: sbatch build_kmer_pa_from_index.sh <KMERS_TSV_GZ> <OUT_DIR>

KMERS=$1
OUT_DIR=$2

BASE=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CHR=Chr1
VCF=$BASE/panel/pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa

mkdir -p $OUT_DIR
OUT_PREFIX=$OUT_DIR/cn_${CHR}
[ -s "$KMERS" ] || { echo "ERROR: missing index $KMERS"; exit 1; }
[ -s "$VCF" ]   || { echo "ERROR: missing $VCF"; exit 1; }
[ ! -s "${OUT_PREFIX}.kmer_pa.npz" ] || { echo "exists, skipping"; exit 0; }

echo "[$(date)] kmer_pa from index=$KMERS -> $OUT_PREFIX (--treat-missing-as-n)"
$PY -u $BASE/src/build_kmer_pa.py \
    --kmers "$KMERS" --vcf "$VCF" --ref "$REF" \
    --chrom "$CHR" --out "$OUT_PREFIX" \
    --treat-missing-as-n
echo "[$(date)] DONE"; ls -lh ${OUT_PREFIX}.kmer_pa.npz ${OUT_PREFIX}.meta.npz

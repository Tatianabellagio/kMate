#!/bin/bash
#SBATCH --job-name=kmer_pa_v3qc_v3
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=logs/kmer_pa_v3qc_v3_%j.out
#SBATCH --error=logs/kmer_pa_v3qc_v3_%j.err
# Build kmer_pa for v3qc_v3, Chr1 only first (full genome later if Chr1 looks good).
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CHR=Chr1
KMERS=$BASE/panel/pangenie_genotyping/data/pang_135_pangenie_index_${CHR}_kmers.tsv.gz
VCF=$BASE/panel/pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
OUT_DIR=$BASE/data/kmer_pa_231_v3qc_v3
mkdir -p $OUT_DIR
OUT_PREFIX=$OUT_DIR/cn_${CHR}
[ -s "$VCF" ] || { echo "ERROR: missing $VCF"; exit 1; }
[ ! -s "${OUT_PREFIX}.kmer_pa.npz" ] || { echo "exists, skipping"; exit 0; }
echo "[$(date)] $CHR kmer_pa v3qc_v3 (--treat-missing-as-n)"
$PY -u $BASE/src/build_kmer_pa.py \
    --kmers "$KMERS" --vcf "$VCF" --ref "$REF" \
    --chrom "$CHR" --out "$OUT_PREFIX" \
    --treat-missing-as-n
echo "[$(date)] $CHR DONE"
ls -lh ${OUT_PREFIX}.kmer_pa.npz ${OUT_PREFIX}.meta.npz

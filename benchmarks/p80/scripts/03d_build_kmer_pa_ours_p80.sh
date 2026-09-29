#!/bin/bash
#SBATCH --job-name=p80_ours
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=logs/03d_ours_%j.out
#SBATCH --error=logs/03d_ours_%j.out

# =============================================================================
# p80 IN-HOUSE-INDEX kmer_pa — matches PRODUCTION construction (in-house `ours`
# k-mer dictionary + --treat-missing-as-n), so p80 vs the new p231 differ ONLY by
# panel composition (apples-to-apples). Builds TWO arms:
#   raw       → data/kmer_pa_p80_ours/kmer_pa_Chr1
#   filt2inv  → data/kmer_pa_p80_ours_filt2inv/kmer_pa_Chr1  (--filter-production)
# The private (ac=1) drop was a guard against the UNEVEN p231 panel; p80 is
# homogeneous, so we benchmark both to confirm it's ~neutral here.
# =============================================================================
mkdir -p logs
set -euo pipefail
CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p80
BASE=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python

KMERS=$BASE/panel/pangenie_index/pang_135_haploid/ours_Chr1_kmers.tsv.gz   # in-house index (prod)
VCF=$CTRL/data/pangenome_p80_chr1.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
for f in "$KMERS" "$VCF" "$REF"; do [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done

RAW=$CTRL/data/kmer_pa_p80_ours
FILT=$CTRL/data/kmer_pa_p80_ours_filt2inv
mkdir -p "$RAW" "$FILT"

echo "[$(date)] p80 in-house kmer_pa — RAW (N-on, no filter)"
$PY -u $BASE/src/build_kmer_pa.py \
    --kmers "$KMERS" --vcf "$VCF" --ref "$REF" --chrom Chr1 \
    --treat-missing-as-n \
    --out "$RAW/kmer_pa_Chr1"

echo "[$(date)] p80 in-house kmer_pa — FILT2INV (--filter-production, 2<=ac<=F-1)"
$PY -u $BASE/src/build_kmer_pa.py \
    --kmers "$KMERS" --vcf "$VCF" --ref "$REF" --chrom Chr1 \
    --treat-missing-as-n --filter-production --min-ac 2 --invariant-margin 1 \
    --out "$FILT/kmer_pa_Chr1"

echo "[$(date)] DONE"; ls -lh "$RAW"/kmer_pa_Chr1.kmer_pa.npz "$FILT"/kmer_pa_Chr1.kmer_pa.npz

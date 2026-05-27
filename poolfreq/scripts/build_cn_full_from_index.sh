#!/bin/bash
#SBATCH --job-name=cnfull_idx
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/cnfull_idx_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/cnfull_idx_%j.err
set -euo pipefail

# Build cn_full Chr1 from a given kmers.tsv.gz index, against the 231 haploid panel.
# Everything downstream of the index is held IDENTICAL to the production PG-index
# build (build_cn_full_v3qc_v3_chr1.sh): same VCF, same ref, --treat-missing-as-n.
# The ONLY thing that varies between runs is the k-mer index -> isolates the
# index builder's effect on cn_full (and hence on the EM h estimate).
#
# Usage: sbatch build_cn_full_from_index.sh <KMERS_TSV_GZ> <OUT_DIR>

KMERS=$1
OUT_DIR=$2

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
CHR=Chr1
VCF=$BASE/pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz
REF=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa

mkdir -p $OUT_DIR
OUT_PREFIX=$OUT_DIR/cn_${CHR}
[ -s "$KMERS" ] || { echo "ERROR: missing index $KMERS"; exit 1; }
[ -s "$VCF" ]   || { echo "ERROR: missing $VCF"; exit 1; }
[ ! -s "${OUT_PREFIX}.cn.npz" ] || { echo "exists, skipping"; exit 0; }

echo "[$(date)] cn_full from index=$KMERS -> $OUT_PREFIX (--treat-missing-as-n)"
$PY -u $BASE/poolfreq/src/build_kmer_cn.py \
    --kmers "$KMERS" --vcf "$VCF" --ref "$REF" \
    --chrom "$CHR" --out "$OUT_PREFIX" \
    --treat-missing-as-n
echo "[$(date)] DONE"; ls -lh ${OUT_PREFIX}.cn.npz ${OUT_PREFIX}.meta.npz

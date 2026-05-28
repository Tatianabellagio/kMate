#!/bin/bash
#SBATCH --job-name=kmidx135d
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=8:00:00
#SBATCH --output=logs/pang_135_diploid_%j.out
#SBATCH --error=logs/pang_135_diploid_%j.err
mkdir -p logs
set -euo pipefail

# Full pang_135 PG-equivalence check: run our build_kmers_tsv.py on the
# diploidized vcfbub-filtered pang_135 VCF (same input PG-index ate on
# 2026-05-01) and produce per-chrom kmers.tsv.gz to byte-compare against
# pangenie_genotyping/data/pang_135_pangenie_index_Chr{1..5}_kmers.tsv.gz.

source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh
conda activate pangenie

WORK=/global/scratch/users/tbellg/hapfire_sv/kmer_index
VCF=/global/scratch/users/tbellg/hapfire_sv/pangenie_genotyping/data/pang_1001gplus_all.dipl.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa
OUT_DIR=$WORK/pang_135_diploid
OUT_PREFIX=$OUT_DIR/ours

mkdir -p $OUT_DIR

echo "[$(date)] start"
echo "  vcf=$VCF"
echo "  ref=$REF"
echo "  out=$OUT_PREFIX"
echo

/usr/bin/time -v python -u $WORK/scripts/build_kmers_tsv.py \
    --vcf $VCF \
    --ref $REF \
    --out $OUT_PREFIX \
    -k 31 \
    --jellyfish-threads 8 \
    --jellyfish-hash 3000000000

echo
echo "[$(date)] DONE"
ls -lh $OUT_DIR/

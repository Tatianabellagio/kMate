#!/bin/bash
#SBATCH --job-name=kmidx135h2x
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=10:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/kmer_index/logs/pang_135_haploid_cap2x_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/kmer_index/logs/pang_135_haploid_cap2x_%j.err
set -euo pipefail

# Raised-cap variant of the haploid index (for the cn_full h-equivalence test).
# Doubles the per-allele unique-kmer caps: biallelic 16->32, multiallelic 32->64.
# Motivation: our index under-emits k-mers on missing-GT-heavy bubbles because
# we drop missing paths (classify biallelic, cap 16) where PG materializes an
# N-allele (classify multiallelic, cap 32). This variant tests whether emitting
# MORE k-mers per allele improves the EM h estimate.

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh
conda activate pangenie

WORK=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/kmer_index
VCF=/home/tbellagio/scratch/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz
REF=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa
OUT_DIR=$WORK/pang_135_haploid_cap2x
OUT_PREFIX=$OUT_DIR/ours

mkdir -p $OUT_DIR
echo "[$(date)] start cap2x (biallelic=32 multiallelic=64) haploid"
echo "  vcf=$VCF"; echo "  out=$OUT_PREFIX"; echo

/usr/bin/time -v python -u $WORK/scripts/build_kmers_tsv.py \
    --vcf $VCF \
    --ref $REF \
    --out $OUT_PREFIX \
    -k 31 \
    --haploid \
    --cap-biallelic 32 \
    --cap-multiallelic 64 \
    --jellyfish-threads 8 \
    --jellyfish-hash 3000000000

echo; echo "[$(date)] DONE"; ls -lh $OUT_DIR/

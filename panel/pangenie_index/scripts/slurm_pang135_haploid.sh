#!/bin/bash
#SBATCH --job-name=kmidx135h
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=8:00:00
#SBATCH --output=logs/pang_135_haploid_%j.out
#SBATCH --error=logs/pang_135_haploid_%j.err
mkdir -p logs
set -euo pipefail

# Parallel sister to slurm_pang135_diploid.sh. Same panel, same reference,
# but input is the RAW haploid-GT vcfbub output (skips the inline diploidize
# step). With --haploid, each inbred sample contributes 1 path of allele X
# instead of 2 paths of X|X. Output should be byte-identical because the SET
# of distinct path-tuples through each bubble is the same in either case.
#
# If diploid vs haploid outputs match: --haploid is production-ready, and
# downstream pipelines can skip the awk diploidize step (build_pangenie_index.sh
# line 41-44) entirely.

source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh
conda activate pangenie

WORK=/global/scratch/users/tbellg/kmate/panel/pangenie_index
VCF=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa
OUT_DIR=$WORK/pang_135_haploid
OUT_PREFIX=$OUT_DIR/ours

mkdir -p $OUT_DIR

echo "[$(date)] start"
echo "  vcf=$VCF (haploid GT, raw vcfbub)"
echo "  ref=$REF"
echo "  out=$OUT_PREFIX"
echo

/usr/bin/time -v python -u $WORK/scripts/build_kmers_tsv.py \
    --vcf $VCF \
    --ref $REF \
    --out $OUT_PREFIX \
    -k 31 \
    --haploid \
    --jellyfish-threads 8 \
    --jellyfish-hash 3000000000

echo
echo "[$(date)] DONE"
ls -lh $OUT_DIR/

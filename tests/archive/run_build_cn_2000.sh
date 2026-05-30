#!/bin/bash
#SBATCH --job-name=build_cn_2k
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/build_cn_2k_%j.out
#SBATCH --error=logs/build_cn_2k_%j.err

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate

/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python src/build_kmer_pa.py \
    --kmers /global/scratch/users/tbellg/kmate/pangenie_test/pangenie_idx_rawv2_Chr1_kmers.tsv.gz \
    --vcf   /global/scratch/users/tbellg/kmate/pangenie_test/raw_vcfbub_lv0_diploid.vcf.gz \
    --ref   /global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa \
    --chrom Chr1 \
    --max-bubbles 2000 \
    --out data/test_chr1_first2000

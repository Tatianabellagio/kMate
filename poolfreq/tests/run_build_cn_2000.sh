#!/bin/bash
#SBATCH --job-name=build_cn_2k
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/build_cn_2k_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/build_cn_2k_%j.err

set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq

/home/tbellagio/miniforge3/envs/hapfm/bin/python src/build_kmer_cn.py \
    --kmers /carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_test/pangenie_idx_rawv2_Chr1_kmers.tsv.gz \
    --vcf   /carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_test/raw_vcfbub_lv0_diploid.vcf.gz \
    --ref   /home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa \
    --chrom Chr1 \
    --max-bubbles 2000 \
    --out data/test_chr1_first2000

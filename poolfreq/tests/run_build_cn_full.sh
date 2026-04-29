#!/bin/bash
#SBATCH --job-name=build_cn_full
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/build_cn_full_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/build_cn_full_%j.err

set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq

# Build cn for full Chr1 first (most informative single chromosome, fastest to validate)
for chrom in Chr1 Chr2 Chr3 Chr4 Chr5; do
    echo "[$(date)] Building cn for $chrom"
    /home/tbellagio/miniforge3/envs/hapfm/bin/python src/build_kmer_cn.py \
        --kmers /carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_test/pangenie_idx_rawv2_${chrom}_kmers.tsv.gz \
        --vcf   /carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_test/raw_vcfbub_lv0_diploid.vcf.gz \
        --ref   /home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa \
        --chrom $chrom \
        --out data/cn_full_${chrom}
done
echo "[$(date)] DONE"
ls -la data/cn_full_*

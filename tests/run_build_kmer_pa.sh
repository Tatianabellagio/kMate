#!/bin/bash
#SBATCH --job-name=build_kmer_pa
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=logs/build_kmer_pa_%j.out
#SBATCH --error=logs/build_kmer_pa_%j.err

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate

# Build kmer_pa for full Chr1 first (most informative single chromosome, fastest to validate)
for chrom in Chr1 Chr2 Chr3 Chr4 Chr5; do
    echo "[$(date)] Building kmer_pa for $chrom"
    /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python src/build_kmer_pa.py \
        --kmers /global/scratch/users/tbellg/kmate/pangenie_test/pangenie_idx_rawv2_${chrom}_kmers.tsv.gz \
        --vcf   /global/scratch/users/tbellg/kmate/pangenie_test/raw_vcfbub_lv0_diploid.vcf.gz \
        --ref   /global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa \
        --chrom $chrom \
        --out data/kmer_pa_${chrom}
done
echo "[$(date)] DONE"
ls -la data/kmer_pa_*

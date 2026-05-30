#!/bin/bash
#SBATCH --job-name=sm_perchrom
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/sm_perchrom_%j.out
#SBATCH --error=logs/sm_perchrom_%j.err

# Validate per_sample_per_chrom.py against the genome-wide v2 result on SEEDMIX_S1.
# 64 GB allocation here is the whole point — proves per-chrom mode fits where the
# genome-wide driver needed 192 GB.

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p results/seedmix_231_v2_perchrom

/usr/bin/time -v /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --kmer-pa-prefix data/kmer_pa_231_v2/kmer_pa \
    --var-pa data/var_pa_231_v2.var_pa.npz \
    --var-meta data/var_pa_231_v2.meta.npz \
    --reads /global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz \
            /global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz \
    --sample SEEDMIX_S1 \
    --out results/seedmix_231_v2_perchrom/SEEDMIX_S1.tsv \
    --threads 8

echo "[$(date)] DONE"
ls -lh results/seedmix_231_v2_perchrom/

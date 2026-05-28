#!/bin/bash
#SBATCH --job-name=sm_perchrom_v3
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/sm_perchrom_v3_%j.out
#SBATCH --error=logs/sm_perchrom_v3_%j.err

# v3 panel test on SEEDMIX_S1 — slope should sit at ~1.0 (no Beagle bias),
# vs v2's 1.43× scale that needed post-hoc calibration.
# Default --block-mode global — SEEDMIX is homogeneous F0 (no recomb).

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p results/seedmix_231_v3_perchrom

/usr/bin/time -v /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3/cn \
    --cn-var data/cn_var_231_v3.cn_var.npz \
    --cn-var-meta data/cn_var_231_v3.meta.npz \
    --reads /global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz \
            /global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz \
    --sample SEEDMIX_S1 \
    --out results/seedmix_231_v3_perchrom/SEEDMIX_S1.tsv \
    --threads 8

echo "[$(date)] DONE"
ls -lh results/seedmix_231_v3_perchrom/

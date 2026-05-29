#!/bin/bash
#SBATCH --job-name=seedmix_82
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/seedmix_82_%j.out
#SBATCH --error=logs/seedmix_82_%j.err

# Run per_sample_driver on all 8 SEEDMIX replicates with the 82-founder cn.
# Outputs go to results/seedmix_82/. Restartable: skips samples whose output
# already exists.

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p results/seedmix_82

/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/batch_runner.py \
    --manifest data/seedmix_manifest.tsv \
    --out-dir results/seedmix_82 \
    --cn-kmer-prefix data/cn_full \
    --cn-var data/cn_var_82.cn_var.npz \
    --cn-var-meta data/cn_var_82.meta.npz \
    --threads 8 \
    --workers 1 \
    --block-mode global   # SEEDMIX is F0 (single mixture) → global is correct

#!/bin/bash
#SBATCH --job-name=seedmix_231
#SBATCH --account=fc_moilab
#SBATCH --partition=savio3_bigmem
#SBATCH --qos=savio_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=12:00:00
#SBATCH --output=logs/seedmix_231_%j.out
#SBATCH --error=logs/seedmix_231_%j.err

# Run per_sample_driver on all 8 SEEDMIX replicates with the 231-founder cn.
# Submit with: sbatch --dependency=afterok:<cn_231_jobid> run_seedmix_231.sh

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate/poolfreq
mkdir -p results/seedmix_231

/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/batch_runner.py \
    --manifest data/seedmix_manifest.tsv \
    --out-dir results/seedmix_231 \
    --cn-kmer-prefix data/cn_full_231/cn \
    --cn-var data/cn_var_231.cn_var.npz \
    --cn-var-meta data/cn_var_231.meta.npz \
    --threads 8 \
    --workers 1 \
    --block-mode global   # SEEDMIX is F0 (single mixture) → global is correct

#!/bin/bash
#SBATCH --job-name=site04_cem
#SBATCH --account=fc_moilab
#SBATCH --partition=savio3_bigmem
#SBATCH --qos=savio_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=logs/site04_cem_%j.out
#SBATCH --error=logs/site04_cem_%j.err

# Run cactus_em (231-kmer_pa v2) on site 04 samples (57 samples).
# block-mode window — evolved samples have local mosaic ancestry, so window > global.
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p results/site04_231_v2 tests/logs

/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/batch_runner.py \
    --manifest data/site04_manifest.tsv \
    --out-dir results/site04_231_v2 \
    --kmer-pa-prefix data/kmer_pa_231_v2/kmer_pa \
    --var-pa data/var_pa_231_v2.var_pa.npz \
    --var-meta data/var_pa_231_v2.meta.npz \
    --threads 8 \
    --workers 4 \
    --block-mode window

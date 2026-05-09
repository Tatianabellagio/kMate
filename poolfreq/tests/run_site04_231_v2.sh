#!/bin/bash
#SBATCH --job-name=site04_cem
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/site04_cem_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/site04_cem_%j.err

# Run cactus_em (231-cn v2) on site 04 samples (57 samples).
# block-mode window — evolved samples have local mosaic ancestry, so window > global.
set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
mkdir -p results/site04_231_v2 tests/logs

/home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/batch_runner.py \
    --manifest data/site04_manifest.tsv \
    --out-dir results/site04_231_v2 \
    --cn-kmer-prefix data/cn_full_231_v2/cn \
    --cn-var data/cn_var_231_v2.cn_var.npz \
    --cn-var-meta data/cn_var_231_v2.meta.npz \
    --threads 8 \
    --workers 4 \
    --block-mode window

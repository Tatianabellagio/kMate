#!/bin/bash
#SBATCH --job-name=build_hybrid
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=logs/build_hybrid_%j.out
#SBATCH --error=logs/build_hybrid_%j.err

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate

PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

echo "[$(date)] building hybrid panel"
$PYTHON -u poolfreq/scripts/build_block_hybrid_panel.py \
    --clean-npz poolfreq/data/block_haplotype_cn/chr1_full_clean_cc5.npz \
    --cn-full-prefix poolfreq/data/cn_full_231_v2/cn \
    --chrom Chr1 \
    --out poolfreq/data/block_haplotype_cn/chr1_full_hybrid_cc5.npz

echo "[$(date)] DONE"
ls -lh poolfreq/data/block_haplotype_cn/chr1_full_hybrid_cc5.npz

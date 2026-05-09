#!/bin/bash
#SBATCH --job-name=build_hybrid
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/block_haplotype_cn/build_hybrid_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/block_haplotype_cn/build_hybrid_%j.err

set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv

PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python

echo "[$(date)] building hybrid panel"
$PYTHON -u poolfreq/scripts/build_block_hybrid_panel.py \
    --clean-npz poolfreq/data/block_haplotype_cn/chr1_full_clean_cc5.npz \
    --cn-full-prefix poolfreq/data/cn_full_231_v2/cn \
    --chrom Chr1 \
    --out poolfreq/data/block_haplotype_cn/chr1_full_hybrid_cc5.npz

echo "[$(date)] DONE"
ls -lh poolfreq/data/block_haplotype_cn/chr1_full_hybrid_cc5.npz

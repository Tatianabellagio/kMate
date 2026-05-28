#!/bin/bash
#SBATCH --job-name=build_bldhap
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=logs/build_chr1_full_%j.out
#SBATCH --error=logs/build_chr1_full_%j.err

# Build per-block haplotype-level cn matrix for full Chr1.
# One-time job per panel; produces poolfreq/data/block_haplotype_cn/chr1_full.npz.
mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate

PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
echo "[$(date)] starting build_block_haplotype_cn for Chr1 (numpy 2-bit hashing, 8 threads)"
$PYTHON -u poolfreq/scripts/build_block_haplotype_cn.py \
    --block-index sims/visor_freqk/chr1_only_panel/hapfire_block_index_chr1.npz \
    --fastas-dir sims/visor_freqk/founder_fastas_231 \
    --out poolfreq/data/block_haplotype_cn/chr1_full.npz \
    --chrom-filter Chr1 --threads 8

echo "[$(date)] DONE"
ls -lh poolfreq/data/block_haplotype_cn/chr1_full.npz

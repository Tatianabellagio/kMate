#!/bin/bash
#SBATCH --job-name=blocks_p80
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/blocks_p80_%j.out
#SBATCH --error=logs/blocks_p80_%j.out
# LD-block partition on the p80 (80 long-read founder) panel, Chr1 — same recipe as
# the 231 build (CLQcut 0.9, min-called-frac 0.9 => >=72/80 founders), so we can test
# whether thin-LD-block h-retrieval problems are fundamental or a 231-panel-imbalance
# artifact. Coarse CompleteLDPartition (corr 0.2) + fine BigLD (gpart). r_env + LD_LIBRARY_PATH.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/grenenet_selection/blocks_mcf90_p80
export LD_LIBRARY_PATH=/usr/lib64:${LD_LIBRARY_PATH:-}
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
RSCRIPT=/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript

echo "[$(date)] p80 Chr1 : BigLD CLQcut=0.9, min-called-frac=0.9 (>=72/80)"
$PY analysis/grenenet_selection/blocks/recompute_blocks.py \
    --var-pa-prefix benchmarks/p80/data/var_pa_p80 \
    --chrom Chr1 --corr 0.2 --window 50 --maf 0.05 --clqcut 0.9 --min-called-frac 0.9 \
    --out analysis/grenenet_selection/blocks_mcf90_p80/chr1_clq0.9 \
    --rscript "$RSCRIPT" --bigld-dir /global/scratch/users/tbellg/HapFM/bin
echo "[$(date)] p80 Chr1 done"

#!/bin/bash
#SBATCH --job-name=clqsweep90
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/clqsweep90_%A_%a.out
#SBATCH --error=logs/clqsweep90_%A_%a.out
# Regenerate Parts A/B/C on the NEW non-missing (min-called-frac 0.9) map, Chr1:
#  task0: CLQ0.5 blocks + frontier   task1: CLQ0.7 blocks + frontier   (A & B)
#  task2: founder-vs-evolved on the existing CLQ0.9 mcf90 blocks        (C)
# (CLQ0.9 mcf90 blocks already exist in blocks_mcf90.)

set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/grenenet_gea/blocks_mcf90
export LD_LIBRARY_PATH=/usr/lib64:${LD_LIBRARY_PATH:-}
export MINCF=0.9
export BLOCKS_DIR=analysis/grenenet_gea/blocks_mcf90
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
RSCRIPT=/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript
BIGLD=/global/scratch/users/tbellg/HapFM/bin
VP=panel/arch3/chr1/var_pa_231_arch3_chr1

case $SLURM_ARRAY_TASK_ID in
  0|1)
    CLQ=$([ "$SLURM_ARRAY_TASK_ID" = "0" ] && echo 0.5 || echo 0.7)
    echo "[$(date)] Chr1 CLQ$CLQ blocks (min-called-frac 0.9)"
    $PY analysis/grenenet_gea/recompute_blocks.py --var-pa-prefix "$VP" \
        --chrom Chr1 --corr 0.2 --window 50 --maf 0.05 --clqcut "$CLQ" --min-called-frac 0.9 \
        --out analysis/grenenet_gea/blocks_mcf90/chr1_clq${CLQ} \
        --rscript "$RSCRIPT" --bigld-dir "$BIGLD"
    echo "[$(date)] Chr1 CLQ$CLQ frontier (clustering)"
    $PY -u analysis/grenenet_gea/block_unit_frontier.py --clqcut "$CLQ" --chrom Chr1
    ;;
  2)
    echo "[$(date)] Chr1 founder-vs-evolved on CLQ0.9 mcf90"
    $PY -u analysis/grenenet_gea/block_founder_vs_evolved.py --chrom Chr1 --clqcut 0.9
    ;;
esac
echo "[$(date)] task $SLURM_ARRAY_TASK_ID done"

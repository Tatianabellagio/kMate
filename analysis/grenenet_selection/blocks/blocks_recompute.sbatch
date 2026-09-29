#!/bin/bash
#SBATCH --job-name=blocks_bigld
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/blocks_bigld_%A_%a.out
#SBATCH --error=logs/blocks_bigld_%A_%a.out
# FAITHFUL HapFM recompute on OUR all-class founder var_pa:
# coarse CompleteLDPartition (corr) + fine BigLD (gpart, CLQmode=density), sweep CLQcut.
# gpart::BigLD lives in r_env (installed from stripped source); BigLD.R needs optparse
# (installed) and LD_LIBRARY_PATH=/usr/lib64 (system libxml2 for igraph).

set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/grenenet_selection/blocks/results/blocks_recompute
export LD_LIBRARY_PATH=/usr/lib64:${LD_LIBRARY_PATH:-}

CLQ=(0.5 0.7 0.9)
c=${CLQ[$SLURM_ARRAY_TASK_ID]}
CHR=${CHR:-Chr1}
chrlc=$(echo "$CHR" | tr 'A-Z' 'a-z')
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
RSCRIPT=/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript

echo "[$(date)] $CHR coarse corr=0.2 + BigLD CLQcut=$c"
$PY analysis/grenenet_selection/blocks/recompute_blocks.py \
    --var-pa-prefix panel/arch3/${chrlc}/var_pa_231_arch3_${chrlc} \
    --chrom "$CHR" --corr 0.2 --window 50 --maf 0.05 --clqcut "$c" \
    --out analysis/grenenet_selection/blocks/results/blocks_recompute/${chrlc}_clq${c} \
    --rscript "$RSCRIPT" --bigld-dir /global/scratch/users/tbellg/HapFM/bin
echo "[$(date)] done CLQcut=$c"

#!/bin/bash
#SBATCH --job-name=blocks_mcf90
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/blocks_mcf90_%A_%a.out
#SBATCH --error=logs/blocks_mcf90_%A_%a.out
# LOCKED rebuild: CLQcut=0.9, panel-support filter --min-called-frac 0.9 (keep variants
# genotyped in >=208/231 founders), all 5 chroms (array=chrom). Coarse CompleteLDPartition
# (corr=0.2) + fine BigLD (gpart, density). Writes to a NEW dir so the old min-called-frac
# 0.5 map is preserved for comparison. gpart::BigLD lives in r_env; needs LD_LIBRARY_PATH.

set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs results/grenenet_gea/blocks_mcf90
export LD_LIBRARY_PATH=/usr/lib64:${LD_LIBRARY_PATH:-}

CHRS=(Chr1 Chr2 Chr3 Chr4 Chr5)
CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}
chrlc=$(echo "$CHR" | tr 'A-Z' 'a-z')
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
RSCRIPT=/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript

echo "[$(date)] $CHR : BigLD CLQcut=0.9, min-called-frac=0.9"
$PY analysis/grenenet_gea/recompute_blocks.py \
    --var-pa-prefix panel/arch3/${chrlc}/var_pa_231_arch3_${chrlc} \
    --chrom "$CHR" --corr 0.2 --window 50 --maf 0.05 --clqcut 0.9 --min-called-frac 0.9 \
    --out results/grenenet_gea/blocks_mcf90/${chrlc}_clq0.9 \
    --rscript "$RSCRIPT" --bigld-dir /global/scratch/users/tbellg/HapFM/bin
echo "[$(date)] $CHR done"

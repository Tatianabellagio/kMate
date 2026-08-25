#!/bin/bash
#SBATCH --job-name=gemma_meta
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=logs/gemma_meta_%j.out
#SBATCH --error=logs/gemma_meta_%j.out
# Cross-site meta over the GEMMA per-site Z (Bolormaa et al. 2014 JOINT/GLOBAL/CLIMATE), the
# same formulas class_split_gwas.py uses -- GEMMA cannot do this stage (no multi-phenotype
# association test in GEMMA/GCTA/EMMAX at 30 traits). Run AFTER run_gemma_gwas.sbatch finishes:
#   sbatch --dependency=afterok:<array_jobid> run_gemma_meta.sbatch
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate
mkdir -p logs
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
echo "[$(date)] gemma meta start on $(hostname)"
$PY analysis/grenenet_selection/r3_persite_gwas/gemma_persite_gwas.py meta
echo "[$(date)] done"

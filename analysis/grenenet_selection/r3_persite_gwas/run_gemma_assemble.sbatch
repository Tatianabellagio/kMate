#!/bin/bash
#SBATCH --job-name=gemma_assemble
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=logs/gemma_assemble_%j.out
#SBATCH --error=logs/gemma_assemble_%j.out
# Concatenate the per-chrom GEMMA parts and attach PER-GARDEN significance (BH + Bonferroni
# within each garden x class). There is NO cross-site meta: scope decision 2026-08-25 is
# per-garden only, since the all-sites question is answered by the GEA track. Every statistic
# this produces comes from GEMMA. Run AFTER run_gemma_gwas.sbatch finishes:
#   sbatch --dependency=afterok:<array_jobid> run_gemma_assemble.sbatch
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate
mkdir -p logs
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
echo "[$(date)] gemma assemble start on $(hostname)"
$PY analysis/grenenet_selection/r3_persite_gwas/gemma_persite_gwas.py meta
echo "[$(date)] done"

#!/bin/bash
#SBATCH --job-name=gemma_gwas
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --array=0-14
#SBATCH --output=logs/gemma_gwas_%A_%a.out
#SBATCH --error=logs/gemma_gwas_%A_%a.out
# Per-site class-split founder GWAS run with GEMMA (replaces the in-house EMMAX/P3D scan).
# 15 array tasks = 5 chroms x 3 classes. Each writes its BIMBAM genotype ONCE and reuses it
# across all 30 sites (phenotype is a -n column switch), so the per-site cost is just GEMMA.
# Markers are chunked below 20,000 per run: GEMMA 0.98.5 crashes with
#   "Enforce failed for N>0 in src/fastblas.cpp ... fast_cblas_dgemm"
# on any run at or above 20k markers -- after finishing the scan, before writing output.
# After the array completes, run the meta stage:  sbatch run_gemma_meta.sbatch
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate
mkdir -p logs analysis/grenenet_selection/r3_persite_gwas/results/gemma_gwas

CHROMS=(chr1 chr2 chr3 chr4 chr5)
CLASSES=(snp nonsnp sv)
CH=${CHROMS[$((SLURM_ARRAY_TASK_ID / 3))]}
CLS=${CLASSES[$((SLURM_ARRAY_TASK_ID % 3))]}

PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
export PATH=/global/home/users/tbellg/miniforge3/envs/gwas_tools/bin:$PATH
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4

echo "[$(date)] gemma GWAS task $SLURM_ARRAY_TASK_ID -> $CH / $CLS on $(hostname)"
$PY analysis/grenenet_selection/r3_persite_gwas/gemma_persite_gwas.py scan --chrom "$CH" --class "$CLS"
echo "[$(date)] done $CH / $CLS"

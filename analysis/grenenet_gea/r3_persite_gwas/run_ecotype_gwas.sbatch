#!/bin/bash
#SBATCH --job-name=eco_gwas
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --output=logs/eco_gwas_%j.out
#SBATCH --error=logs/eco_gwas_%j.out
# Tier 2 — 231-ecotype kinship-corrected fitness-GWAS. Two streaming passes over the 5
# per-chrom var_pa npz: (1) SNP GRM K, (2) EMMAX-style vectorized LMM Wald z per variant on
# 8 fitness axes ({relative,census} x {global,cold,mid,hot}). See ecotype_fitness_gwas.py.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/grenenet_gea/ecotype_fitness/gwas
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
echo "[$(date)] ecotype fitness GWAS start on $(hostname)"
$PY analysis/grenenet_gea/r3_persite_gwas/ecotype_fitness_gwas.py
echo "[$(date)] done"

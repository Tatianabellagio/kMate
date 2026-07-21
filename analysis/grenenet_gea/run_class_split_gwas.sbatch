#!/bin/bash
#SBATCH --job-name=class_split_gwas
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=2:00:00
#SBATCH --output=logs/class_split_gwas_%j.out
#SBATCH --error=logs/class_split_gwas_%j.out
# Class-split (SNP/non-SNP/SV) single-marker founder GWAS -- ~2.3M markers x 31 sites.
# Added 2026-07-21: the previous ad-hoc interactive run OOM-killed under the 16G interactive
# session cgroup (only 12.5GB RSS before the kill) -- this needs a real sbatch allocation,
# not the login/interactive-shell cgroup (matches project convention for heavy kMate compute).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/grenenet_gea/varexp
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
echo "[$(date)] class-split GWAS start on $(hostname)"
$PY analysis/grenenet_gea/class_split_gwas.py
echo "[$(date)] done"

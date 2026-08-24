#!/bin/bash
#SBATCH --job-name=class_grms
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=logs/class_grms_%j.out
#SBATCH --error=logs/class_grms_%j.out
# Per-variant-class founder GRMs (SNP / indel / SV / non-SNP / all / MAC-matched-SNP) for the
# ecotype-selection variance partition. Two streaming passes over the 5 per-chrom var_pa npz.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/grenenet_selection/r3_persite_gwas/results/varexp
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
echo "[$(date)] class GRMs start on $(hostname)"
$PY analysis/grenenet_selection/r3_persite_gwas/build_class_grms.py
echo "[$(date)] done"

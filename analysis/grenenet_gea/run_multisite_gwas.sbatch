#!/bin/bash
#SBATCH --job-name=multisite_gwas
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:30:00
#SBATCH --output=logs/multisite_gwas_%j.out
#SBATCH --error=logs/multisite_gwas_%j.out
# Multi-trait founder GWAS across sites (joint/global/climate), gen1 trait, site 33 excluded -> 30 sites.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] multisite founder GWAS on $(hostname) | TRAIT_GENS=${TRAIT_GENS:-1,2,3} EXCLUDE_SITES=${EXCLUDE_SITES:-33}"
$PY analysis/grenenet_gea/founder_gwas_multisite.py
echo "[$(date)] done"

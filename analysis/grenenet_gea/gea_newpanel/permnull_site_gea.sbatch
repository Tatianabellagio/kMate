#!/bin/bash
#SBATCH --job-name=permnull_gea
#SBATCH --partition=savio4_htc
#SBATCH --account=co_moilab
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/results/grenenet_gea/gea_newpanel/permnull/permnull_%j.log

set -euo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd /global/scratch/users/tbellg/kmate
mkdir -p results/grenenet_gea/gea_newpanel/permnull
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python \
  analysis/grenenet_gea/gea_newpanel/permnull_site_gea.py

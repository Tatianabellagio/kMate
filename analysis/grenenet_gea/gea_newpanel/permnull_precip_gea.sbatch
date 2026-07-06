#!/bin/bash
#SBATCH --job-name=permnull_precip
#SBATCH --partition=savio4_htc
#SBATCH --account=co_moilab
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/results/grenenet_gea/gea_newpanel/permnull_precip/permnull_precip_%j.log

set -euo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd /global/scratch/users/tbellg/kmate
mkdir -p results/grenenet_gea/gea_newpanel/permnull_precip

# {snp,nonsnp,sv} x {bio12,bio13,bio16,bio19} x K in {1,3}, 2000 site-permutations
export PERMNULL_NPERM=2000
export PERMNULL_KS=1,3
export PERMNULL_AXES=bio12,bio13,bio16,bio19
export PERMNULL_CLASSES=snp,nonsnp,sv

/global/home/users/tbellg/miniforge3/envs/kmate/bin/python \
  analysis/grenenet_gea/gea_newpanel/permnull_precip_gea.py

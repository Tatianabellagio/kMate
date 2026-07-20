#!/bin/bash
#SBATCH --job-name=seedmix_all
#SBATCH --partition=savio4_htc
#SBATCH --account=co_moilab
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/seedmix_validation/fix_seedmix/allsamples_%j.log

source ~/miniforge3/etc/profile.d/conda.sh
conda activate kmate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# DBs are big (8 x ~2GB); keep on shared scratch so they persist for reruns
export JF_DIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/seedmix_validation/fix_seedmix
cd /global/scratch/users/tbellg/kmate
python analysis/grenenet_gea/seedmix_validation/fix_seedmix/run_seedmix_allsamples.py

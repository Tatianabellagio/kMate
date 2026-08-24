#!/bin/bash
#SBATCH --job-name=emreg_sweep
#SBATCH --partition=savio4_htc
#SBATCH --account=co_moilab
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/qc/seedmix_validation/fix_emreg/sweep_%j.log

source ~/miniforge3/etc/profile.d/conda.sh
conda activate kmate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /global/scratch/users/tbellg/kmate
python analysis/grenenet_gea/qc/seedmix_validation/fix_emreg/emreg_sweep.py

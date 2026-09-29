#!/bin/bash
#SBATCH --job-name=kfw_synth
#SBATCH --partition=savio4_htc
#SBATCH --account=co_moilab
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/qc/seedmix_validation/fix_kfw_fullpanel/synthetic_%j.log
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate kmate
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd /global/scratch/users/tbellg/kmate
python analysis/grenenet_selection/qc/seedmix_validation/fix_kfw_fullpanel/run_synthetic_compare.py

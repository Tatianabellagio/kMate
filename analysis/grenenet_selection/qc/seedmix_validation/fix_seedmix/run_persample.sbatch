#!/bin/bash
#SBATCH --job-name=seedmix_ps
#SBATCH --partition=savio4_htc
#SBATCH --account=co_moilab
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --array=0-7
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/qc/seedmix_validation/fix_seedmix/ps_%A_%a.log
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
conda activate kmate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
SAMPLES=(SEEDMIX_S1 SEEDMIX_S2 SEEDMIX_S3 SEEDMIX_S4 SEEDMIX_S5 SEEDMIX_S6 SEEDMIX_S7 SEEDMIX_S8)
SID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
cd /global/scratch/users/tbellg/kmate
python analysis/grenenet_selection/qc/seedmix_validation/fix_seedmix/run_persample.py $SID

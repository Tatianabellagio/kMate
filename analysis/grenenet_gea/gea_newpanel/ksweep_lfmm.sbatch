#!/bin/bash
#SBATCH --job-name=ksweep_lfmm
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --array=0-2
#SBATCH --output=/global/scratch/users/tbellg/kmate/results/grenenet_gea/gea_newpanel/ksweep_calib/ksweep_%a_%j.log
set -euo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd /global/scratch/users/tbellg/kmate
CLASSES=(snp nonsnp sv)
CLS=${CLASSES[$SLURM_ARRAY_TASK_ID]}
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
$R analysis/grenenet_gea/gea_newpanel/ksweep_lfmm.R "$CLS"

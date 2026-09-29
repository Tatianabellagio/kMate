#!/bin/bash
#SBATCH --job-name=sv_freqrob
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --output=logs/sv_freqrob_%j.out
#SBATCH --error=logs/sv_freqrob_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
export MEMB_TAG=clq90 OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
echo "[$(date)] host=$(hostname) start"
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python analysis/grenenet_selection/extras/_sv_hap_freqrobust.py
echo "[$(date)] done"

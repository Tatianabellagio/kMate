#!/bin/bash
#SBATCH --job-name=sv_hap_axes
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=logs/sv_hap_axes_%j.out
#SBATCH --error=logs/sv_hap_axes_%j.out
# Correct-UNIT (hap-cluster) SV enrichment across JOINT/CLIMATE/TEMPORAL axes,
# MAF sweep (SV_MAC 2..46) + power + positive control. Env: kmate (MEMB_TAG=clq90).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
export MEMB_TAG=clq90
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
echo "[$(date)] host=$(hostname) start"
$PY analysis/grenenet_gea/_sv_haplotype_axes_sweep.py
echo "[$(date)] done"

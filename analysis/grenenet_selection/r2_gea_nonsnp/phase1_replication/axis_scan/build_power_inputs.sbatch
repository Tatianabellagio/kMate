#!/bin/bash
#SBATCH --job-name=power_inputs
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=1:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/axis_scan/power_inputs_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/axis_scan/power_inputs_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python -u analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/axis_scan/build_power_inputs.py
echo DONE_POWER_INPUTS

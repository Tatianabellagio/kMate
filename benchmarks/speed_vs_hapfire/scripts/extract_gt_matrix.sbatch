#!/bin/bash
#SBATCH --job-name=gt_matrix
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/gt_matrix_%j.out
#SBATCH --error=logs/gt_matrix_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/benchmarks/speed_vs_hapfire/scripts
/global/home/users/tbellg/miniforge3/envs/basic/bin/python extract_greneNet_gt_matrix.py

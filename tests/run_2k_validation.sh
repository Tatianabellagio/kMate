#!/bin/bash
#SBATCH --job-name=val_2k
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/val_2k_%j.out
#SBATCH --error=logs/val_2k_%j.err

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python tests/test_2k_validation.py

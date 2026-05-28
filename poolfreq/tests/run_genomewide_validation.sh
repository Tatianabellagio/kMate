#!/bin/bash
#SBATCH --job-name=val_genomewide
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=logs/val_gw_%j.out
#SBATCH --error=logs/val_gw_%j.err

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate/poolfreq
/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u tests/test_genomewide_validation.py

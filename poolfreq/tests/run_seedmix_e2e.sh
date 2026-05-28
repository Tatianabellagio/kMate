#!/bin/bash
#SBATCH --job-name=seedmix_e2e
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/seedmix_e2e_%j.out
#SBATCH --error=logs/seedmix_e2e_%j.err

set -uo pipefail
mkdir -p /global/scratch/users/tbellg/kmate/poolfreq/tests/logs
cd /global/scratch/users/tbellg/kmate/poolfreq
/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python tests/test_seedmix_e2e.py

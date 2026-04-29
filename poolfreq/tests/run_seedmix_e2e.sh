#!/bin/bash
#SBATCH --job-name=seedmix_e2e
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/seedmix_e2e_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/seedmix_e2e_%j.err

set -uo pipefail
mkdir -p /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
/home/tbellagio/miniforge3/envs/hapfm/bin/python tests/test_seedmix_e2e.py

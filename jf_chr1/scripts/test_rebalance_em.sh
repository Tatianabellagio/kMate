#!/bin/bash
#SBATCH --job-name=test_rebalance
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=0:30:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/jf_chr1/logs/rebalance_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/jf_chr1/logs/rebalance_%j.err

set -euo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv
/home/tbellagio/miniforge3/envs/hapfm/bin/python -u jf_chr1/scripts/test_rebalance_em.py

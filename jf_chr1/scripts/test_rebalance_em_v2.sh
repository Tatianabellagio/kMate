#!/bin/bash
#SBATCH --job-name=rebal_v2
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=01:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/jf_chr1/logs/rebal_v2_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/jf_chr1/logs/rebal_v2_%j.err

set -euo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv
/home/tbellagio/miniforge3/envs/hapfm/bin/python -u jf_chr1/scripts/test_rebalance_em_v2.py

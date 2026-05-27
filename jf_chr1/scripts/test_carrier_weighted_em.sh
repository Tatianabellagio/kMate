#!/bin/bash
#SBATCH --job-name=cw_em
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=00:45:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/jf_chr1/logs/cw_em_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/jf_chr1/logs/cw_em_%j.err

set -euo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv
/home/tbellagio/miniforge3/envs/hapfm/bin/python -u jf_chr1/scripts/test_carrier_weighted_em.py

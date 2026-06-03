#!/bin/bash
#SBATCH --job-name=cw_em
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=00:45:00
#SBATCH --output=logs/cw_em_%j.out
#SBATCH --error=logs/cw_em_%j.err

mkdir -p logs
set -euo pipefail
cd /global/scratch/users/tbellg/hapfire_sv
/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u jf_chr1/scripts/test_carrier_weighted_em.py

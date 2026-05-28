#!/bin/bash
#SBATCH --job-name=locicov
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=0:40:00
#SBATCH --output=logs/locicov_%j.out
#SBATCH --error=logs/locicov_%j.err
mkdir -p logs
/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u /global/scratch/users/tbellg/hapfire_sv/poolfreq/tests/loci_cov.py

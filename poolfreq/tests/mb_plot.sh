#!/bin/bash
#SBATCH --job-name=mbplot
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=0:40:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/diag/mbplot_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/diag/mbplot_%j.err
/home/tbellagio/miniforge3/envs/hapfm/bin/python -u /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/mb_plot.py

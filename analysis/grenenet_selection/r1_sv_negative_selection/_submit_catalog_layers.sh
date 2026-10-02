#!/bin/bash
#SBATCH --job-name=cat_layers
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/cat_layers_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python _sv_event_catalog.py

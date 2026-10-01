#!/bin/bash
#SBATCH --job-name=regen_gea_ctx
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/regen_gea_ctx_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python _sv_context_decomposition.py

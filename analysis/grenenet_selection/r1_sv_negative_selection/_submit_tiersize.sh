#!/bin/bash
#SBATCH --job-name=burden_tiersize
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/burden_tiersize_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export STRATUM=tiersize SV_CALL_MIN=0.5 MIN_CAR=2
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python _sv_burden_strata.py

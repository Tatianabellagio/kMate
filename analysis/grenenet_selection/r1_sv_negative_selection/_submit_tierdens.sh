#!/bin/bash
# Which axis carries the sign flip -- how well an event is genotyped, or where it sits?
# Coding events split 2x2: call-rate tier (high >=208/231, low 116-207) x local gene density
# (rich >=212 genes/Mb, poor below). If location is what matters, lo_rich behaves like the high-call
# events; if genotyping quality is what matters, both lo_* bins go positive together.
#SBATCH --job-name=burden_tierdens
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/burden_tierdens_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export STRATUM=tierdens SV_CALL_MIN=0.5 MIN_CAR=2
echo "[$(date)] host=$(hostname) building 2x2 tier x gene-density burdens"
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python _sv_burden_strata.py

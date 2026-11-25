#!/bin/bash
# THE decisive call-rate test. Build SV burdens at the 0.5 floor, but split each context into two tiers:
#   highcall  genotyped in >= 208/231 founders  -- exactly the set the reported 0.9 model uses
#   lowcall   genotyped in 116-207 founders     -- what relaxing the floor newly admits
# Both tiers enter the SAME model, so the question is asked directly: does the reported coding effect
# survive when the newly admitted events sit beside it, or does it move onto them?
#SBATCH --job-name=burden_crtier
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/burden_crtier_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export STRATUM=callrate SV_CALL_MIN=0.5 MIN_CAR=2
echo "[$(date)] host=$(hostname) building call-rate-tiered burdens"
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python _sv_burden_strata.py

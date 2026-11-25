#!/bin/bash
# OVERCOUNT parameterisation: an event contributes to EVERY annotation layer it overlaps, instead of
# being assigned one class by the first-match cascade. Removes the arbitrary TE-vs-promoter ordering;
# each burden becomes "how many SVs touch this class" rather than "touch it and nothing else".
#SBATCH --job-name=overcount
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/overcount_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export CTX_MODE=overcount
echo "[$(date)] build overcount burdens"
SV_UNIT=event /global/home/users/tbellg/miniforge3/envs/kmate/bin/python _sv_burden_build.py
echo "[$(date)] fit all6 + every interaction on the overcount burdens"
export POOL_SET=all6 INTER_SET=all BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all
export FINAL_OUT=final_models_pooled_all BURDEN_MODE=rate CRTAG=_overcount
export EXTRA_FIXED="obio1c + obio1c:bio1c" EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid
for t in 0 1 2 3 4 5 6 7 8; do TASK=$t /global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R; done
echo "[$(date)] DONE"

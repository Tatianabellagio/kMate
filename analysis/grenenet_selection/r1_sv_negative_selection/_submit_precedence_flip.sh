#!/bin/bash
# Does the positive TE / intergenic effect follow the VARIANTS or the LABEL?
# 26% of events called TE also sit in a promoter window, and TE-above-promoter is an arbitrary choice.
# Rebuild the burdens with that one ordering flipped, refit, and compare.
#SBATCH --job-name=prec_flip
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/prec_flip_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export CTX_PRECEDENCE=promoter_first
echo "[$(date)] rebuild burdens with promoter above TE"
SV_UNIT=event /global/home/users/tbellg/miniforge3/envs/kmate/bin/python _sv_burden_build.py
echo "[$(date)] fit all6 + every interaction on the flipped burdens"
export POOL_SET=all6 INTER_SET=all BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all
export FINAL_OUT=final_models_pooled_all BURDEN_MODE=rate CRTAG=_promfirst
export EXTRA_FIXED="obio1c + obio1c:bio1c" EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid
for t in 0 1 2 3 4 5 6 7 8; do TASK=$t /global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R; done
echo "[$(date)] DONE"

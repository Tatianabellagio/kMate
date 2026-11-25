#!/bin/bash
# Permutation null for the burden x garden-temperature interactions of THE model (all6, every context
# interacted, 0.9 floor, carried/called). Each replicate shuffles which temperature belongs to which
# garden and refits: founders keep their burdens and kinship, gardens keep their own mean fitness, and
# only the garden-identity/temperature pairing is destroyed. Gives a calibrated null for the z statistic
# of every interaction without relying on lme4's asymptotics.
#SBATCH --job-name=burden_perm
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=0-99
#SBATCH --output=logs/burden_perm_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=all6 INTER_SET=all BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all
export FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid CRTAG=""
export TASK=1                      # M2p_reml -- the reported fit
for k in $(seq 0 9); do
  export PERM_SEED=$(( SLURM_ARRAY_TASK_ID * 10 + k + 1 ))
  echo "[$(date)] seed $PERM_SEED"
  /global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
done

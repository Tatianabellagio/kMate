#!/bin/bash
# Same permutation null for the cds2 arm (coding vs everything else), so sv_burden_cds_vs_rest.ipynb
# can report permutation p-values for its two interactions too.
#SBATCH --job-name=burden_perm2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=0-99
#SBATCH --output=logs/burden_perm2_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=cds2 BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid CRTAG="" BURDEN_MODE=rate
export TASK=1
for k in $(seq 0 9); do
  export PERM_SEED=$(( SLURM_ARRAY_TASK_ID * 10 + k + 1 ))
  /global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
done

#!/bin/bash
# Permutation nulls for both reported arms, rebuilt under the 1 kb promoter rule.
#SBATCH --job-name=perm_ctx
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=10:00:00
#SBATCH --array=0-199
#SBATCH --output=logs/perm_ctx_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid CRTAG="" BURDEN_MODE=rate TASK=1
i=$SLURM_ARRAY_TASK_ID
if [ "$i" -lt 100 ]; then export POOL_SET=all6 INTER_SET=all; k=$i
else                       export POOL_SET=cds2;              k=$(( i - 100 )); fi
for j in $(seq 0 9); do
  export PERM_SEED=$(( k * 10 + j + 1 ))
  /global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
done

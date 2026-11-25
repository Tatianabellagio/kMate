#!/bin/bash
# Every genomic context as its own burden AND its own garden-temperature interaction, rate vs count at
# the 0.9 floor. Asks whether the promoter x temperature signal holds when UTR and intron are in the
# model with their own terms, instead of being absent entirely.
#SBATCH --job-name=burden_all6i
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-17
#SBATCH --output=logs/burden_all6i_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
i=$SLURM_ARRAY_TASK_ID
if [ "$i" -lt 9 ]; then export BURDEN_MODE=rate; export TASK=$i
else                   export BURDEN_MODE=count; export TASK=$(( i - 9 )); fi
export POOL_SET=all6 INTER_SET=all
export BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid CRTAG=""
echo "[$(date)] host=$(hostname) all6 + all interactions MODE=$BURDEN_MODE TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R

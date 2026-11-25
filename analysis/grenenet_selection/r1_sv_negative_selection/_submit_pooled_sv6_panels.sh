#!/bin/bash
# The FULL model (pooled SVs, no origin filter, reference TE kept, + local-adaptation control) fitted on
# ONE sequencing-technology panel at a time: cactus assembly (80 founders) and short-read (151).
# SUBSET= re-z-scores the burdens within the panel, subsets the kinship matrix, and drops `cactus` from
# the fixed effects (constant within a panel). Outputs: _x_obio1cobio1cbio1c_<panel>.
#SBATCH --job-name=pooled_sv6_panel
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-11
#SBATCH --output=logs/pooled_sv6_panel_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
# array index 0-5 -> cactus tasks 0-5;  6-11 -> shortread tasks 0-5
if [ "$SLURM_ARRAY_TASK_ID" -lt 6 ]; then
  export SUBSET=cactus;    export TASK=$SLURM_ARRAY_TASK_ID
else
  export SUBSET=shortread; export TASK=$(( SLURM_ARRAY_TASK_ID - 6 ))
fi
export POOL_SET=sv6
export BURDEN_ORIGIN=all
export ONESTAGE_DIR=onestage_all
export FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
echo "[$(date)] host=$(hostname) SUBSET=$SUBSET TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R

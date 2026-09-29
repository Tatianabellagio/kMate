#!/bin/bash
# The local-adaptation control on the REPORTED burden model: split ins/del, `shared` variants, no
# reference-TE (SV_SET=sv4 + ONESTAGE_DIR=onestage = final_models_sv4_noTE, the arm behind
# sv_burden_final_models.ipynb), with founder origin climate and its garden-climate interaction added.
# Writes alongside the existing fits with the _x_obio1cobio1cbio1c suffix -- nothing is overwritten.
#SBATCH --job-name=burden_origctrl
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-6
#SBATCH --output=logs/origctrl_split_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export SV_SET=sv4
export ONESTAGE_DIR=onestage
export FINAL_OUT=final_models_sv4_noTE
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export TASK=$SLURM_ARRAY_TASK_ID
echo "[$(date)] host=$(hostname) task=$TASK  SV_SET=$SV_SET ONESTAGE_DIR=$ONESTAGE_DIR EXTRA_FIXED='$EXTRA_FIXED'"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_final_models.R

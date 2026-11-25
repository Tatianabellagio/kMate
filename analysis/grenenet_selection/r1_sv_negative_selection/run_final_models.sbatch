#!/bin/bash
#SBATCH --job-name=burden_final
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --array=0-8
#SBATCH --output=logs/final_models_%A_%a.out
#SBATCH --error=logs/final_models_%A_%a.out
# Final founder burden models with lme4qtl, one fit per task (see _sv_burden_final_models.R):
#   0 Model 1 (REML)  1 Model 2 (REML)  2 Model 2 (ML)  3-8 Model 2 (ML) without one interaction each (6 SV terms incl. TE)
set -eo pipefail
cd /global/home/users/tbellg/kmate/analysis/grenenet_selection/r1_sv_negative_selection
echo "[$(date)] host=$(hostname) task=$SLURM_ARRAY_TASK_ID"
TASK=$SLURM_ARRAY_TASK_ID /global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_final_models.R
echo "[$(date)] done"

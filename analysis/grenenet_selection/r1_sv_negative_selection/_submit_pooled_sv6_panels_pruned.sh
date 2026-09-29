#!/bin/bash
# The FULL model per sequencing-technology panel, with the boundary (singular) random garden slopes pruned.
# From the unpruned panel fits: cactus had sv_CDS slope var = 0 and fid var ~ 6e-7; short-read had
# sv_CDS and sv_TE slope vars = 0. A slope pinned at zero still costs a parameter and leaves that term's
# interaction SE anti-conservative, so it is removed rather than left at the boundary.
# Outputs: _x_obio1cobio1cbio1c_<panel>_re<dropped>.
#SBATCH --job-name=pooled_sv6_pruned
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-11
#SBATCH --output=logs/pooled_sv6_pruned_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
if [ "$SLURM_ARRAY_TASK_ID" -lt 6 ]; then
  export SUBSET=cactus;    export RE_DROP="fid:sv_CDS";   export TASK=$SLURM_ARRAY_TASK_ID
else
  export SUBSET=shortread; export RE_DROP="sv_CDS:sv_TE"; export TASK=$(( SLURM_ARRAY_TASK_ID - 6 ))
fi
export POOL_SET=sv6
export BURDEN_ORIGIN=all
export ONESTAGE_DIR=onestage_all
export FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
echo "[$(date)] host=$(hostname) SUBSET=$SUBSET RE_DROP=$RE_DROP TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R

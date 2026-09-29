#!/bin/bash
# Re-run the three arms that share final_models_pooled_all: the main sv6 fit (whose files were overwritten
# when cds2/all6 were first run before SUF carried POOL_SET), plus cds2 and all6 under their new names.
#SBATCH --job-name=burden_rebuild
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-16
#SBATCH --output=logs/burden_rebuild_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
i=$SLURM_ARRAY_TASK_ID
if   [ "$i" -lt 6 ];  then export POOL_SET=sv6;  export TASK=$i
elif [ "$i" -lt 11 ]; then export POOL_SET=cds2; export TASK=$(( i - 6 ))
else                       export POOL_SET=all6; export TASK=$(( i - 11 ))
fi
export BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export RE_DROP=fid
echo "[$(date)] host=$(hostname) POOL_SET=$POOL_SET TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R

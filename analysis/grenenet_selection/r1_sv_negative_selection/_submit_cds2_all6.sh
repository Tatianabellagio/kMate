#!/bin/bash
# Two burden parameterisations that isolate the coding contrast, on the main model's covariate set
# (no origin filter, local-adaptation control, fid pruned, SNP burden + assembly route + kinship):
#   cds2  sv_CDS + sv_nonCDS, where nonCDS pools UTR + intron + TE + promoter + intergenic
#   all6  every genomic context as its own term -- closes the gap that UTR and intron (20% of passing
#         events, 24% of a founder's carried events) entered NO burden term in the sv4/sv6 arms
#SBATCH --job-name=burden_ctx
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-9
#SBATCH --output=logs/burden_ctx_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
# cds2 has 2 interaction terms -> 5 tasks (0-4); all6 has 3 -> 6 tasks (0-5)
if [ "$SLURM_ARRAY_TASK_ID" -lt 5 ]; then
  export POOL_SET=cds2; export TASK=$SLURM_ARRAY_TASK_ID
else
  export POOL_SET=all6; export TASK=$(( SLURM_ARRAY_TASK_ID - 5 ))
fi
export BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export RE_DROP=fid
echo "[$(date)] host=$(hostname) POOL_SET=$POOL_SET TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R

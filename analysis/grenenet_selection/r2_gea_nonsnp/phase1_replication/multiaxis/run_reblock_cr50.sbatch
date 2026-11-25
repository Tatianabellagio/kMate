#!/bin/bash
#SBATCH --job-name=reblock_cr50
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --array=0-3
#SBATCH --output=/global/scratch/projects/fc_moilab/tbellg/kmate/logs/reblock_cr50_%A_%a.out
#SBATCH --error=/global/scratch/projects/fc_moilab/tbellg/kmate/logs/reblock_cr50_%A_%a.out
# GEA step 5 (_submit_gea.sh): the production tables, tile_cr50_arch3v2/ (reblock_cr50.py). One task per class,
# 22 axes each.
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection
CLS=(snp smallindel sv nonsnp)
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python -u \
    r2_gea_nonsnp/phase1_replication/multiaxis/reblock_cr50.py \
    --classes ${CLS[$SLURM_ARRAY_TASK_ID]}
echo "== ALLDONE_REBLOCK_CR50 ${CLS[$SLURM_ARRAY_TASK_ID]} =="

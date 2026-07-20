#!/bin/bash
#SBATCH --job-name=temp_s1
#SBATCH --partition=savio4_htc
#SBATCH --account=co_moilab
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --array=0-1
#SBATCH --output=temporal_stage1_%A_%a.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
CLS=(snp nonsnp)
C=${CLS[$SLURM_ARRAY_TASK_ID]}
echo "host=$(hostname) class=$C"
$PY temporal_stage1.py --class "$C" \
    --out /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/results/temporal_twostage

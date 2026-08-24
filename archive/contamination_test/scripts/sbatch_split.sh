#!/bin/bash
#SBATCH --job-name=panel_split
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-5
#SBATCH --requeue
#SBATCH --output=logs/split_chr%a.out
#SBATCH --error=logs/split_chr%a.err

mkdir -p logs
set -eo pipefail
CH="${SLURM_ARRAY_TASK_ID}"
/global/scratch/users/tbellg/hapfire_sv/contamination_test/scripts/split_panel_per_chrom.sh "${CH}"

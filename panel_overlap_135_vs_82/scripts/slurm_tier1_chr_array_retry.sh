#!/bin/bash
#SBATCH --job-name=overlap_t1_retry
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --array=4-5
#SBATCH --requeue
#SBATCH --output=logs/tier1_chr%a_retry_%j.out
#SBATCH --error=logs/tier1_chr%a_retry_%j.err
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv
CHROM="Chr${SLURM_ARRAY_TASK_ID}"
bash "$BASE/panel_overlap_135_vs_82/scripts/tier1_compute_ac.sh" "$CHROM"

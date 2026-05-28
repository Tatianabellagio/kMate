#!/bin/bash
#SBATCH --job-name=overlap_t1_ac
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --array=2-5
#SBATCH --requeue
#SBATCH --output=logs/tier1_chr%a_%j.out
#SBATCH --error=logs/tier1_chr%a_%j.err
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv
CHROM="Chr${SLURM_ARRAY_TASK_ID}"

bash "$BASE/panel_overlap_135_vs_82/scripts/tier1_compute_ac.sh" "$CHROM"

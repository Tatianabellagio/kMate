#!/bin/bash
#SBATCH --job-name=overlap_t1_ac
#SBATCH --partition=bse
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --array=2-5
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_135_vs_82/logs/tier1_chr%a_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_135_vs_82/logs/tier1_chr%a_%j.err
set -euo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
CHROM="Chr${SLURM_ARRAY_TASK_ID}"

bash "$BASE/panel_overlap_135_vs_82/scripts/tier1_compute_ac.sh" "$CHROM"

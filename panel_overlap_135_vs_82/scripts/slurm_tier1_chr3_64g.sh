#!/bin/bash
#SBATCH --job-name=overlap_t1_chr3
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=logs/tier1_chr3_retry_%j.out
#SBATCH --error=logs/tier1_chr3_retry_%j.err
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv
bash "$BASE/panel_overlap_135_vs_82/scripts/tier1_compute_ac.sh" Chr3

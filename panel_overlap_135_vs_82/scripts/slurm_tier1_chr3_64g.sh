#!/bin/bash
#SBATCH --job-name=overlap_t1_chr3
#SBATCH --partition=bse
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_135_vs_82/logs/tier1_chr3_retry_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_135_vs_82/logs/tier1_chr3_retry_%j.err
set -euo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
bash "$BASE/panel_overlap_135_vs_82/scripts/tier1_compute_ac.sh" Chr3

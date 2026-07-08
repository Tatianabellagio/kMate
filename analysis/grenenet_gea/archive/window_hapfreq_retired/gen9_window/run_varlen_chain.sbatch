#!/bin/bash
#SBATCH --job-name=varlen_gea
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=0:50:00
#SBATCH --output=logs/varlen_gea_%j.out
#SBATCH --error=logs/varlen_gea_%j.err
# Variable-length trajectories: Step 1+2 -> GEA (Step 3+4) -> block-WZA, all on the
# pipelineB_varlen dir, so the 19-site persistent results stay intact for comparison.
set -euo pipefail
hostname
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
export N_PERM=${N_PERM:-5000}
export TRAJ_MODE=varlen
export PB_DIR=results/grenenet_gea/hapfreq/pipelineB_varlen
echo "=== Step 1+2 (varlen trajectories) ==="
python analysis/grenenet_gea/build_hap_trajectories.py
echo "=== Step 3+4 (per-block GEA) N_PERM=$N_PERM ==="
python analysis/grenenet_gea/build_hap_gea.py
echo "=== block-WZA N_PERM=$N_PERM ==="
python analysis/grenenet_gea/build_hap_wza.py
echo "DONE"

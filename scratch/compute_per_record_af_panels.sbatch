#!/bin/bash
#SBATCH --job-name=af_panels
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=30:00
#SBATCH --array=1-3
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/kmate/scratch/logs/af_panels_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/kmate/scratch/logs/af_panels_%A_%a.err

set -euo pipefail
PYBIN=/home/tbellagio/miniforge3/envs/hapfm/bin/python
ROOT=/carnegie/nobackup/scratch/tbellagio/kmate
mkdir -p $ROOT/scratch/logs
echo "=== task ${SLURM_ARRAY_TASK_ID} on $SLURMD_NODENAME, $(date) ==="
$PYBIN -u $ROOT/scratch/compute_per_record_af_panels.py $SLURM_ARRAY_TASK_ID
echo "=== task ${SLURM_ARRAY_TASK_ID} done $(date) ==="

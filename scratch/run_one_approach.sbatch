#!/bin/bash
#SBATCH --job-name=approach
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=30:00
#SBATCH --array=1-12
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/kmate/scratch/logs/approach_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/kmate/scratch/logs/approach_%A_%a.err

set -euo pipefail
PYBIN=/home/tbellagio/miniforge3/envs/hapfm/bin/python
ROOT=/carnegie/nobackup/scratch/tbellagio/kmate
mkdir -p $ROOT/scratch/logs $ROOT/scratch/three_approaches
echo "=== task ${SLURM_ARRAY_TASK_ID} on $SLURMD_NODENAME, $(date) ==="
$PYBIN -u $ROOT/scratch/run_one_approach.py $SLURM_ARRAY_TASK_ID
echo "=== task ${SLURM_ARRAY_TASK_ID} done $(date) ==="

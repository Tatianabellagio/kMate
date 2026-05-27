#!/bin/bash
#SBATCH --job-name=uniform_anchor
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=20:00
#SBATCH --array=1-6
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/logs/uniform_anchor_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/logs/uniform_anchor_%A_%a.err

set -euo pipefail
PYBIN=/home/tbellagio/miniforge3/envs/hapfm/bin/python
ROOT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
mkdir -p $ROOT/scratch/logs
echo "=== task ${SLURM_ARRAY_TASK_ID} on $SLURMD_NODENAME, $(date) ==="
$PYBIN -u $ROOT/scratch/run_uniform_anchor.py $SLURM_ARRAY_TASK_ID
echo "=== task ${SLURM_ARRAY_TASK_ID} done $(date) ==="

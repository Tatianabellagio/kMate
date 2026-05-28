#!/bin/bash
#SBATCH --job-name=three_approaches
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/kmate/scratch/three_approaches.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/kmate/scratch/three_approaches.err

set -euo pipefail
PYBIN=/home/tbellagio/miniforge3/envs/hapfm/bin/python
ROOT=/carnegie/nobackup/scratch/tbellagio/kmate

echo "=== node $SLURMD_NODENAME, started $(date) ==="
$PYBIN -u $ROOT/scratch/three_approaches.py
echo "=== finished $(date) ==="

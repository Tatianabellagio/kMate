#!/bin/bash
#SBATCH --job-name=atom_mech
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/logs/atom_mech_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/logs/atom_mech_%j.err

set -euo pipefail
PYBIN=/home/tbellagio/miniforge3/envs/hapfm/bin/python
ROOT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv

echo "=== atomization-mechanism test on $SLURMD_NODENAME, $(date) ==="
$PYBIN -u $ROOT/scratch/test_atomization_mechanism.py
echo "=== done $(date) ==="

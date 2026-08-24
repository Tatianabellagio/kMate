#!/bin/bash
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --mem=56G
#SBATCH --time=6:00:00
# Usage: sbatch --job-name=NAME -o LOG run_one.sh <script_relpath> [args...]
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
export MEMB_TAG=clq90
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "HOST $(hostname)  START $(date)  SCRIPT $*"
$PY "$@"
echo "DONE $(date) EXIT $?"

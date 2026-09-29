#!/bin/bash
#SBATCH --job-name=cargoctrl
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --requeue
#SBATCH --output=/global/home/users/tbellg/kmate/analysis/grenenet_selection/genes/convergence/results/_cargo_null_control_%j.out
#SBATCH --error=/global/home/users/tbellg/kmate/analysis/grenenet_selection/genes/convergence/results/_cargo_null_control_%j.out
# Length-matched signal-free control: is the dinucleotide-shuffle null calibrated at length?
set -uo pipefail
D=/global/home/users/tbellg/kmate/analysis/grenenet_selection/genes/convergence
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
cd $D && $PY -u cargo_null_control.py 32

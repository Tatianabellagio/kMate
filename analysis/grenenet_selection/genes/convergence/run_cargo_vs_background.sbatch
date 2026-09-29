#!/bin/bash
#SBATCH --job-name=cargobg
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --requeue
#SBATCH --output=/global/home/users/tbellg/kmate/analysis/grenenet_selection/genes/convergence/results/_cargo_vs_background_%j.out
#SBATCH --error=/global/home/users/tbellg/kmate/analysis/grenenet_selection/genes/convergence/results/_cargo_vs_background_%j.out
# Candidate cargo vs length/class/frequency-matched NON-HIT panel insertions, presence-based.
# Both arms are real insertion sequence, so the dinucleotide-shuffle length artefact cannot arise.
set -uo pipefail
D=/global/home/users/tbellg/kmate/analysis/grenenet_selection/genes/convergence
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
cd $D && $PY -u cargo_vs_background.py 32

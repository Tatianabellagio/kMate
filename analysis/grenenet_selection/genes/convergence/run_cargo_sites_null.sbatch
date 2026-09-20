#!/bin/bash
#SBATCH --job-name=cargosites
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --requeue
#SBATCH --output=/global/home/users/tbellg/kmate/analysis/grenenet_selection/genes/convergence/results/_cargo_sites_null_%j.out
#SBATCH --error=/global/home/users/tbellg/kmate/analysis/grenenet_selection/genes/convergence/results/_cargo_sites_null_%j.out
# Cargo enrichment recomputed on merged SITES rather than motif hits (audit 2026-09-20).
set -uo pipefail
D=/global/home/users/tbellg/kmate/analysis/grenenet_selection/genes/convergence
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
cd $D && $PY -u cargo_sites_null.py 16

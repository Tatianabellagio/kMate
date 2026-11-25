#!/bin/bash
#SBATCH --job-name=evmat_cr50
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/genes/convergence/event_gea/cr50/_build_%j.out
#SBATCH --error=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/genes/convergence/event_gea/cr50/_build_%j.out
# GEA step 3 (_submit_gea.sh): event-arm units, MAF >= 0.05 and founder call rate >= 0.50 (build_event_matrix.py
# defaults) -> event_gea/cr50_arch3v2/
set -euo pipefail
D=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/genes/convergence/event_gea
V=${KMATE_PANEL:-arch3v2}; S=""; [ "$V" = arch3 ] || S="_$V"     # panel version (lib.PANEL): outputs under *_arch3v2 dirs
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python -u $D/build_event_matrix.py --out $D/cr50$S
echo "EVMAT_CR50_DONE"

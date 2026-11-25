#!/bin/bash
#SBATCH --job-name=nb_gea50
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --requeue
#SBATCH --output=/global/scratch/projects/fc_moilab/tbellg/kmate/logs/nb_gea50_%A_%a.out
#SBATCH --error=/global/scratch/projects/fc_moilab/tbellg/kmate/logs/nb_gea50_%A_%a.out
# Rebuild the GEA notebooks on the 2026-09-30 re-run (MAF 0.05 + founder call rate 0.50,
# SV/non-SNP on truvari events; source results/multiaxis/tile_cr50).
#   0  raw_manhattan_snp_vs_{nonsnp,sv}_lfmm_tile.ipynb
#   1  peak_alignment_manhattan_gea.ipynb
#   2  class_signal_gain.ipynb          (needs gea_prep50 task 0 and 1)
#   3  gea_hit_landscape.ipynb          (needs gea_prep50 task 2)
# The pre-filter versions of all four are archived in notebooks/old_lfmm/.
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
case $SLURM_ARRAY_TASK_ID in
  0) $PY -u r2_gea_nonsnp/phase1_replication/multiaxis/_build_raw_manhattan_lfmm_tile_nb.py ;;
  1) $PY -u genes/convergence/_build_peak_alignment_nb.py pc1 ;;
  2) $PY -u genes/convergence/_build_class_signal_nb.py ;;
  3) $PY -u genes/convergence/_build_gea_hit_landscape_nb.py ;;
esac
echo "== ALLDONE_NB_GEA50 $SLURM_ARRAY_TASK_ID =="

#!/bin/bash
#SBATCH --job-name=c3_kendall
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/logs/c3_kendall_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/logs/c3_kendall_%j.out
# Consolidated 3-class GEA: fresh Kendall-tau per-record on the current-AF (Jul-8)
# class matrices, for snp / sv / smallindel, gen9 vs bio1. kendall_snp/smallindel
# were stale (Jun-16, pre-clq90-AF); this re-runs all three on the same matrices.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_gea/phase1_replication
echo "== host $(hostname) =="
for CLS in snp sv smallindel; do
  echo "===== kendall $CLS gen9 bio1 ====="
  $PY -u $PR/run_kendall.py --class "$CLS" --gen 9 --climate bio1 --threads "${SLURM_CPUS_PER_TASK:-8}"
done
echo "== ALLDONE_C3_KENDALL =="

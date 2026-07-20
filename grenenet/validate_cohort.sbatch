#!/bin/bash
#SBATCH --job-name=validate_cohort
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:30:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/logs/reconcile/validate_cohort_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/logs/reconcile/validate_cohort_%j.out
# Post-cohort completeness check (guards Finding 1: preemption-truncated per-chrom
# TSVs). --fix removes any short final + per-chrom TSV and writes a redo-manifest.
mkdir -p /global/scratch/users/tbellg/kmate/logs/reconcile
source ~/.bashrc 2>/dev/null || true
mamba activate kmate 2>/dev/null || true
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
echo "host=$(hostname)  $(date)"
python grenenet/validate_cohort.py \
  --manifest data/sample_manifest_usesample.tsv \
  --out-dir analysis/grenenet_gea/rerun_kfw_hb/evolved \
  --fix
echo "[$(date)] DONE"

#!/bin/bash
#SBATCH --job-name=clq90_binom
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/logs/clq90_binom_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/logs/clq90_binom_%j.out
# clq0.9 replication: effective-N (site-deflated) binomial GEA, last-gen (gen9),
# for the two variant classes SNP and non-SNP (SV+small-indel together). The
# per-record binomial is block-independent, so this is the ONLY model that
# re-runs vs the phase-1-block outputs (kendall/lfmm are reused + re-blocked).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python   # statsmodels in basic
PR=analysis/grenenet_gea/phase1_replication
OUT=analysis/grenenet_gea/phase1_replication/results/clq90/binomial

echo "== host $(hostname) | cpus=${SLURM_CPUS_PER_TASK:-16} =="
for CLS in snp nonsnp; do
  echo "===== binomial effective-N=site: $CLS gen9 bio1 ====="
  $PY -u $PR/run_binomial.py --class "$CLS" --gen 9 --climate bio1 \
      --effective-n site --threads "${SLURM_CPUS_PER_TASK:-16}" --out "$OUT"
done
echo "== ALLDONE_CLQ90_BINOM =="

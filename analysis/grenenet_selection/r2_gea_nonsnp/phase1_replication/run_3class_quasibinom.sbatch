#!/bin/bash
#SBATCH --job-name=c3_qbinom
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/logs/c3_qbinom_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/logs/c3_qbinom_%j.out
# Consolidated 3-class GEA: PRODUCTION binomial = per-variant quasi-binomial with
# K=16 LFMM latent factors (STATUS_clq90.md §0B; pval_quasi). For snp / sv /
# smallindel, gen9 vs bio1. statsmodels lives in the `basic` env. After each class
# we emit a phase-1-shaped binomial_<cls>_gen9_bio1.csv (pval = pval_quasi) so
# reblock.py / run_wza.py consume it exactly like kendall/lfmm.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
PR=analysis/grenenet_selection/phase1_replication
QDIR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/clq90/quasibinom
BDIR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/clq90/binomial
mkdir -p "$QDIR" "$BDIR"
K=16
echo "== host $(hostname) | cpus=${SLURM_CPUS_PER_TASK:-16} =="
for CLS in snp sv smallindel; do
  echo "===== quasibinom K=$K: $CLS gen9 bio1 ====="
  $PY -u $PR/run_quasibinom.py --class "$CLS" --gen 9 --climate bio1 --k $K \
      --threads "${SLURM_CPUS_PER_TASK:-16}" --out "$QDIR"
  echo "== emit phase-1-shaped binomial_${CLS}_gen9_bio1.csv (pval=pval_quasi) =="
  $PY - "$CLS" "$K" <<'EOF'
import sys, pandas as pd
cls, k = sys.argv[1], sys.argv[2]
QDIR="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/clq90/quasibinom"
BDIR="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/clq90/binomial"
q=pd.read_csv(f"{QDIR}/quasibinom_lf{k}_{cls}_gen9_bio1.csv")
q=q.rename(columns={"pval_quasi":"pval"})
out=f"{BDIR}/binomial_{cls}_gen9_bio1.csv"
q[["chrom","pos","ref_len","alt_len","MAF","block","slope","pval"]].to_csv(out,index=False)
print(f"  wrote {out} | records={len(q):,} min pval={q.pval.min():.2e}")
EOF
done
echo "== ALLDONE_C3_QBINOM =="

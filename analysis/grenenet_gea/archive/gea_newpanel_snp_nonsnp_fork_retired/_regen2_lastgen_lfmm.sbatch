#!/bin/bash
#SBATCH --job-name=lastgen_lfmm_re
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=10:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/lastgen_lfmm_%j.log
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/lastgen_lfmm_%j.log
# REGEN of the per-plot last-gen LFMM (run_lfmm_lastgen.R) for snp + nonsnp on the
# FRESH 352-pool class matrices. Mirrors phase1_replication/run_lfmm_gea.sh's
# build_lfmm_input.py -> run_lfmm_lastgen.R -> join steps (no phase1 WZA; the clq0.9
# WZA is done by run_wza_clq09.py). Writes lfmm_{cls}_gen9_bio1.csv that run_wza_clq09
# consumes.
set -euo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_gea/phase1_replication
LDIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/results/lfmm
CMDIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/results/class_matrices
K=16

echo "== host $(hostname) =="
for CLS in nonsnp snp; do
  echo "===== $CLS ====="
  echo "== build LFMM input (Delta p, impute, env) =="
  $PY -u $PR/build_lfmm_input.py --class $CLS --gen 9 --climate bio1

  echo "== LFMM ridge K=$K + gif =="
  $R $PR/run_lfmm_lastgen.R "$LDIR/lfmm_${CLS}_gen9" $K "$LDIR/lfmm_${CLS}_gen9_bio1_calibp.csv"

  echo "== join calibrated p -> per-record table (WZA input) =="
  $PY - "$CLS" <<'EOF'
import sys, pandas as pd
cls=sys.argv[1]
CMDIR="analysis/grenenet_gea/phase1_replication/results/class_matrices"
LDIR="analysis/grenenet_gea/phase1_replication/results/lfmm"
rec=pd.read_csv(f"{CMDIR}/{cls}_gen9.records.csv")
pv=pd.read_csv(f"{LDIR}/lfmm_{cls}_gen9_bio1_calibp.csv")["pval"].to_numpy()
assert len(pv)==len(rec),(len(pv),len(rec))
rec=rec.assign(pval=pv).rename(columns={"maf":"MAF"})
out=f"{LDIR}/lfmm_{cls}_gen9_bio1.csv"
rec[["chrom","pos","ref_len","alt_len","MAF","block","pval"]].to_csv(out,index=False)
print("  wrote",out,"| n=%d min p=%.2e"%(len(rec),rec.pval.min()))
EOF
done
echo "== ALLDONE_LASTGEN_LFMM_RE =="

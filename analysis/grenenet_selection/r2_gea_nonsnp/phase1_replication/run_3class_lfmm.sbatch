#!/bin/bash
#SBATCH --job-name=c3_lfmm
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/logs/c3_lfmm_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/logs/c3_lfmm_%j.out
# Consolidated 3-class GEA: fresh LFMM-ridge K=16 + gif-calibrated p per record,
# for snp / sv / smallindel, gen9 vs bio1, on the current-AF (Jul-8) class matrices.
# Writes lfmm_<cls>_gen9_bio1.csv (chrom,pos,ref_len,alt_len,MAF,block,pval) for reblock.
# This is the long-pole job (genome-wide K=16 ridge). No WZA here (done in run_3class_wza).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_selection/phase1_replication
LDIR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/lfmm
CMDIR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/class_matrices
mkdir -p "$LDIR"
K=16
echo "== host $(hostname) =="
for CLS in snp sv smallindel; do
  echo "===== LFMM $CLS gen9 bio1 (K=$K) ====="
  echo "== build LFMM input (delta-p, impute, env) =="
  $PY -u $PR/build_lfmm_input.py --class "$CLS" --gen 9 --climate bio1
  echo "== lfmm_ridge K=$K + gif calibrate =="
  $R $PR/run_lfmm_lastgen.R "$LDIR/lfmm_${CLS}_gen9" $K "$LDIR/lfmm_${CLS}_gen9_bio1_calibp.csv"
  echo "== join calibrated p -> per-record WZA input =="
  $PY - "$CLS" <<'EOF'
import sys, pandas as pd
cls=sys.argv[1]
CMDIR="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/class_matrices"
LDIR="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/lfmm"
rec=pd.read_csv(f"{CMDIR}/{cls}_gen9.records.csv")
pv=pd.read_csv(f"{LDIR}/lfmm_{cls}_gen9_bio1_calibp.csv")["pval"].to_numpy()
assert len(pv)==len(rec),(len(pv),len(rec))
rec=rec.assign(pval=pv).rename(columns={"maf":"MAF"})
out=f"{LDIR}/lfmm_{cls}_gen9_bio1.csv"
rec[["chrom","pos","ref_len","alt_len","MAF","block","pval"]].to_csv(out,index=False)
print(f"  wrote {out} | min p={rec.pval.min():.2e}")
EOF
done
echo "== ALLDONE_C3_LFMM =="

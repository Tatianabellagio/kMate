#!/bin/bash
#SBATCH --job-name=lastgen_lfmm
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=10:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/lastgen_lfmm_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/lastgen_lfmm_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication
WZA=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/wza/wza_script.py
LDIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/lfmm
KDIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/kendall
WDIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/wza
CMDIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/class_matrices
mkdir -p "$LDIR" "$WDIR"
K=16

for CLS in snp smallindel sv; do
  echo "===== $CLS ====="
  echo "== build LFMM input (Δp, impute, env) =="
  $PY -u $PR/build_lfmm_input.py --class $CLS --gen 9 --climate bio1

  echo "== LFMM ridge K=$K + gif =="
  $R $PR/run_lfmm_lastgen.R "$LDIR/lfmm_${CLS}_gen9" $K "$LDIR/lfmm_${CLS}_gen9_bio1_calibp.csv"

  echo "== join calibrated p -> WZA input =="
  $PY - "$CLS" <<'EOF'
import sys, pandas as pd
cls=sys.argv[1]
CMDIR="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/class_matrices"
LDIR="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/lfmm"
KDIR="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/kendall"
rec=pd.read_csv(f"{CMDIR}/{cls}_gen9.records.csv")
pv=pd.read_csv(f"{LDIR}/lfmm_{cls}_gen9_bio1_calibp.csv")["pval"].to_numpy()
assert len(pv)==len(rec),(len(pv),len(rec))
rec=rec.assign(pval=pv).rename(columns={"maf":"MAF"})
out=f"{LDIR}/lfmm_{cls}_gen9_bio1.csv"
rec[["chrom","pos","ref_len","alt_len","MAF","block","pval"]].to_csv(out,index=False)
print("  wrote",out,"| min p=%.2e"%rec.pval.min())
EOF

  echo "== WZA deg7 nocap + cap2000 on LFMM p =="
  for spec in "deg7nocap 7 10 0" "deg7cap2000 7 10 2000"; do
    set -- $spec; tag=$1; deg=$2; minE=$3; cap=$4
    d=$(mktemp -d /tmp/wzaLFMM_${CLS}_${tag}.XXXX)
    ( cd "$d" && $PY "$WZA" --correlations "$LDIR/lfmm_${CLS}_gen9_bio1.csv" \
        --summary_stat pval --window block --MAF MAF --maf_filter 0.05 --sep "," \
        --retain chrom pos --poly_deg $deg --min_entries $minE --sample_snps $cap \
        --output "$WDIR/wza_lfmm_${CLS}_gen9_bio1_${tag}.csv" )
    echo "  wrote wza_lfmm_${CLS}_gen9_bio1_${tag}.csv"
  done
done

echo "== CAM5 + cross-class (LFMM->WZA deg7nocap) =="
$PY - <<'EOF'
import pandas as pd, numpy as np, scipy.stats as st, os
WDIR="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/wza"; CAM5="2_1265"
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n)
    q[o]=(p[o]*n)/(np.arange(n)+1); q[o]=np.minimum.accumulate(q[o][::-1])[::-1]; return np.clip(q,0,1)
def gif(p): p=np.clip(p,1e-300,1); return np.median(st.chi2.isf(p,1))/st.chi2.isf(0.5,1)
sig={}
print(f"{'class':11s} | {'blk':>5} {'Bonf':>4} {'BHq<.05':>7} {'GIF':>4} | CAM5 2_1265 (LFMM->WZA deg7)")
for cls in ["snp","smallindel","sv"]:
    f=f"{WDIR}/wza_lfmm_{cls}_gen9_bio1_deg7nocap.csv"
    if not os.path.exists(f): print(f"{cls}: MISSING"); continue
    w=pd.read_csv(f); gc="gene" if "gene" in w.columns else w.columns[0]
    w[gc]=w[gc].astype(str).str.replace(r"\.0$","",regex=True)
    ws=w.dropna(subset=["Z_pVal"]).sort_values("Z_pVal").reset_index(drop=True)
    n=len(ws); p=ws["Z_pVal"].to_numpy(); q=bh(p); r=ws.index[ws[gc]==CAM5]
    sig[cls]=set(ws.loc[q<0.05,gc])
    cam=(f"rank {int(r[0])+1}/{n} p={ws.loc[r[0],'Z_pVal']:.2e} q={q[int(r[0])]:.3f}" if len(r) else "ABSENT")
    print(f"{cls:11s} | {n:5d} {int((p<0.05/n).sum()):4d} {int((q<0.05).sum()):7d} {gif(p):4.2f} | {cam}")
print("\n=== LFMM->WZA BH-sig overlap with Kendall->WZA (snp) ===")
for cls in ["snp","smallindel","sv"]:
    kf=f"{WDIR}/wza_kendall_{cls}_gen9_bio1_deg7nocap.csv"
    if cls in sig and os.path.exists(kf):
        kw=pd.read_csv(kf); gc="gene" if "gene" in kw.columns else kw.columns[0]
        kw[gc]=kw[gc].astype(str).str.replace(r"\.0$","",regex=True)
        kw=kw.dropna(subset=["Z_pVal"]); ksig=set(kw.loc[bh(kw["Z_pVal"].to_numpy())<0.05,gc])
        print(f"  {cls}: LFMM BH-sig={len(sig[cls])}  Kendall BH-sig={len(ksig)}  shared={len(sig[cls]&ksig)}")
EOF
echo "== ALLDONE_LFMM =="

#!/bin/bash
#SBATCH --job-name=lastgen_svindel
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/phase1_replication/lastgen_svindel_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/phase1_replication/lastgen_svindel_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_gea/phase1_replication
WZA=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/wza_script.py
KDIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/phase1_replication/results/kendall
WDIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/phase1_replication/results/wza
mkdir -p "$WDIR"

echo "== host $(hostname) =="
echo "== STEP 1: class split (sv + smallindel, gen9 = last_gen) =="
$PY -u $PR/build_class_matrices.py --classes sv smallindel --gens 9

for CLS in sv smallindel; do
  echo "== STEP 2 ($CLS): Kendall vs bio1 (gen9) =="
  $PY -u $PR/run_kendall.py --class $CLS --gen 9 --climate bio1 --threads 8
  echo "== STEP 3 ($CLS): WZA deg7 nocap + cap2000 =="
  for spec in "deg7nocap 7 10 0" "deg7cap2000 7 10 2000"; do
    set -- $spec; tag=$1; deg=$2; minE=$3; cap=$4
    d=$(mktemp -d /tmp/wzaLG_${CLS}_${tag}.XXXX)
    ( cd "$d" && $PY "$WZA" --correlations "$KDIR/kendall_${CLS}_gen9_bio1.csv" \
        --summary_stat pval --window block --MAF MAF --maf_filter 0.05 --sep "," \
        --retain chrom pos --poly_deg $deg --min_entries $minE --sample_snps $cap \
        --output "$WDIR/wza_kendall_${CLS}_gen9_bio1_${tag}.csv" )
    echo "  wrote wza_kendall_${CLS}_gen9_bio1_${tag}.csv"
  done
done

echo "== STEP 4: cross-class comparison (CAM5 + block overlap, deg7nocap) =="
$PY - <<'EOF'
import pandas as pd, numpy as np, scipy.stats as st, os
WDIR="analysis/grenenet_gea/r2_gea_nonsnp/phase1_replication/results/wza"; CAM5="2_1265"
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n)
    q[o]=(p[o]*n)/(np.arange(n)+1); q[o]=np.minimum.accumulate(q[o][::-1])[::-1]; return np.clip(q,0,1)
def gif(p): p=np.clip(p,1e-300,1); return np.median(st.chi2.isf(p,1))/st.chi2.isf(0.5,1)
sig={}
print(f"{'class':11s} | {'blk':>5} {'Bonf':>4} {'BHq<.05':>7} {'GIF':>4} | CAM5 2_1265 (deg7nocap)")
for cls in ["snp","smallindel","sv"]:
    f=f"{WDIR}/wza_kendall_{cls}_gen9_bio1_deg7nocap.csv"
    if not os.path.exists(f): print(f"{cls}: MISSING"); continue
    w=pd.read_csv(f); gc="gene" if "gene" in w.columns else w.columns[0]
    w[gc]=w[gc].astype(str).str.replace(r"\.0$","",regex=True)
    ws=w.dropna(subset=["Z_pVal"]).sort_values("Z_pVal").reset_index(drop=True)
    n=len(ws); p=ws["Z_pVal"].to_numpy(); q=bh(p); r=ws.index[ws[gc]==CAM5]
    sig[cls]=set(ws.loc[q<0.05,gc])
    cam=(f"rank {int(r[0])+1}/{n} p={ws.loc[r[0],'Z_pVal']:.2e} q={q[int(r[0])]:.3f}" if len(r) else "ABSENT")
    print(f"{cls:11s} | {n:5d} {int((p<0.05/n).sum()):4d} {int((q<0.05).sum()):7d} {gif(p):4.2f} | {cam}")
print("\n=== BH-significant block overlap across classes (deg7nocap, last_gen) ===")
ks=[k for k in ["snp","smallindel","sv"] if k in sig]
for i in range(len(ks)):
    for j in range(i+1,len(ks)):
        a,b=ks[i],ks[j]; print(f"  {a} ({len(sig[a])}) ∩ {b} ({len(sig[b])}) = {len(sig[a]&sig[b])} shared")
if len(ks)==3:
    only_sv=sig["sv"]-sig["snp"]-sig["smallindel"]
    only_ind=sig["smallindel"]-sig["snp"]
    print(f"  SV-ONLY (not in snp or smallindel): {len(only_sv)}  {sorted(only_sv)[:15]}")
    print(f"  smallindel-not-snp: {len(only_ind)}  {sorted(only_ind)[:15]}")
EOF
echo "== ALLDONE_SVINDEL =="

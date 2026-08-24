#!/bin/bash
#SBATCH --job-name=lastgen_gea
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/lastgen_gea_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/lastgen_gea_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_selection/phase1_replication
WZA=analysis/grenenet_selection/r2_gea_nonsnp/wza_script.py
KEN=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/kendall/kendall_snp_gen9_bio1.csv
WDIR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/wza

echo "== host $(hostname) | mem ${SLURM_MEM_PER_NODE:-?} =="

echo "== STEP 1: class split (snp gen9 = last_gen merge) =="
$PY -u $PR/build_class_matrices.py --classes snp --gens 9

echo "== STEP 2: Kendall vs bio1 (snp gen9) =="
$PY -u $PR/run_kendall.py --class snp --gen 9 --climate bio1 --threads 8

echo "== STEP 3: WZA three regimes on last_gen Kendall =="
run_wza () {  # tag deg minE cap
  d=$(mktemp -d "$WDIR/lastgen_$1.XXXX")
  ( cd "$d" && $PY "/global/scratch/users/tbellg/kmate/$WZA" \
      --correlations "/global/scratch/users/tbellg/kmate/$KEN" \
      --summary_stat pval --window block --MAF MAF --maf_filter 0.05 --sep "," \
      --retain chrom pos --poly_deg $2 --min_entries $3 --sample_snps $4 \
      --output "$WDIR/wza_kendall_snp_gen9_bio1_$1.csv" )
  echo "  wrote wza_kendall_snp_gen9_bio1_$1.csv"
}
run_wza deg7nocap   7 10 0
run_wza deg7cap2000 7 10 2000
run_wza deg2nocap   2 40 0

echo "== STEP 4: CAM5 (block 2_1265) report =="
$PY - <<'EOF'
import pandas as pd, numpy as np, scipy.stats as st, glob, os
WDIR="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/wza"; CAM5="2_1265"
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n)
    q[o]=(p[o]*n)/(np.arange(n)+1); q[o]=np.minimum.accumulate(q[o][::-1])[::-1]; return np.clip(q,0,1)
def gif(p): p=np.clip(p,1e-300,1); return np.median(st.chi2.isf(p,1))/st.chi2.isf(0.5,1)
print("phase-1 last_gen CAM5 (published): Z=10.84  Z_pVal=1.03e-07 (~3e-8 in STATUS)")
print(f"{'regime':12s} | {'blk':>5} {'p<1e-12':>7} {'Bonf':>4} {'BHq<.05':>7} {'GIF':>4} | CAM5 2_1265")
for tag in ["deg7nocap","deg7cap2000","deg2nocap"]:
    f=f"{WDIR}/wza_kendall_snp_gen9_bio1_{tag}.csv"
    if not os.path.exists(f): print(tag,"MISSING"); continue
    w=pd.read_csv(f); gc="gene" if "gene" in w.columns else w.columns[0]
    w[gc]=w[gc].astype(str).str.replace(r"\.0$","",regex=True)
    ws=w.dropna(subset=["Z_pVal"]).sort_values("Z_pVal").reset_index(drop=True)
    n=len(ws); p=ws["Z_pVal"].to_numpy(); q=bh(p)
    r=ws.index[ws[gc]==CAM5]
    cam=(f"rank {int(r[0])+1}/{n} p={ws.loc[r[0],'Z_pVal']:.2e} q={q[int(r[0])]:.3f} Z={ws.loc[r[0],'Z']:.2f}" if len(r) else "ABSENT")
    print(f"{tag:12s} | {n:5d} {int((p<1e-12).sum()):7d} {int((p<0.05/n).sum()):4d} {int((q<0.05).sum()):7d} {gif(p):4.2f} | {cam}")
EOF
echo "== ALLDONE_LASTGEN =="

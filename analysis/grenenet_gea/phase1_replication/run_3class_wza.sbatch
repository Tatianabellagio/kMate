#!/bin/bash
#SBATCH --job-name=c3_wza
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/logs/c3_wza_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/logs/c3_wza_%j.out
# Consolidated 3-class WZA: re-block every per-record output (kendall/lfmm/binomial
# x snp/sv/smallindel) onto the clq0.9 (blocks_mcf90) LD islands, drop gap records,
# then aggregate with WZA in two regimes:
#   deg2   (PRIMARY)  = canonical Booker deg-2, SD-floor hack removed, per-class SNP
#                       cap (1000 snp / 350 sv,smallindel) -> STATUS_clq90 §0 settled.
#   deg7cap2000 (SENS) = phase-1 deg-7 poly, cap 2000 (unaffected by the SD-floor bug).
# Ends with a NaN / negative-SD audit per output (should be 0 in the capped regimes).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_gea/phase1_replication
CLQ=analysis/grenenet_gea/phase1_replication/results/clq90

echo "== host $(hostname) =="
echo "== STEP 1: reblock 3 models x 3 classes onto clq0.9 (blocks_mcf90) =="
$PY -u $PR/reblock.py --models kendall lfmm binomial \
    --classes snp sv smallindel --gen 9 --climate bio1 --r2 0.9

echo "== STEP 2: WZA deg2-capped (primary) + deg7cap2000 (sensitivity) =="
for M in kendall lfmm binomial; do
  for C in snp sv smallindel; do
    IN=$CLQ/wza_in/${M}_${C}_gen9_bio1.csv
    [ -f "$IN" ] || { echo "  SKIP $M $C (no input $IN)"; continue; }
    CAP=350; [ "$C" = "snp" ] && CAP=1000
    echo "  --- $M $C (deg2 cap=$CAP) ---"
    $PY -u $PR/run_wza.py --model $M --class $C --gen 9 --climate bio1 \
        --indir $CLQ/wza_in --out $CLQ/wza --regime deg2 \
        --poly-deg 2 --min-entries 40 --sample-snps $CAP
    echo "  --- $M $C (deg7cap2000) ---"
    $PY -u $PR/run_wza.py --model $M --class $C --gen 9 --climate bio1 \
        --indir $CLQ/wza_in --out $CLQ/wza --regime deg7cap2000 \
        --poly-deg 7 --min-entries 10 --sample-snps 2000
  done
done

echo "== STEP 3: NaN / negative-SD audit per regime =="
$PY - <<'EOF'
import pandas as pd, numpy as np, glob, os
W="analysis/grenenet_gea/phase1_replication/results/clq90/wza"
print(f"{'file':52s} {'blocks':>7} {'NaN_p':>6} {'Bonf':>5} {'BHq<.05':>7}")
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n)
    q[o]=(p[o]*n)/(np.arange(n)+1); q[o]=np.minimum.accumulate(q[o][::-1])[::-1]; return np.clip(q,0,1)
for f in sorted(glob.glob(f"{W}/wza_*_gen9_bio1_*.csv")):
    w=pd.read_csv(f)
    if "Z_pVal" not in w.columns:
        print(f"{os.path.basename(f):52s} (no Z_pVal)"); continue
    p=w["Z_pVal"]; n=len(p); nan=int(p.isna().sum()); ok=w[p.notna()]
    pp=ok["Z_pVal"].to_numpy()
    bonf=int((pp<0.05/max(len(pp),1)).sum()); fdr=int((bh(pp)<0.05).sum()) if len(pp) else 0
    print(f"{os.path.basename(f):52s} {n:7d} {nan:6d} {bonf:5d} {fdr:7d}")
EOF
echo "== ALLDONE_C3_WZA =="

#!/bin/bash
#SBATCH --job-name=clq90_wza
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/logs/clq90_wza_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/logs/clq90_wza_%j.out
# clq0.9 replication WZA: re-block all per-record outputs onto clq0.9 blocks, then
# run WZA at deg-2 canonical (PRIMARY) + deg-7cap2000 (phase-1 sensitivity) for
# every model x class. The clq0.9 blocks are small (SNP max ~5168, only 3 >2000;
# nonsnp max 963), so the deg-2 tail-extrapolation NaN that plagued the coarse
# hapFIRE blocks should be gone -> deg-2 is expected to run clean. This job proves
# that by counting NaN / negative-SD windows per regime.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_gea/phase1_replication
CLQ=results/grenenet_gea/phase1_replication/clq90

echo "== host $(hostname) =="
echo "== STEP 1: re-block all 3 models x 2 classes onto clq0.9 =="
$PY -u $PR/reblock.py --models kendall lfmm binomial --classes snp nonsnp --gen 9 --climate bio1 --r2 0.9

echo "== STEP 2: WZA deg-2 (primary) + deg-7cap2000 (sensitivity) =="
for M in kendall lfmm binomial; do
  for C in snp nonsnp; do
    IN=$CLQ/wza_in/${M}_${C}_gen9_bio1.csv
    [ -f "$IN" ] || { echo "  SKIP $M $C (no input $IN)"; continue; }
    echo "  --- $M $C ---"
    $PY -u $PR/run_wza.py --model $M --class $C --gen 9 --climate bio1 \
        --indir $CLQ/wza_in --out $CLQ/wza --regime deg2
    $PY -u $PR/run_wza.py --model $M --class $C --gen 9 --climate bio1 \
        --indir $CLQ/wza_in --out $CLQ/wza --regime deg7cap2000 \
        --poly-deg 7 --min-entries 10 --sample-snps 2000
  done
done

echo "== STEP 3: NaN / neg-SD audit per regime (does deg-2 run clean on clq0.9?) =="
$PY - <<'EOF'
import pandas as pd, numpy as np, glob, os
W="results/grenenet_gea/phase1_replication/clq90/wza"
print(f"{'file':45s} {'blocks':>7} {'NaN_p':>6} {'Bonf':>5}")
for f in sorted(glob.glob(f"{W}/wza_*_gen9_bio1_*.csv")):
    w=pd.read_csv(f)
    if "Z_pVal" not in w.columns:
        print(f"{os.path.basename(f):45s} (no Z_pVal)"); continue
    p=w["Z_pVal"]; n=len(p); nan=int(p.isna().sum()); ok=w[p.notna()]
    bonf=int((ok["Z_pVal"]<0.05/max(len(ok),1)).sum())
    print(f"{os.path.basename(f):45s} {n:7d} {nan:6d} {bonf:5d}")
EOF
echo "== ALLDONE_CLQ90_WZA =="

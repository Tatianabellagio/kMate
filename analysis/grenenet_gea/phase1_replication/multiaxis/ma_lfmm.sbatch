#!/bin/bash
#SBATCH --job-name=ma_lfmm
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=6:00:00
#SBATCH --array=0-59
#SBATCH --output=/global/scratch/users/tbellg/kmate/logs/ma_lfmm_%A_%a.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/logs/ma_lfmm_%A_%a.out
# Multi-axis LFMM ridge K=16 + gif-calibrated test, gen9. One array task per
# (axis, class). Each task works in its OWN dir (builds the Delta-p Y matrix +
# axis env there), so parallel axes never collide on the un-tagged env filename;
# the big Y.f64 is deleted after R finishes.
# Output -> multiaxis/lfmm/lfmm_{cls}_gen9_{axis}.csv  (chrom,pos,ref_len,alt_len,MAF,block,pval)
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
source analysis/grenenet_gea/phase1_replication/multiaxis/axes.sh
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_gea/phase1_replication
CMDIR=analysis/grenenet_gea/phase1_replication/results/class_matrices
LDIR=analysis/grenenet_gea/phase1_replication/results/multiaxis/lfmm
i=${SLURM_ARRAY_TASK_ID}; AXIS=$(axis_for $i); C=$(cls_for $i)
WORK=$LDIR/work_${C}_${AXIS}; mkdir -p "$WORK"
echo "== host $(hostname) | lfmm $C gen9 $AXIS (work=$WORK) =="

# 1) build Delta-p Y + dims + z-scored axis env into the isolated workdir
$PY -u $PR/build_lfmm_input.py --class "$C" --gen 9 --climate "$AXIS" --out "$WORK"

# 2) LFMM ridge K=16, NO gif calibration -> RAW p (gif logged per axis, not applied)
$R $PR/multiaxis/run_lfmm_nogif.R "$WORK/lfmm_${C}_gen9" 16 \
    "$LDIR/lfmm_${C}_gen9_${AXIS}_calibp.csv" "$LDIR/lfmm_${C}_gen9_${AXIS}_gif.txt"

# 3) join calibrated p onto the record metadata -> WZA-input schema
$PY - "$C" "$AXIS" <<'EOF'
import sys, pandas as pd
cls, axis = sys.argv[1], sys.argv[2]
CMDIR="analysis/grenenet_gea/phase1_replication/results/class_matrices"
LDIR="analysis/grenenet_gea/phase1_replication/results/multiaxis/lfmm"
rec=pd.read_csv(f"{CMDIR}/{cls}_gen9.records.csv")
pv=pd.read_csv(f"{LDIR}/lfmm_{cls}_gen9_{axis}_calibp.csv")["pval"].to_numpy()
assert len(pv)==len(rec),(len(pv),len(rec))
rec=rec.assign(pval=pv).rename(columns={"maf":"MAF"})
out=f"{LDIR}/lfmm_{cls}_gen9_{axis}.csv"
rec[["chrom","pos","ref_len","alt_len","MAF","block","pval"]].to_csv(out,index=False)
print("  wrote",out,"| min p=%.2e"%rec.pval.min())
EOF

# 4) drop the multi-GB Y matrix (climate-independent, rebuildable) but keep the workdir env for provenance
rm -f "$WORK/lfmm_${C}_gen9_Y.f64"
echo "== ALLDONE_MA_LFMM $C $AXIS =="

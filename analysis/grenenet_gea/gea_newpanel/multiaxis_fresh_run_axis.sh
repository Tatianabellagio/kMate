#!/bin/bash
# Per-axis fresh GEA regen: binomial + kendall + lfmm (K=16) x snp/nonsnp/sv, then WZA.
# Usage: multiaxis_fresh_run_axis.sh <AXIS>   (e.g. bio5)
# LFMM intermediates are isolated per-axis (build_lfmm_input writes fixed lfmm_{cls}_gen9_*
# names with no axis tag, so parallel axes would clobber a shared dir).
set -uo pipefail
AXIS="$1"
ROOT=/global/scratch/users/tbellg/kmate
cd "$ROOT"
BPY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
PYK=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_gea/phase1_replication
NP=analysis/grenenet_gea/gea_newpanel
CMDIR=$ROOT/analysis/grenenet_gea/phase1_replication/results/class_matrices
MA=$ROOT/analysis/grenenet_gea/gea_newpanel/results/multiaxis_fresh
LTMP=$MA/_lfmm_tmp/$AXIS
mkdir -p "$MA/binomial" "$MA/kendall" "$MA/lfmm" "$LTMP"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

echo "==== AXIS $AXIS on $(hostname) ===="
for CLS in snp nonsnp sv; do
  echo "-- [$AXIS $CLS] binomial --"
  $BPY -u $PR/run_binomial.py --class $CLS --gen 9 --climate $AXIS --out $MA/binomial --threads 8 \
       || echo "FAIL binomial $AXIS $CLS"
  echo "-- [$AXIS $CLS] kendall --"
  $PYK -u $PR/run_kendall.py --class $CLS --gen 9 --climate $AXIS --out $MA/kendall --threads 8 \
       || echo "FAIL kendall $AXIS $CLS"
  echo "-- [$AXIS $CLS] lfmm K=16 --"
  if $PYK -u $PR/build_lfmm_input.py --class $CLS --gen 9 --climate $AXIS --out $LTMP; then
    $R $PR/run_lfmm_lastgen.R "$LTMP/lfmm_${CLS}_gen9" 16 "$LTMP/lfmm_${CLS}_gen9_${AXIS}_calibp.csv" \
       && $PYK - "$CLS" "$AXIS" "$CMDIR" "$LTMP" "$MA" <<'EOF'
import sys, pandas as pd
cls, axis, cmdir, ltmp, ma = sys.argv[1:6]
rec = pd.read_csv(f"{cmdir}/{cls}_gen9.records.csv")
pv = pd.read_csv(f"{ltmp}/lfmm_{cls}_gen9_{axis}_calibp.csv")["pval"].to_numpy()
assert len(pv) == len(rec), (len(pv), len(rec))
rec = rec.assign(pval=pv).rename(columns={"maf": "MAF"})
out = f"{ma}/lfmm/lfmm_{cls}_gen9_{axis}.csv"
rec[["chrom", "pos", "ref_len", "alt_len", "MAF", "block", "pval"]].to_csv(out, index=False)
print("  wrote", out, "min p=%.2e" % rec.pval.min())
EOF
  else
    echo "FAIL lfmm-input $AXIS $CLS"
  fi
done

echo "-- [$AXIS] WZA (3 models x 3 classes) --"
$PYK -u $NP/run_wza_multiaxis.py --axis $AXIS || echo "FAIL wza $AXIS"
rm -rf "$LTMP"
echo "==== DONE AXIS $AXIS ===="

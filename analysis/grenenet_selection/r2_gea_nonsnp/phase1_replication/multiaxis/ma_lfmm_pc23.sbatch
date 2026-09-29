#!/bin/bash
#SBATCH --job-name=ma_lfmm_pc23
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=6:00:00
#SBATCH --array=0-7
#SBATCH --output=/global/scratch/users/tbellg/kmate/logs/ma_lfmm_pc23_%A_%a.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/logs/ma_lfmm_pc23_%A_%a.out
# Extension of ma_lfmm.sbatch to the 2nd and 3rd climate PCs (2026-08-25).
#
# The 20-axis production set covers bio1-19 + pc1 (46.8% of bioclim variance). pc2
# (23.4%) and pc3 (11.7%) lift that to 82.0% cumulative, and are orthogonal to pc1
# by construction, so they test climate structure the existing axes cannot express.
# Sign convention (see build_class_matrices.py): pc1 -> +bio1, pc2 -> +bio7
# (continentality), pc3 -> +bio10 (warmest-quarter temp). Sign does not change the
# two-sided p, only the direction an effect is read in.
#
# 2 axes x 4 classes = 8 tasks. Array index i: axis = AXES[i/4], cls = CLS[i%4].
# NOTE the pooled `nonsnp` class is included here (ma_lfmm.sbatch runs the 3-class
# split and ma_nonsnp_lfmm.sbatch runs nonsnp separately); wza_in_clq09_tile carries
# all four, so all four are needed for the axis set to stay rectangular.
#
# Output -> results/multiaxis/lfmm/lfmm_{cls}_gen9_{axis}.csv
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication
LDIR=$PR/results/multiaxis/lfmm

AXES=(pc2 pc3)
CLS=(snp sv smallindel nonsnp)
i=${SLURM_ARRAY_TASK_ID}
AXIS=${AXES[$(( i / 4 ))]}
C=${CLS[$(( i % 4 ))]}
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
CMDIR="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/class_matrices"
LDIR="analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/multiaxis/lfmm"
rec=pd.read_csv(f"{CMDIR}/{cls}_gen9.records.csv")
pv=pd.read_csv(f"{LDIR}/lfmm_{cls}_gen9_{axis}_calibp.csv")["pval"].to_numpy()
assert len(pv)==len(rec),(len(pv),len(rec))
rec=rec.assign(pval=pv).rename(columns={"maf":"MAF"})
out=f"{LDIR}/lfmm_{cls}_gen9_{axis}.csv"
rec[["chrom","pos","ref_len","alt_len","MAF","block","pval"]].to_csv(out,index=False)
print("  wrote",out,"| min p=%.2e"%rec.pval.min())
EOF

# 4) drop the multi-GB Y matrix (climate-independent, rebuildable), keep env for provenance
rm -f "$WORK/lfmm_${C}_gen9_Y.f64"
echo "== ALLDONE_MA_LFMM $C $AXIS =="

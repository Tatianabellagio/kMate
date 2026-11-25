#!/bin/bash
#SBATCH --job-name=ma_lfmm_cr50
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=8:00:00
#SBATCH --array=0-87%12
#SBATCH --requeue
#SBATCH --output=/global/scratch/projects/fc_moilab/tbellg/kmate/logs/ma_lfmm_cr50_%A_%a.out
#SBATCH --error=/global/scratch/projects/fc_moilab/tbellg/kmate/logs/ma_lfmm_cr50_%A_%a.out
# Multi-axis LFMM ridge K=16, gen9, WITH the founder call-rate floor of 0.50.
#
# Identical to the pre-floor run (ma_lfmm.sbatch, archived to archive_gea/ on 2026-10-20) in every
# respect except the input class matrices (class_matrices_cr50, built by run_class_matrices_cr50.sbatch)
# and the output directory. 4 classes x 22 axes = 88 tasks, throttled to 12 concurrent because each
# holds ~180 GB. Production uses the snp and smallindel tables from here; sv and nonsnp come from the
# event arm (genes/convergence/event_gea/run_event_lfmm_cr50.sbatch), see reblock_cr50.py.
#
# Why the floor: records called in few founders are a different estimand, not a noisier one
# -- kMate's MAR projection returns the frequency WITHIN the called founders, which varies
# less across gardens than the population frequency, so the same slope comes with less
# residual scatter and draws an inflated threshold. Measured within matched MAF strata the
# 78-only stratum runs lambda 2.15-2.27 vs 1.67-1.80, ~1.7-2x the hit rate, while 0.50-0.90
# is indistinguishable from fully called. See multiaxis/callrate_calibration.py and
# notebooks/flc_case.ipynb.
#
# Out -> results/multiaxis/lfmm_cr50/lfmm_{cls}_gen9_{axis}.csv
#        (chrom,pos,ref_len,alt_len,MAF,block,pval; block = hapFIRE LD block, the clq0.9
#         tiling is re-assigned afterwards by reblock_cr50.py into results/multiaxis/tile_cr50/)
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication
V=${KMATE_PANEL:-arch3v2}; S=""; [ "$V" = arch3 ] || S="_$V"     # panel version (lib.PANEL): outputs under *_arch3v2 dirs
CMDIR=$PR/results/class_matrices_cr50$S
LDIR=$PR/results/multiaxis/lfmm_cr50$S
mkdir -p "$LDIR" /global/scratch/projects/fc_moilab/tbellg/kmate/logs

AXES=(bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 pc1 pc2 pc3)
CLS=(snp sv smallindel nonsnp)
i=${SLURM_ARRAY_TASK_ID}
AXIS=${AXES[$(( i / 4 ))]}
C=${CLS[$(( i % 4 ))]}
WORK=$LDIR/work_${C}_${AXIS}; mkdir -p "$WORK"
echo "== host $(hostname) | lfmm_cr50 $C gen9 $AXIS (work=$WORK) =="

# 1) Delta-p Y + z-scored axis env, from the call-rate-filtered class matrices
$PY -u $PR/build_lfmm_input.py --class "$C" --gen 9 --climate "$AXIS" \
    --cmdir "$CMDIR" --out "$WORK"

# 2) LFMM ridge K=16, raw p (gif recorded per axis, never applied)
$R $PR/multiaxis/run_lfmm_nogif.R "$WORK/lfmm_${C}_gen9" 16 \
    "$LDIR/lfmm_${C}_gen9_${AXIS}_calibp.csv" "$LDIR/lfmm_${C}_gen9_${AXIS}_gif.txt"

# 3) join p onto the record metadata (POSITIONAL -- records.csv is row-aligned with the
#    model's output; never join on chrom/pos, 14.1% of non-SNP records share such a key)
$PY - "$C" "$AXIS" "$CMDIR" "$LDIR" <<'EOF'
import sys, pandas as pd
cls, axis, CMDIR, LDIR = sys.argv[1:5]
rec = pd.read_csv(f"{CMDIR}/{cls}_gen9.records.csv")
pv = pd.read_csv(f"{LDIR}/lfmm_{cls}_gen9_{axis}_calibp.csv")["pval"].to_numpy()
assert len(pv) == len(rec), (len(pv), len(rec))
rec = rec.assign(pval=pv).rename(columns={"maf": "MAF"})
out = f"{LDIR}/lfmm_{cls}_gen9_{axis}.csv"
rec[["chrom", "pos", "ref_len", "alt_len", "MAF", "block", "pval"]].to_csv(out, index=False)
print("  wrote", out, "| n=%d | min p=%.2e" % (len(rec), rec.pval.min()))
EOF

rm -f "$WORK/lfmm_${C}_gen9_Y.f64"
echo "== ALLDONE_MA_LFMM_CR50 $C $AXIS =="

#!/bin/bash
#SBATCH --job-name=ma_lfmm_cr50
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=8:00:00
#SBATCH --array=0-43%12
#SBATCH --requeue
#SBATCH --output=/global/scratch/projects/fc_moilab/tbellg/kmate/logs/ma_lfmm_cr50_%A_%a.out
#SBATCH --error=/global/scratch/projects/fc_moilab/tbellg/kmate/logs/ma_lfmm_cr50_%A_%a.out
# GEA record arm (step 2 of _submit_gea.sh): LFMM ridge K = 16, raw p, for snp and smallindel x 22 climate axes
# (44 tasks, 12 at a time: a SNP task holds ~180 GB). SVs and the pooled non-SNP scan are the event arm
# (genes/convergence/event_gea/run_event_lfmm_cr50.sbatch).
# In:  results/class_matrices_cr50_arch3v2/   (step 1, run_class_matrices_cr50.sbatch)
# Out: results/multiaxis/lfmm_cr50_arch3v2/lfmm_{cls}_gen9_{axis}.csv  (chrom,pos,ref_len,alt_len,MAF,block,pval)
#      and lfmm_{cls}_gen9_{axis}_gif.txt (lambda, recorded only)
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication
V=${KMATE_PANEL:-arch3v2}; S=""; [ "$V" = arch3 ] || S="_$V"     # panel version (lib.PANEL): outputs under *_arch3v2 dirs
CMDIR=$PR/results/class_matrices_cr50$S
LDIR=$PR/results/multiaxis/lfmm_cr50$S
mkdir -p "$LDIR"

AXES=(bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 pc1 pc2 pc3)
CLS=(snp smallindel)
AXIS=${AXES[$(( SLURM_ARRAY_TASK_ID / 2 ))]}; C=${CLS[$(( SLURM_ARRAY_TASK_ID % 2 ))]}
WORK=$LDIR/work_${C}_${AXIS}; mkdir -p "$WORK"
echo "== host $(hostname) | record arm $C $AXIS =="

$PY -u $PR/build_lfmm_input.py --class "$C" --climate "$AXIS" --cmdir "$CMDIR" --out "$WORK"
$R $PR/multiaxis/run_lfmm_nogif.R "$WORK/lfmm_${C}_gen9" 16 \
    "$LDIR/lfmm_${C}_gen9_${AXIS}_calibp.csv" "$LDIR/lfmm_${C}_gen9_${AXIS}_gif.txt"

# p back onto the records by POSITION (row i of records.csv is column i of Y; chrom:pos is not a unique key)
$PY - "$C" "$AXIS" "$CMDIR" "$LDIR" <<'EOF'
import sys, pandas as pd
cls, axis, CMDIR, LDIR = sys.argv[1:5]
rec = pd.read_csv(f"{CMDIR}/{cls}_gen9.records.csv")
pv = pd.read_csv(f"{LDIR}/lfmm_{cls}_gen9_{axis}_calibp.csv")["pval"].to_numpy()
assert len(pv) == len(rec), (len(pv), len(rec))
rec = rec.assign(pval=pv).rename(columns={"maf": "MAF"})
out = f"{LDIR}/lfmm_{cls}_gen9_{axis}.csv"
rec[["chrom", "pos", "ref_len", "alt_len", "MAF", "block", "pval"]].to_csv(out, index=False)
print(f"  wrote {out} | n={len(rec):,} | min p={rec.pval.min():.2e}")
EOF

rm -f "$WORK/lfmm_${C}_gen9_Y.f64"
echo "== ALLDONE_MA_LFMM_CR50 $C $AXIS =="

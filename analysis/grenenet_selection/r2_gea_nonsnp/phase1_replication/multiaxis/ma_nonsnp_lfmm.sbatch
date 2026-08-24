#!/bin/bash
#SBATCH --job-name=man_lfmm
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=6:00:00
#SBATCH --array=0-19
#SBATCH --output=/global/scratch/users/tbellg/kmate/logs/man_lfmm_%A_%a.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/logs/man_lfmm_%A_%a.out
# Pooled non-SNP LFMM ridge K=16 (raw p, no gif — multiaxis convention), gen9,
# one array task per axis, isolated workdir (parallel axes never collide on env).
# LFMM MUST be run on the nonsnp class matrix (joint latent factors) — it is NOT
# the pool of sv+smallindel per-record LFMM. Output schema matches reblock_multiaxis.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
source analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/multiaxis/axes.sh
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication
CMDIR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/class_matrices
LDIR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/multiaxis/lfmm
AXIS=${AXES[$SLURM_ARRAY_TASK_ID]}; C=nonsnp
WORK=$LDIR/work_${C}_${AXIS}; mkdir -p "$WORK"
echo "== host $(hostname) | lfmm nonsnp gen9 $AXIS (work=$WORK) =="

$PY -u $PR/build_lfmm_input.py --class "$C" --gen 9 --climate "$AXIS" --out "$WORK"
$R $PR/multiaxis/run_lfmm_nogif.R "$WORK/lfmm_${C}_gen9" 16 \
    "$LDIR/lfmm_${C}_gen9_${AXIS}_calibp.csv" "$LDIR/lfmm_${C}_gen9_${AXIS}_gif.txt"
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
rm -f "$WORK/lfmm_${C}_gen9_Y.f64"
echo "== ALLDONE_MAN_LFMM nonsnp $AXIS =="

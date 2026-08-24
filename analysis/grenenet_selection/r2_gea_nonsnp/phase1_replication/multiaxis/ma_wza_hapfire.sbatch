#!/bin/bash
#SBATCH --job-name=ma_wza_hapfire
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --array=0-79
#SBATCH --output=logs/ma_wza_hapfire_%A_%a.out
#SBATCH --error=logs/ma_wza_hapfire_%A_%a.out
# ---------------------------------------------------------------------------
# OLD (COARSE) BLOCK DEFINITION TEST (2026-07-27)
#
# Hypothesis: A. thaliana is a heavy selfer -> LD extends far -> the production
# clq0.9 blocks are too FINE. Measured on current non-SNP data:
#     hapFIRE : 703,789 records (100%) in 14,043 blocks, median 8 rec/block
#     clq0.9  : 417,692 records ( 59%) in 42,749 blocks, median 3 rec/block
# clq blocks are gap-y intervals, so clq0.9 DISCARDS 40.7% of records outright,
# while also giving less evidence per test and ~3x more tests.
#
# Runs, on the SAME current (post-Kf_w) per-record data:
#   isotonic  (current production SD correction)   -> _hf_isotonic
#   deg7      (what phase-1 actually used WITH these blocks; author's advice)
#             -> _hf_deg7    == the true "old configuration" on corrected data
# Cap 2000 (phase-1's cap; hapFIRE blocks are far larger than clq0.9's).
#
# Array index i (0..79): axis = AXES[i/4], cls = CLS4[i%4].
# ---------------------------------------------------------------------------
set -euo pipefail
# Repo root: BASH_SOURCE is useless under SLURM (script is spooled). Use
# SLURM_SUBMIT_DIR, walk up to the tree, and VERIFY -- never cd blindly.
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
while [ "$ROOT" != "/" ] && [ ! -d "$ROOT/analysis/grenenet_gea/phase1_replication" ]; do
  ROOT="$(dirname "$ROOT")"
done
if [ ! -d "$ROOT/analysis/grenenet_gea/phase1_replication" ]; then
  echo "FATAL: cannot locate repo root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-unset}" >&2; exit 1
fi
cd "$ROOT"
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_gea/phase1_replication
MA=$PR/results/multiaxis
IND=$MA/wza_in_hapfire

AXES=(bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 pc1)
CLS4=(snp sv smallindel nonsnp)
i=${SLURM_ARRAY_TASK_ID}
AXIS=${AXES[$(( i / 4 ))]}
C=${CLS4[$(( i % 4 ))]}
CAP=2000

echo "== host $(hostname) | hapfire $C $AXIS (cap=$CAP) | root=$ROOT =="

echo "-- reblock onto hapFIRE blocks --"
$PY -u $PR/multiaxis/reblock_blockdef.py --axis "$AXIS" --cls "$C" --blockdef hapfire

RAN=0
for M in kendall lfmm binomial; do
  IN=$IND/${M}_${C}_gen9_${AXIS}.csv
  [ -f "$IN" ] || { echo "FATAL: missing reblocked input $IN (cwd=$PWD)" >&2; exit 2; }
  RAN=$((RAN+1))

  echo "-- $M hf_isotonic --"
  $PY -u $PR/run_wza.py --model "$M" --class "$C" --gen 9 --climate "$AXIS" \
      --indir "$IND" --out "$MA/wza_hapfire" \
      --regime hf_isotonic --sd-fit isotonic --sample-snps $CAP

  echo "-- $M hf_deg7 (phase-1's own configuration) --"
  $PY -u $PR/run_wza.py --model "$M" --class "$C" --gen 9 --climate "$AXIS" \
      --indir "$IND" --out "$MA/wza_hapfire" \
      --regime hf_deg7 --sd-fit poly --poly-deg 7 --min-entries 10 --sample-snps $CAP
done
[ "$RAN" -eq 3 ] || { echo "FATAL: only $RAN/3 models ran for $C $AXIS" >&2; exit 3; }
echo "== ALLDONE_HAPFIRE $C $AXIS (models run: $RAN) =="

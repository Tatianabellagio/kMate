#!/bin/bash
#SBATCH --job-name=ma_wza_degsweep
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --array=0-79
#SBATCH --output=logs/ma_wza_degsweep_%A_%a.out
#SBATCH --error=logs/ma_wza_degsweep_%A_%a.out
# ---------------------------------------------------------------------------
# SD-correction DEGREE SWEEP (2026-07-27), to adjudicate the worry that the
# isotonic fix is over-conservative and destroyed real climate signal.
#
# Re-runs WZA on the CURRENT (post-Kf_w) per-record data under two additional
# regimes alongside the production isotonic one, so DEGREE is the only thing
# varying and can be compared against the already-present *_isotonic.csv:
#   deg2cur : canonical Booker deg-2, class-specific cap  (--poly-deg 2)
#   deg7cur : author-advised deg-7,   class-specific cap  (--poly-deg 7 --min-entries 10)
# Caps are the SAME class-specific values production uses (snp 1000, else 350),
# so cap is held constant and only the polynomial degree changes.
#
# Rationale for deg-7 (wza_investigation/RESULTS.md): deg-7 fits the empirical
# SNP-count->SD curve far better (RMSE 0.94 vs deg-2 2.06) and specifically
# lifts the 11-50-SNP window range -- exactly where the disputed climate genes
# (HSFA2 N=24, COR15A N=9, HSBP N=16) live. Its known danger is TAIL
# extrapolation, which the cap removes.
#
# Array index i (0..79): axis = AXES[i/4], cls = CLS4[i%4].  4 classes x 20 axes.
# Reblocking is NOT redone -- wza_in/ already holds the current reblocked input.
# ---------------------------------------------------------------------------
set -euo pipefail
# Repo root: SLURM copies the script to a spool dir, so BASH_SOURCE is USELESS here
# (it silently resolved to "/" on the first attempt and every task SKIPped to a
# 1-second "COMPLETED"). Use SLURM_SUBMIT_DIR (submit from the repo root) and walk
# up to the dir that actually contains the tree, then VERIFY -- never cd blindly.
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
while [ "$ROOT" != "/" ] && [ ! -d "$ROOT/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication" ]; do
  ROOT="$(dirname "$ROOT")"
done
if [ ! -d "$ROOT/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication" ]; then
  echo "FATAL: cannot locate repo root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-unset}" >&2
  exit 1
fi
cd "$ROOT"
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication
MA=$PR/results/multiaxis

AXES=(bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 pc1)
CLS4=(snp sv smallindel nonsnp)
i=${SLURM_ARRAY_TASK_ID}
AXIS=${AXES[$(( i / 4 ))]}
C=${CLS4[$(( i % 4 ))]}
CAP=350; [ "$C" = "snp" ] && CAP=1000

echo "== host $(hostname) | degsweep $C $AXIS (cap=$CAP) | root=$ROOT =="

RAN=0
for M in kendall lfmm binomial; do
  IN=$MA/wza_in/${M}_${C}_gen9_${AXIS}.csv
  # HARD failure, not a silent skip -- a missing input means the task did no work,
  # which previously looked identical to success in sacct.
  [ -f "$IN" ] || { echo "FATAL: missing input $IN (cwd=$PWD)" >&2; exit 2; }
  RAN=$((RAN+1))

  echo "-- $M deg2cur --"
  $PY -u $PR/run_wza.py --model "$M" --class "$C" --gen 9 --climate "$AXIS" \
      --indir "$MA/wza_in" --out "$MA/wza" \
      --regime deg2cur --sd-fit poly --poly-deg 2 --sample-snps $CAP

  echo "-- $M deg7cur --"
  $PY -u $PR/run_wza.py --model "$M" --class "$C" --gen 9 --climate "$AXIS" \
      --indir "$MA/wza_in" --out "$MA/wza" \
      --regime deg7cur --sd-fit poly --poly-deg 7 --min-entries 10 --sample-snps $CAP
done
[ "$RAN" -eq 3 ] || { echo "FATAL: only $RAN/3 models ran for $C $AXIS" >&2; exit 3; }
echo "== ALLDONE_DEGSWEEP $C $AXIS (models run: $RAN) =="

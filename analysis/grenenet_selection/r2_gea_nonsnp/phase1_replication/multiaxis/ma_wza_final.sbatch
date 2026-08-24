#!/bin/bash
#SBATCH --job-name=ma_wza_final
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --array=0-79
#SBATCH --output=logs/ma_wza_final_%A_%a.out
#SBATCH --error=logs/ma_wza_final_%A_%a.out
# ---------------------------------------------------------------------------
# REGENERATE THE PRODUCTION MULTIAXIS WZA SET ON THE SETTLED REGIME (2026-07-28)
#
#   blocks   clq0.9 BigLD islands, TILING assignment   (0.0% of records dropped;
#            strict interval containment dropped 39.7% of snp / 49.3% of SV)
#   cap      NONE  (--sample-snps 0; upstream's own default, never passed in
#            TBooker/WZA's own invocations)
#   SD fit   isotonic (monotone non-decreasing, flat past support; cannot go <=0,
#            which deg-2 does on 14-16 blocks -> NaN / fabricated p==0)
#   mean fit poly_clamp degree 5 (upstream fits the mean with a polynomial too;
#            clamping at the rolling-support edge is what stops the -289 blow-up)
# Full derivation: notebooks/cap_poly_decision.ipynb,
# wza_investigation/fit_inspection_{lfmm,kendall,binomial}.ipynb,
# wza_investigation/wza_fit_decision.py, and STATUS_clq90 "SETTLED 2026-07-28".
#
# Output: results/multiaxis/wza/wza_{model}_{cls}_gen9_{axis}_final.csv
#         20 axes x 4 classes x 3 models = 240 files.
#
# ⚠ SUPERSEDES the _isotonic / _deg2 / _deg2cur / _deg7cur sets already in that
#   directory. Those used strict (lossy) assignment and/or a cap and/or the
#   retired interp mean -- DO NOT CITE THEM. `_final` is the only current set.
#
# ⚠ Block *definition* (r2=0.9 vs 0.5) is still open; ma_wza_blockarms.sbatch
#   answers that. This regenerates the CURRENT production definition so the
#   downstream comparison/candidate-gene scripts have a coherent set to read.
#
# ⚠ MUST NOT RUN CONCURRENTLY WITH ma_wza_blockarms.sbatch. Both invoke
#   reblock_blockdef.py for clq0.9/tiling and therefore write the SAME
#   wza_in_clq09_tile/{model}_{cls}_gen9_{axis}.csv files, with identical array-index
#   -> (axis, class) mappings. df.to_csv is not atomic, so a concurrent run can read a
#   truncated per-record file and emit silently wrong WZA output. Submit with
#     sbatch --dependency=afterany:<blockarms_jobid> ma_wza_final.sbatch
#
# NOTE: paths are relative to the repo root -- the repo moved once and the old
# ma_wza.sbatch still hardcodes /global/scratch/users/tbellg/kmate.
# ---------------------------------------------------------------------------
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
while [ "$ROOT" != "/" ] && [ ! -d "$ROOT/analysis/grenenet_selection/phase1_replication" ]; do
  ROOT="$(dirname "$ROOT")"
done
if [ ! -d "$ROOT/analysis/grenenet_selection/phase1_replication" ]; then
  echo "FATAL: cannot locate repo root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-unset}" >&2; exit 1
fi
cd "$ROOT"
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_selection/phase1_replication
MA=$PR/results/multiaxis
IND=$MA/wza_in_clq09_tile

AXES=(bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 pc1)
CLS4=(snp sv smallindel nonsnp)
i=${SLURM_ARRAY_TASK_ID}
AXIS=${AXES[$(( i / 4 ))]}
C=${CLS4[$(( i % 4 ))]}

echo "== host $(hostname) | ma_wza_final $C $AXIS | root=$ROOT =="
echo "== regime: clq0.9 TILING | no cap | isotonic SD | poly_clamp deg5 mean =="

# tiling reblock (idempotent -- rewrites the same per-record file the arms run uses)
$PY -u $PR/multiaxis/reblock_blockdef.py --axis "$AXIS" --cls "$C" --blockdef clq0.9 --how tiling

for M in kendall lfmm binomial; do
  IN=$IND/${M}_${C}_gen9_${AXIS}.csv
  [ -f "$IN" ] || { echo "FATAL: missing $IN" >&2; exit 2; }
  $PY -u $PR/run_wza.py --model "$M" --class "$C" --gen 9 --climate "$AXIS" \
      --indir "$IND" --out "$MA/wza" --regime final \
      --sd-fit isotonic --mean-fit poly_clamp --mean-poly-deg 5 --sample-snps 0
done

echo "== ALLDONE_MA_WZA_FINAL $C $AXIS =="

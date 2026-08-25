#!/bin/bash
#SBATCH --job-name=ma_wza_blockarms
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --array=0-79
#SBATCH --output=logs/ma_wza_blockarms_%A_%a.out
#SBATCH --error=logs/ma_wza_blockarms_%A_%a.out
# ---------------------------------------------------------------------------
# FOUR-ARM BLOCK-PARTITION TEST (2026-07-27)
#
# Separates two confounded problems in the production GEA:
#   (a) strict interval containment DISCARDED 40-49% of all records
#       (snp 39.7%, SV 49.3%, nonsnp 40.7%) -- a methodological error; HapFM's
#       convert_fine_genomewide_breakpoints exists precisely to avoid it.
#   (b) clq0.9 blocks may simply be too FINE for a heavy selfer (long LD).
#
# Arms (all on the SAME current post-Kf_w per-record data, WZA isotonic):
#   clq0.9 strict  = production baseline           -> wza_arms/..._a09s
#   clq0.9 tiling  = same LD boundaries, 0% lost   -> wza_arms/..._a09t   (isolates (a))
#   clq0.5 strict  = coarser, still lossy          -> wza_arms/..._a05s
#   clq0.5 tiling  = coarser AND 0% lost           -> wza_arms/..._a05t   (isolates (b))
#
# Array index i (0..79): axis = AXES[i/4], cls = CLS4[i%4]; all 4 arms per task.
# ---------------------------------------------------------------------------
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
while [ "$ROOT" != "/" ] && [ ! -d "$ROOT/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication" ]; do
  ROOT="$(dirname "$ROOT")"
done
if [ ! -d "$ROOT/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication" ]; then
  echo "FATAL: cannot locate repo root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-unset}" >&2; exit 1
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

echo "== host $(hostname) | blockarms $C $AXIS | root=$ROOT =="

# ---------------------------------------------------------------------------
# REGIME UPDATE (2026-07-28) -- CAP REMOVED, so the arms are now a CLEAN 2x2.
#
# The previous version scaled the cap with the arm (strict 350/1000, tiling
# 700/2000), which meant a09s -> a09t changed BOTH the assignment rule AND the cap.
# The tiling effect could therefore not be attributed. Settled in
# notebooks/cap_poly_decision.ipynb (see STATUS_clq90 "SETTLED 2026-07-28"):
#   * no cap        -- isotonic cannot give a negative SD, only 8-15 blocks/class sit
#                      past rolling support, and the flat clip is conservative there;
#                      capping instead swaps the block's Z for an unseeded
#                      resample-mean (wza_script.py:204-207) => irreproducible.
#   * isotonic SD   -- deg7 out-of-sample RMSE blows up (17.7 vs 0.13 on sv/clq0.5)
#   * poly_clamp mean (deg 5) -- upstream fits the mean with an unclamped deg-2; const was
#                      WORST of four at the large blocks (|err| 2.81 vs 1.60)
# Every arm now uses the identical correction, so the ONLY thing varying across arms
# is blockdef x assignment. Outputs overwrite the confounded pre-2026-07-28 ones by
# design -- those must not be cited (see the stale-results banners in STATUS_clq90).
# ---------------------------------------------------------------------------
run_arm () {
  local BD="$1" HOW="$2" SUF="$3"
  local TAG="${BD/./}"; [ "$HOW" = "tiling" ] && TAG="${TAG}_tile"
  local IND=$MA/wza_in_${TAG}
  echo "---- arm $SUF ($BD/$HOW, no cap, isotonic SD + deg5-clamped mean) ----"
  $PY -u $PR/multiaxis/reblock_blockdef.py --axis "$AXIS" --cls "$C" --blockdef "$BD" --how "$HOW"
  for M in kendall lfmm binomial; do
    local IN=$IND/${M}_${C}_gen9_${AXIS}.csv
    [ -f "$IN" ] || { echo "FATAL: missing $IN" >&2; exit 2; }
    $PY -u $PR/run_wza.py --model "$M" --class "$C" --gen 9 --climate "$AXIS" \
        --indir "$IND" --out "$MA/wza_arms" \
        --regime "$SUF" --sd-fit isotonic --mean-fit poly_clamp --mean-poly-deg 5 --sample-snps 0
  done
}

run_arm clq0.9 strict a09s
run_arm clq0.9 tiling a09t
run_arm clq0.5 strict a05s
run_arm clq0.5 tiling a05t

echo "== ALLDONE_BLOCKARMS $C $AXIS =="

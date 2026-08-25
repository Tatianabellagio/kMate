#!/bin/bash
#SBATCH --job-name=wza_meanfit_strat
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/permnull_regime_%j.out
#SBATCH --error=logs/permnull_regime_%j.out
# ---------------------------------------------------------------------------
# Permutation-null CALIBRATION of the settled WZA regime (2026-07-28):
#   tiling blocks + isotonic SD + isotonic_auto mean + NO CAP
# against upstream deg-2/no-cap, the old kMate capped+interp production, and the
# same proposed fit WITH a cap (to isolate the cap's contribution).
#
# bio1 is site-permuted across the 31 sites -> signal-free, so any Bonferroni/FDR hit
# is a correction artifact. This is the check the fit diagnostics in
# notebooks/cap_poly_decision.ipynb could NOT do: whether isotonic's flat SD
# extrapolation past the support edge is anti-conservative for the largest blocks.
#
# NB the per-record test uses n=352 pools but there are only 31 independent climates,
# so some null heavy-tailedness is expected in EVERY regime. The comparison across
# regimes is what matters, not the absolute count.
#
# NOTE: paths are relative to the repo root (repo has moved once; do not hardcode).
# ---------------------------------------------------------------------------
set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
while [ "$ROOT" != "/" ] && [ ! -d "$ROOT/analysis/grenenet_selection/wza/investigation" ]; do
  ROOT="$(dirname "$ROOT")"
done
if [ ! -d "$ROOT/analysis/grenenet_selection/wza/investigation" ]; then
  echo "FATAL: cannot locate repo root from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-unset}" >&2
  exit 1
fi
cd "$ROOT/analysis/grenenet_selection/wza/investigation"
mkdir -p logs
echo "== host $(hostname) | root=$ROOT =="
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python -u mean_fit_stratified.py \
    
echo "== ALLDONE_MEANFIT_STRAT =="

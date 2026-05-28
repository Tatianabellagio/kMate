#!/bin/bash
# =============================================================================
# Re-run p80 truth + cactus_em after the MAR-semantics patch (2026-05-21).
#
# Truth: compute_recomb_truth.py now does (Σw·cn_var) / (Σw·cn_var_called)
#        and emits an `info` column.
# cactus_em window mode now uses cn_var_called for the called-mask projection
# (previously it silently treated `.` as REF). Output TSV gains info, n_called,
# se columns.
#
# Deletes old truth files + result TSVs so the dependency chain re-runs cleanly.
# =============================================================================
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/benchmarks/p80

PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
TRUTH_SCRIPT=/global/scratch/users/tbellg/kmate/sims/visor_freqk/scripts/compute_recomb_truth.py

CN_VAR=$PWD/data/cn_var_p80.cn_var.npz
CN_VAR_META=$PWD/data/cn_var_p80.meta.npz
for f in "$CN_VAR" "$CN_VAR_META" "${CN_VAR%.cn_var.npz}.cn_var_called.npz"; do
    [ -s "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

# Remove old truth + result TSVs (force regeneration with the new semantics)
echo "[$(date)] removing old truth + result outputs..."
for REG in n50_g1 n50_g3; do
    rm -fv sims/cov10_${REG}_s42_hotspots_p80_chr1/recomb_truth.tsv.gz
done
for METHOD in global star2; do
    for REG in n50_g1 n50_g3; do
        rm -fv results/cactus_em_${METHOD}/${REG}/p80_${REG}_cov10_s42.tsv
    done
done

# Recompute truth (fast — ~1-2 min/regime)
echo
echo "[$(date)] === Recomputing truth (MAR convention) ==="
for REG in n50_g1 n50_g3; do
    WORK=$PWD/sims/cov10_${REG}_s42_hotspots_p80_chr1
    echo "[$(date)] truth $REG"
    $PY $TRUTH_SCRIPT \
        --ancestry $WORK/ancestry.tsv \
        --weights  $WORK/pool_weights.tsv \
        --cn-var   $CN_VAR \
        --cn-var-meta $CN_VAR_META \
        --out      $WORK/recomb_truth.tsv.gz
done

# Resubmit cactus_em — depend on this script finishing
echo
echo "[$(date)] === Submitting cactus_em runs ==="
C1G=$(sbatch --parsable scripts/07_run_cactus_em_p80.sh n50_g1 global)
C1S=$(sbatch --parsable scripts/07_run_cactus_em_p80.sh n50_g1 star2)
C3G=$(sbatch --parsable scripts/07_run_cactus_em_p80.sh n50_g3 global)
C3S=$(sbatch --parsable scripts/07_run_cactus_em_p80.sh n50_g3 star2)
cat <<EOF
Submitted cactus_em jobs:
  $C1G   global  n50_g1
  $C1S   star2   n50_g1
  $C3G   global  n50_g3
  $C3S   star2   n50_g3
Monitor: squeue -u \$USER | grep p80
EOF

#!/bin/bash
# Rebuild everything that carries a genomic-context label, after the project-wide rule change:
# promoter window 2 kb -> 1 kb (_sv_sfs_context.PROM_BP) with TE still above promoter.
# CDS / UTR / intron / TE assignments are unaffected -- only the promoter/intergenic boundary moves --
# but the burden tables and the event catalogue both store the label, so both must be rebuilt.
#SBATCH --job-name=regen_ctx
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=0-6
#SBATCH --output=logs/regen_ctx_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
case $SLURM_ARRAY_TASK_ID in
  0) echo "[catalog]";              $PY _sv_event_catalog.py ;;
  1) echo "[burden 0.9]";           SV_UNIT=event $PY _sv_burden_build.py ;;
  2) echo "[burden 0.7]";           SV_UNIT=event SV_CALL_MIN=0.7 $PY _sv_burden_build.py ;;
  3) echo "[burden 0.5]";           SV_UNIT=event SV_CALL_MIN=0.5 $PY _sv_burden_build.py ;;
  4) echo "[strata freq]";          STRATUM=freq $PY _sv_burden_strata.py ;;
  5) echo "[strata size]";          STRATUM=size $PY _sv_burden_strata.py ;;
  6) echo "[strata callrate]";      STRATUM=callrate SV_CALL_MIN=0.5 MIN_CAR=2 $PY _sv_burden_strata.py ;;
esac

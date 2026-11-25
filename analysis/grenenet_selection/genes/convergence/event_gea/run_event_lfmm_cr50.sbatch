#!/bin/bash
#SBATCH --job-name=event_cr50
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=1:30:00
#SBATCH --requeue
#SBATCH --array=0-43
#SBATCH --output=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/genes/convergence/event_gea/cr50/logs/event_gea_%A_%a.out
#SBATCH --error=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/genes/convergence/event_gea/cr50/logs/event_gea_%A_%a.out
# GEA event arm (step 4 of _submit_gea.sh): LFMM ridge K = 16, raw p, for sv and nonsnp Truvari-event units x 22
# climate axes (44 tasks), with the record arm's input writer and LFMM script.
# In:  event_gea/cr50_arch3v2/{cls}_gen9_{dp.npy,units.csv} (step 3); the pools' climate from class_matrices_cr50_arch3v2
# Out: event_gea/cr50_arch3v2/lfmm_event_{cls}_gen9_{axis}.csv.gz  (units.csv + pval, B, z)
set -euo pipefail
G=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection
V=${KMATE_PANEL:-arch3v2}; S=""; [ "$V" = arch3 ] || S="_$V"     # panel version (lib.PANEL): outputs under *_arch3v2 dirs
D=$G/genes/convergence/event_gea/cr50$S
PR=$G/r2_gea_nonsnp/phase1_replication
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
AXES=(bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 pc1 pc2 pc3)
CLS=(nonsnp sv)
AXIS=${AXES[$(( SLURM_ARRAY_TASK_ID / 2 ))]}; C=${CLS[$(( SLURM_ARRAY_TASK_ID % 2 ))]}
WORK=/global/scratch/projects/fc_moilab/tbellg/kmate/tmp_event_gea_cr50$S/work_${C}_${AXIS}
echo "== host $(hostname) | event arm $C $AXIS =="

$PY -u $PR/build_lfmm_input.py --class "$C" --climate "$AXIS" --cmdir "$PR/results/class_matrices_cr50$S" \
    --out "$WORK" --dp "$D/${C}_gen9_dp.npy"
$R $PR/multiaxis/run_lfmm_nogif.R "$WORK/lfmm_${C}_gen9" 16 "$WORK/raw.csv" -

# p, B, z back onto the units by POSITION (row i of units.csv is column i of Y)
$PY - "$C" "$AXIS" "$WORK" "$D" <<'PYEOF'
import sys, pandas as pd
cls, ax, work, D = sys.argv[1:]
T = pd.read_csv(f"{D}/{cls}_gen9_units.csv"); r = pd.read_csv(f"{work}/raw.csv")
assert len(T) == len(r)
T[["pval", "B", "z"]] = r[["pval", "B", "z"]].to_numpy()
T.to_csv(f"{D}/lfmm_event_{cls}_gen9_{ax}.csv.gz", index=False)
thr = 0.05 / len(T)
print(f"  {cls} {ax}: {len(T):,} units, Bonferroni {thr:.2e}: {int((T.pval < thr).sum())} hits "
      f"({int((T.is_sv & (T.pval < thr)).sum())} SV events)", flush=True)
PYEOF
rm -rf "$WORK"
echo "== ALLDONE_EVENT_GEA $C $AXIS =="

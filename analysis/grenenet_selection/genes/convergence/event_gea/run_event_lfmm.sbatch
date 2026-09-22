#!/bin/bash
#SBATCH --job-name=event_gea
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=1:30:00
#SBATCH --requeue
#SBATCH --array=0-43
#SBATCH --output=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/genes/convergence/event_gea/logs/event_gea_%A_%a.out
#SBATCH --error=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/genes/convergence/event_gea/logs/event_gea_%A_%a.out
# Event-level non-SNP GEA: the GEA's LFMM (ridge, K = 16, raw p) on truvari-event units.
# Array index i: axis = AXES[i/2], class = CLS[i%2]. Inputs from build_event_matrix.py.
set -euo pipefail
D=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/genes/convergence/event_gea
CM=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/class_matrices
AXES=(bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 pc1 pc2 pc3)
CLS=(nonsnp sv)
AXIS=${AXES[$(( SLURM_ARRAY_TASK_ID / 2 ))]}; C=${CLS[$(( SLURM_ARRAY_TASK_ID % 2 ))]}
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
WORK=/global/scratch/projects/fc_moilab/tbellg/kmate/tmp_event_gea/work_${C}_${AXIS}; mkdir -p "$WORK"
echo "== host $(hostname) | event GEA $C $AXIS =="
$PY - "$C" "$AXIS" "$WORK" "$D" "$CM" <<'PYEOF'
import sys, numpy as np, pandas as pd
cls, ax, work, D, CM = sys.argv[1:]
dp = np.load(f"{D}/{cls}_gen9_dp.npy")
pools = pd.read_csv(f"{CM}/gen9.pools.csv")
assert dp.shape[0] == len(pools)
stem = f"{work}/ev_{cls}_gen9"
dp.astype(np.float64).tofile(f"{stem}_Y.f64")
open(f"{stem}_dims.txt", "w").write(f"{dp.shape[0]} {dp.shape[1]}\n")
x = pools[ax].to_numpy(float); z = (x - x.mean()) / x.std(ddof=0)
pd.DataFrame({ax: z}).to_csv(f"{stem}_env.csv", index=False)
print(f"  Y [{dp.shape[0]} x {dp.shape[1]:,}], env z-scored {ax}", flush=True)
PYEOF
$R $D/../lfmm_signed/run_lfmm_signed.R "$WORK/ev_${C}_gen9" 16 "$WORK/raw.csv"
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

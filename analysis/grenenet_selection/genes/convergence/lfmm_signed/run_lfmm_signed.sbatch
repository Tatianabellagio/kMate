#!/bin/bash
#SBATCH --job-name=lfmm_sign
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=6:00:00
#SBATCH --requeue
#SBATCH --array=0-21
#SBATCH --output=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/genes/convergence/lfmm_signed/logs/lfmm_sign_%A_%a.out
#SBATCH --error=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/genes/convergence/lfmm_signed/logs/lfmm_sign_%A_%a.out
# Pooled non-SNP LFMM, gen9, K = 16 -- the GEA's own model -- re-run per axis to recover
# the sign of the association (B, z). Inputs rebuilt exactly as ma_nonsnp_lfmm.sbatch does.
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate
AXES=(bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 pc1 pc2 pc3)
AXIS=${AXES[$SLURM_ARRAY_TASK_ID]}; C=nonsnp
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
GS=analysis/grenenet_selection
PR=$GS/r2_gea_nonsnp/phase1_replication
D=$GS/genes/convergence/lfmm_signed
WORK=/global/scratch/projects/fc_moilab/tbellg/kmate/tmp_lfmm_signed/work_${C}_${AXIS}; mkdir -p "$WORK"
echo "== host $(hostname) | signed lfmm nonsnp gen9 $AXIS =="
$PY -u $PR/build_lfmm_input.py --class "$C" --gen 9 --climate "$AXIS" --out "$WORK"
$R $D/run_lfmm_signed.R "$WORK/lfmm_${C}_gen9" 16 "$WORK/raw_${AXIS}.csv"
$PY - "$AXIS" <<'PYEOF'
import sys, numpy as np, pandas as pd
ax = sys.argv[1]
GS = "analysis/grenenet_selection"
rec = pd.read_csv(f"{GS}/r2_gea_nonsnp/phase1_replication/results/class_matrices/nonsnp_gen9.records.csv")
new = pd.read_csv(f"/global/scratch/projects/fc_moilab/tbellg/kmate/tmp_lfmm_signed/work_nonsnp_{ax}/raw_{ax}.csv")
assert len(new) == len(rec), (len(new), len(rec))
old = pd.read_csv(f"{GS}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile/lfmm_nonsnp_gen9_{ax}.csv")
out = rec.rename(columns={"maf": "MAF"})[["chrom", "pos", "ref_len", "alt_len", "MAF"]].copy()
out[["pval", "B", "z"]] = new[["pval", "B", "z"]].to_numpy()
# the saved GEA file is in the same record order; check the p-values reproduce
assert len(old) == len(out) and (old.pos.to_numpy() == out.pos.to_numpy()).all(), "record order differs"
d = np.abs(np.log10(np.clip(old.pval.to_numpy(), 1e-300, 1)) - np.log10(np.clip(out.pval.to_numpy(), 1e-300, 1)))
print(f"  reproduction vs saved GEA p: max |dlog10 p| = {np.nanmax(d):.2e}, median {np.nanmedian(d):.2e}")
if np.nanmax(d) > 1e-6:
    raise SystemExit("p-values do NOT reproduce the saved GEA -- not writing a signed table")
out["block"] = old.block.to_numpy()
out.to_csv(f"{GS}/genes/convergence/lfmm_signed/lfmm_nonsnp_gen9_{ax}_signed.csv.gz", index=False)
print(f"  wrote signed table; share of records with B > 0: {(out.B > 0).mean():.3f}")
PYEOF
rm -rf "$WORK"
echo "== ALLDONE_LFMM_SIGNED $AXIS =="

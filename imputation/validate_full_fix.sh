#!/bin/bash
#SBATCH --job-name=val_v2
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/val_v2_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/val_v2_%j.err

# =============================================================================
# validate_full_fix.sh
# After job 57794 (full_fix) finishes, validate the corrected 231-founder
# cn matrices against the SEEDMIX recipe truth across all 8 replicates.
# Should give slope ≈ 1.0 with R² ≈ 0.99 (confirming Chr4 result on full genome).
# =============================================================================
set -euo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq

PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python

# Run per_sample_driver in global mode on SEEDMIX_S1 with corrected cn
mkdir -p results/seedmix_231_v2

$PYTHON -u src/per_sample_driver.py \
    --cn-kmer-prefix data/cn_full_231_v2/cn \
    --cn-var data/cn_var_231_v2.cn_var.npz \
    --cn-var-meta data/cn_var_231_v2.meta.npz \
    --reads /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz \
            /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz \
    --sample SEEDMIX_S1_v2 \
    --out results/seedmix_231_v2/SEEDMIX_S1.tsv \
    --threads 8 \
    --block-mode global

echo ""
echo "[$(date)] === slope diagnostic on SEEDMIX_S1, full-genome corrected cn ==="
$PYTHON -u <<'PYEOF'
import numpy as np, pandas as pd
from scipy.sparse import load_npz

cn_var = load_npz("data/cn_var_231_v2.cn_var.npz")
meta = np.load("data/cn_var_231_v2.meta.npz", allow_pickle=True)
founders = list(meta["founders"])
recipe = pd.read_csv("../data/seedmix_recipe_normalized.tsv", sep="\t")
rd = dict(zip(recipe.ID.astype(str), recipe.seed_prop))
h_truth = np.array([rd.get(str(f), 0.0) for f in founders]); h_truth /= h_truth.sum()
af_truth = (h_truth.astype(np.float32) @ cn_var.toarray()).astype(np.float64)

df = pd.read_csv("results/seedmix_231_v2/SEEDMIX_S1.tsv", sep="\t")
af_pred = df["alt_freq"].to_numpy()

slope, intercept = np.polyfit(af_truth, af_pred, 1)
r = np.corrcoef(af_pred, af_truth)[0,1]
r2 = 1 - ((af_pred-af_truth)**2).sum() / ((af_truth-af_truth.mean())**2).sum()
rmse = float(np.sqrt(np.mean((af_pred-af_truth)**2)))

print(f"  slope     = {slope:.4f}  (was 1.43 with bug; Chr4-only test gave 1.003)")
print(f"  intercept = {intercept:.5f}")
print(f"  R² (raw)  = {r2:.4f}")
print(f"  Pearson r = {r:.4f}")
print(f"  RMSE      = {rmse:.4f}")
print(f"\n  → if slope ≈ 1.0 and R² ≈ 0.99, the full-genome fix is confirmed.")
PYEOF
echo "[$(date)] DONE"

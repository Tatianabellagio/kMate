#!/bin/bash
#SBATCH --job-name=chr1_hdiag
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1/E2_hdiag_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1/E2_hdiag_%j.err
set -euo pipefail

# Why does arch3 AF (h_v3 @ cn_var) diverge from recipe_fixed (uniform_h @ cn_var)?
# Hypothesis: the h vector is non-uniform; at outlier records, carrier set correlates
# with the deviation of h from uniform.
# This is signal (real seed-mix composition), not a bug.

cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python

$PY -u <<'PYEOF'
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
import json

print('=== Load h vector + cn_var + meta ===')
h_data = np.load('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.h_per_chrom.npz', allow_pickle=True)
h = h_data['Chr1'].astype(np.float64)
h_founders = list(h_data['founders'])
F = len(h)
print(f'h shape: {h.shape}, sum: {h.sum():.6f}, mean: {h.mean():.6f}, expected uniform: {1.0/F:.6f}')

# Founder split: cactus vs PG
with open('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/founder_split_cactus_pg.json') as f:
    split = json.load(f)
cactus_set = set(map(str, split['cactus']))
pg_set = set(map(str, split['PG']))
is_cactus = np.array([str(s) in cactus_set for s in h_founders])
is_pg = np.array([str(s) in pg_set for s in h_founders])

print(f'\n=== h vector inspection ===')
print(f'min: {h.min():.6e}  max: {h.max():.4e}  std: {h.std():.4e}')
print(f'CV (std/mean): {h.std()/h.mean():.3f}  (0 = perfectly uniform)')
print(f'fraction of founders within ±10% of uniform: {np.mean(np.abs(h - 1/F) < 0.1/F):.3f}')
print(f'top 10 founders by h:')
order = np.argsort(h)[::-1]
for i in order[:10]:
    cls = 'cact' if is_cactus[i] else ('PG' if is_pg[i] else '???')
    print(f'  {h_founders[i]:>8}  h={h[i]:.6f}  ({h[i]/(1/F):.3f}× uniform)  class={cls}')
print(f'bottom 10 founders by h:')
for i in order[-10:]:
    cls = 'cact' if is_cactus[i] else ('PG' if is_pg[i] else '???')
    print(f'  {h_founders[i]:>8}  h={h[i]:.6e}  ({h[i]/(1/F):.3f}× uniform)  class={cls}')

print(f'\n=== h by founder class ===')
print(f'cactus mean h: {h[is_cactus].mean():.6f}  ({h[is_cactus].mean()/(1/F):.3f}× uniform)')
print(f'PG     mean h: {h[is_pg].mean():.6f}  ({h[is_pg].mean()/(1/F):.3f}× uniform)')
print(f'cactus sum h:  {h[is_cactus].sum():.4f}  (expected uniform: {is_cactus.sum()/F:.4f})')
print(f'PG     sum h:  {h[is_pg].sum():.4f}  (expected uniform: {is_pg.sum()/F:.4f})')
print(f'h-bias (cact per-founder / PG per-founder): {(h[is_cactus].mean())/(h[is_pg].mean()):.3f}')

# === Now find outlier records (arch3 vs recipe deviation) and look at their carrier composition ===
print('\n=== Load cn_var (raw arch3) + project both arch3 + recipe ===')
cn_var = load_npz('cn_var_231_arch3_chr1.cn_var.npz').tocsr()
cn_var_called = load_npz('cn_var_231_arch3_chr1.cn_var_called.npz').tocsr()
meta = np.load('cn_var_231_arch3_chr1.meta.npz', allow_pickle=True)
N = cn_var.shape[1]

h_uniform = np.full(F, 1.0/F)
recipe_n = h_uniform @ cn_var
recipe_d = h_uniform @ cn_var_called
recipe_af = np.where(recipe_d > 0, recipe_n / recipe_d, np.nan)

arch3_n = h @ cn_var
arch3_d = h @ cn_var_called
arch3_af = np.where(arch3_d > 0, arch3_n / arch3_d, np.nan)

# Outlier records (|arch3 - recipe| > 0.10)
diff = arch3_af - recipe_af
outl_mask = np.abs(diff) > 0.10
n_outl = outl_mask.sum()
print(f'\nrecords with |arch3 - recipe| > 0.10: {n_outl:,} ({100*n_outl/N:.3f}%)')

print(f'\n=== Direction split among outliers ===')
above = (diff > 0.10).sum(); below = (diff < -0.10).sum()
print(f'  arch3 > recipe (high-h founders carry more): {above:,} ({100*above/n_outl:.1f}%)')
print(f'  arch3 < recipe (low-h founders carry more):  {below:,} ({100*below/n_outl:.1f}%)')

# At each outlier record, compute the "h-bias index" of its carrier set:
#   bias(r) = (mean h of carriers) / mean h overall - 1
#  Positive: carriers tend to be high-h founders → arch3 > recipe
#  Negative: carriers tend to be low-h founders  → arch3 < recipe
print('\n=== h-bias of carrier set at outlier records ===')
# Sample 1000 outliers for speed
outl_idx = np.where(outl_mask)[0]
sample = np.random.RandomState(42).choice(outl_idx, size=min(5000, len(outl_idx)), replace=False)
h_bias_carriers = []
fmiss_outl = []
ac_outl = []
for i in sample:
    col = cn_var[:, i].toarray().ravel()
    col_k = cn_var_called[:, i].toarray().ravel()
    if col.sum() == 0:
        h_bias_carriers.append(np.nan)
    else:
        mean_h_car = (col * h).sum() / col.sum()
        h_bias_carriers.append(mean_h_car / h.mean() - 1)
    fmiss_outl.append(1 - col_k.sum()/F)
    ac_outl.append(col.sum())
h_bias_carriers = np.array(h_bias_carriers)
print(f'sampled outliers: {len(sample):,}')
print(f'h-bias of carriers (mean_h_carriers / mean_h - 1):')
print(f'  mean: {np.nanmean(h_bias_carriers):+.4f}')
print(f'  median: {np.nanmedian(h_bias_carriers):+.4f}')
print(f'  P10/P90: {np.nanpercentile(h_bias_carriers, [10, 90])}')

# Correlation between carrier h-bias and AF deviation
diff_sample = diff[sample]
finite = np.isfinite(h_bias_carriers) & np.isfinite(diff_sample)
if finite.sum() > 100:
    r = np.corrcoef(h_bias_carriers[finite], diff_sample[finite])[0,1]
    print(f'\nPearson r (carrier h-bias vs arch3-recipe deviation): {r:.4f}')
    print(f'  → if r ≈ 1, deviation is fully explained by h-bias of carriers (signal, not bug)')
    print(f'  → if r ≈ 0, deviation is random (artifact)')

# Check whether outliers correlate with cactus/PG founder asymmetry
print('\n=== Cactus vs PG carrier composition at outliers ===')
n_c_only = 0; n_p_only = 0; n_both = 0
for i in sample[:1500]:
    col = cn_var[:, i].toarray().ravel() > 0
    has_c = col[is_cactus].any()
    has_p = col[is_pg].any()
    if has_c and not has_p: n_c_only += 1
    elif has_p and not has_c: n_p_only += 1
    else: n_both += 1
print(f'  carriers in BOTH classes:    {n_both:>4,}  ({100*n_both/1500:.1f}%)')
print(f'  carriers in CACTUS only:     {n_c_only:>4,}  ({100*n_c_only/1500:.1f}%)')
print(f'  carriers in PG only:         {n_p_only:>4,}  ({100*n_p_only/1500:.1f}%)')

# Spot-check 10 outliers: who carries the variant, what's their h?
print('\n=== Spot-check: 10 random outliers — carrier identity ===')
spot = np.random.RandomState(7).choice(outl_idx, size=10, replace=False)
for i in spot:
    col = cn_var[:, i].toarray().ravel() > 0
    car_idx = np.where(col)[0]
    car_h = h[car_idx]
    car_names = [h_founders[j] for j in car_idx]
    car_cls = ['cact' if is_cactus[j] else ('PG' if is_pg[j] else '?') for j in car_idx]
    print(f'\nrecord {i}: pos={meta["pos"][i]} {meta["ref"][i]}>{meta["alt"][i]}')
    print(f'  recipe_af={recipe_af[i]:.4f}  arch3_af={arch3_af[i]:.4f}  Δ={arch3_af[i]-recipe_af[i]:+.4f}')
    print(f'  n_carriers={len(car_idx)}, mean h: {car_h.mean():.4e}, ratio to uniform: {car_h.mean()/h.mean():.3f}')
    top3 = np.argsort(car_h)[::-1][:3]
    print(f'  top-h carriers: {[(car_names[j], f"{car_h[j]:.2e}", car_cls[j]) for j in top3]}')
PYEOF

echo
echo "[$(date)] DONE E2"

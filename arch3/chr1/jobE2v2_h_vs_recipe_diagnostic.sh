#!/bin/bash
#SBATCH --job-name=chr1_hdiag2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --output=logs/E2v2_hdiag_%j.out
#SBATCH --error=logs/E2v2_hdiag_%j.err
mkdir -p logs
set -euo pipefail

# Fast version of the h-vs-recipe diagnostic — fully vectorized.
# At each record r, compute mean_h(carriers at r) = (h @ cn_var)[r] / AC[r]
# carrier_hbias[r] = mean_h(carriers) / mean(h) - 1
# Then correlate carrier_hbias with (arch3_af - recipe_af) at outlier records.

cd /global/scratch/users/tbellg/kmate/arch3/chr1
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

$PY -u <<'PYEOF'
import numpy as np
from scipy.sparse import load_npz
import json

print('=== Load h + cn_var ===')
h_data = np.load('/global/scratch/users/tbellg/kmate/scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.h_per_chrom.npz', allow_pickle=True)
h = h_data['Chr1'].astype(np.float64)
h_founders = list(h_data['founders'])
F = len(h)

with open('/global/scratch/users/tbellg/kmate/data/founder_split_cactus_pg.json') as f:
    split = json.load(f)
is_cactus = np.array([str(s) in set(map(str, split['cactus'])) for s in h_founders])
is_pg = np.array([str(s) in set(map(str, split['PG'])) for s in h_founders])

cn_var = load_npz('cn_var_231_arch3_chr1.cn_var.npz').tocsr()
cn_var_called = load_npz('cn_var_231_arch3_chr1.cn_var_called.npz').tocsr()
meta = np.load('cn_var_231_arch3_chr1.meta.npz', allow_pickle=True)
N = cn_var.shape[1]
print(f'cn_var: {cn_var.shape}, {cn_var.nnz:,} nnz; {N:,} records')

# Project both
h_uni = np.full(F, 1.0/F)
arch3_n = h @ cn_var          # (N,) — h-weighted carrier sum per record
arch3_d = h @ cn_var_called   # (N,) — h-weighted called sum per record
recipe_n = h_uni @ cn_var
recipe_d = h_uni @ cn_var_called
arch3_af  = np.where(arch3_d > 0,  arch3_n  / arch3_d,  np.nan)
recipe_af = np.where(recipe_d > 0, recipe_n / recipe_d, np.nan)

# Vectorized per-record carrier stats
ac = np.asarray(cn_var.sum(axis=0)).ravel()    # (N,) carrier count per record
an = np.asarray(cn_var_called.sum(axis=0)).ravel()  # (N,) called count per record
mean_h_carriers = np.where(ac > 0, arch3_n / ac, np.nan)
mean_h_all = h.mean()
carrier_hbias = mean_h_carriers / mean_h_all - 1.0  # >0: carriers tend to be high-h founders

# Per-class carrier counts (vectorized)
cn_var_c = cn_var[is_cactus, :]
cn_var_p = cn_var[~is_cactus, :]
ac_c = np.asarray(cn_var_c.sum(axis=0)).ravel()
ac_p = np.asarray(cn_var_p.sum(axis=0)).ravel()

# Outlier mask
diff = arch3_af - recipe_af
out_mask = np.abs(diff) > 0.10
above_mask = diff > 0.10
below_mask = diff < -0.10
n_out = out_mask.sum()

print(f'\n=== Recipe-vs-arch3 outliers (|Δ|>0.10) ===')
print(f'records: {n_out:,} ({100*n_out/N:.3f}%)')
print(f'  arch3 > recipe (above):  {above_mask.sum():,} ({100*above_mask.sum()/n_out:.1f}%)')
print(f'  arch3 < recipe (below):  {below_mask.sum():,} ({100*below_mask.sum()/n_out:.1f}%)')

# === Pearson r: carrier_hbias vs diff at outliers ===
print(f'\n=== Pearson r between carrier h-bias and (arch3 − recipe) deviation ===')
for mask, label in [(out_mask, 'outliers only'),
                     (np.isfinite(diff) & np.isfinite(carrier_hbias), 'ALL records'),
                     (above_mask, 'above outliers'),
                     (below_mask, 'below outliers')]:
    m = mask & np.isfinite(diff) & np.isfinite(carrier_hbias)
    if m.sum() < 5: continue
    r = np.corrcoef(carrier_hbias[m], diff[m])[0,1]
    print(f'  {label:<24} (n={m.sum():>8,}):  r = {r:+.4f}')

# === carrier h-bias distribution at outliers vs non-outliers ===
print(f'\n=== Distribution of carrier h-bias ===')
nout_mask = ~out_mask & np.isfinite(carrier_hbias)
out_finite = out_mask & np.isfinite(carrier_hbias)
for nm, mk in [('non-outlier records', nout_mask), ('outlier records', out_finite),
               ('  above outliers',     above_mask & np.isfinite(carrier_hbias)),
               ('  below outliers',     below_mask & np.isfinite(carrier_hbias))]:
    if mk.sum() == 0: continue
    v = carrier_hbias[mk]
    print(f'  {nm:<24}  median={np.median(v):+.3f}  P25={np.percentile(v,25):+.3f}  P75={np.percentile(v,75):+.3f}  P10={np.percentile(v,10):+.3f}  P90={np.percentile(v,90):+.3f}')

# === Cactus vs PG carrier composition at outlier records ===
# At outlier records, what fraction of carriers come from each class?
print(f'\n=== Cactus vs PG carrier composition at outliers ===')
# Among outliers: classify by carrier composition
total_ac = ac_c + ac_p
cact_frac = np.where(total_ac > 0, ac_c / total_ac, np.nan)
for nm, mk, n_total_class in [('all records',  np.isfinite(cact_frac), is_cactus.sum()/(is_cactus.sum()+is_pg.sum())),
                                ('non-outliers', nout_mask & np.isfinite(cact_frac), is_cactus.sum()/(is_cactus.sum()+is_pg.sum())),
                                ('outliers',     out_finite & np.isfinite(cact_frac), is_cactus.sum()/(is_cactus.sum()+is_pg.sum())),
                                ('  above',      above_mask & np.isfinite(cact_frac), is_cactus.sum()/(is_cactus.sum()+is_pg.sum())),
                                ('  below',      below_mask & np.isfinite(cact_frac), is_cactus.sum()/(is_cactus.sum()+is_pg.sum())) ]:
    if mk.sum() == 0: continue
    v = cact_frac[mk]
    print(f'  {nm:<20}  (n={mk.sum():>8,})  cactus fraction of carriers: median={np.median(v):.3f}  P25={np.percentile(v,25):.3f}  P75={np.percentile(v,75):.3f}  (uniform={n_total_class:.3f})')

# === Spot-check top 10 above + below outliers ===
print(f'\n=== Top 10 ABOVE outliers (arch3 >> recipe) ===')
print(f"{'idx':>10} {'pos':>10} {'ref':>4} {'alt':>4} {'recipe':>7} {'arch3':>7} {'Δ':>7} {'AC':>5} {'AC_c':>5} {'AC_p':>5} {'mean_h_car/avg':>16}")
above_idx = np.argsort(-(diff * above_mask))[:10]  # most positive deviations
for i in above_idx:
    if not above_mask[i]: continue
    print(f'  {i:>8,} {int(meta["pos"][i]):>10,} {str(meta["ref"][i])[:4]:>4} {str(meta["alt"][i])[:4]:>4} {recipe_af[i]:>7.3f} {arch3_af[i]:>7.3f} {diff[i]:>+7.3f} {int(ac[i]):>5} {int(ac_c[i]):>5} {int(ac_p[i]):>5} {carrier_hbias[i]+1:>16.3f}')

print(f'\n=== Top 10 BELOW outliers (arch3 << recipe) ===')
print(f"{'idx':>10} {'pos':>10} {'ref':>4} {'alt':>4} {'recipe':>7} {'arch3':>7} {'Δ':>7} {'AC':>5} {'AC_c':>5} {'AC_p':>5} {'mean_h_car/avg':>16}")
below_idx = np.argsort(diff * below_mask)[:10]  # most negative deviations
for i in below_idx:
    if not below_mask[i]: continue
    print(f'  {i:>8,} {int(meta["pos"][i]):>10,} {str(meta["ref"][i])[:4]:>4} {str(meta["alt"][i])[:4]:>4} {recipe_af[i]:>7.3f} {arch3_af[i]:>7.3f} {diff[i]:>+7.3f} {int(ac[i]):>5} {int(ac_c[i]):>5} {int(ac_p[i]):>5} {carrier_hbias[i]+1:>16.3f}')
PYEOF

echo
echo "[$(date)] DONE E2v2"

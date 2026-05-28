#!/usr/bin/env python3
"""
At Chr1 biallelic SNP records with F_MISSING_post >= 0.10, decompose every missing
PG cell into one of:
  n_hetmask_pg(r)   = AC_Het_pre(r)                  # 0/1 cells turned to ./.
  n_premissing_pg(r) = (306 - AN_pre)/2              # already ./. pre-mask (GQ<20, PanGenie no-call, cactus-only)
  n_cactus_missing(r) = (231 - AN_post) - n_hetmask_pg - n_premissing_pg  # cactus haploid `.`

Then check whether the hapFIRE-cactus_em residual scales with the het-mask share
of total missingness.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')

J = pd.read_csv(ROOT / 'scratch/hetmask_diag/joined_hetmask_chr1.tsv', sep='\t')
fm = pd.read_csv(
    ROOT / 'scratch/hetmask_diag/panel231_chr1_snp_fmissing.tsv',
    sep='\t', header=None,
    names=['chrom','pos','ref','alt','AC_post','AN_post','F_MISSING_post']
)
for c in ('AC_post','AN_post','F_MISSING_post'):
    fm[c] = pd.to_numeric(fm[c], errors='coerce')

J = J.merge(fm, on=['chrom','pos','ref','alt'], how='inner')
J['residual_hf'] = J['hapfire_af'] - J['alt_freq']
J['abs_res']    = J['residual_hf'].abs()

# Decompose missing-cell count per record
# AN_pre comes from the pre-mask PG file (153 diploid samples, max AN = 306)
J['n_hetmask_pg']     = J['AC_Het'].astype(int)
J['n_premissing_pg']  = ((306 - J['AN_pre']) / 2).astype(int)
J['n_missing_total']  = (231 - J['AN_post']).astype(int)         # all-panel missing (haploid)
J['n_cactus_missing'] = J['n_missing_total'] - J['n_hetmask_pg'] - J['n_premissing_pg']

# Sanity: n_cactus_missing should be >= 0
neg = (J['n_cactus_missing'] < 0).sum()
print(f'rows with n_cactus_missing < 0 (math sanity): {neg:,}  ({100*neg/len(J):.4f}%)')
# If negative ever appears it would mean post-mask had MORE missing on PG side than expected;
# this can happen at the AC=0 cleanup boundary (rare).

# Restrict to high-missing records: F_MISSING_post >= 0.10
HI = J[J['F_MISSING_post'] >= 0.10].copy()
print(f'\nrecords with F_MISSING_post >= 0.10 : {len(HI):,}')
print(f'of which n_missing_total > 0       : {(HI["n_missing_total"]>0).sum():,}')

# For each record: what fraction of the total missingness is from each bucket?
HI['share_hetmask']    = HI['n_hetmask_pg']    / HI['n_missing_total'].clip(lower=1)
HI['share_premissing'] = HI['n_premissing_pg'] / HI['n_missing_total'].clip(lower=1)
HI['share_cactus']     = HI['n_cactus_missing']/ HI['n_missing_total'].clip(lower=1)

print('\nMean missingness-source composition at F_MISSING_post >= 0.10 records:')
print(f'  share_hetmask     : {HI["share_hetmask"].mean():.3f}  (mean per-record fraction)')
print(f'  share_premissing  : {HI["share_premissing"].mean():.3f}')
print(f'  share_cactus      : {HI["share_cactus"].mean():.3f}')
print(f'  sum (should ~1)   : '
      f'{(HI["share_hetmask"]+HI["share_premissing"]+HI["share_cactus"]).mean():.3f}')

# Mean missing cells per bucket
print('\nMean missing cells per bucket at F_MISSING_post >= 0.10 records:')
print(f'  n_hetmask_pg      mean = {HI["n_hetmask_pg"].mean():.2f}')
print(f'  n_premissing_pg   mean = {HI["n_premissing_pg"].mean():.2f}')
print(f'  n_cactus_missing  mean = {HI["n_cactus_missing"].mean():.2f}')
print(f'  n_missing_total   mean = {HI["n_missing_total"].mean():.2f}')

# Correlations: residual vs each bucket COUNT (this is what biases the projection)
print('\nResidual vs missing-cell counts (high-missing records only):')
for col in ['n_hetmask_pg', 'n_premissing_pg', 'n_cactus_missing', 'n_missing_total']:
    rp = pearsonr(HI[col], HI['residual_hf'])
    rs = spearmanr(HI[col], HI['residual_hf'])
    print(f'  {col:18s}  Pearson r={rp.statistic:+.4f}  Spearman r={rs.statistic:+.4f}')

# Stratify by het_mask_share bins; check mean_bias + MAE
HI['het_share_bin'] = pd.cut(
    HI['share_hetmask'], bins=[-1e-9, 0.05, 0.2, 0.5, 0.8, 1.01],
    labels=['<0.05','0.05-0.2','0.2-0.5','0.5-0.8','>=0.8'],
)
print('\nResidual (hapFIRE-cactus_em) stratified by HET-MASK share of total missing:')
print(HI.groupby('het_share_bin', observed=True).agg(
    n=('residual_hf','size'),
    mean_residual=('residual_hf','mean'),
    median_residual=('residual_hf','median'),
    mean_abs_res=('abs_res','mean'),
    out05=('abs_res', lambda v: (v>0.05).mean()*100),
    out10=('abs_res', lambda v: (v>0.10).mean()*100),
).round(5).to_string())

HI['premiss_share_bin'] = pd.cut(
    HI['share_premissing'], bins=[-1e-9, 0.05, 0.2, 0.5, 0.8, 1.01],
    labels=['<0.05','0.05-0.2','0.2-0.5','0.5-0.8','>=0.8'],
)
print('\nResidual stratified by PRE-EXISTING-MISSING share:')
print(HI.groupby('premiss_share_bin', observed=True).agg(
    n=('residual_hf','size'),
    mean_residual=('residual_hf','mean'),
    median_residual=('residual_hf','median'),
    mean_abs_res=('abs_res','mean'),
    out05=('abs_res', lambda v: (v>0.05).mean()*100),
    out10=('abs_res', lambda v: (v>0.10).mean()*100),
).round(5).to_string())

HI['cactus_share_bin'] = pd.cut(
    HI['share_cactus'], bins=[-1e-9, 0.05, 0.2, 0.5, 0.8, 1.01],
    labels=['<0.05','0.05-0.2','0.2-0.5','0.5-0.8','>=0.8'],
)
print('\nResidual stratified by CACTUS-MISSING share:')
print(HI.groupby('cactus_share_bin', observed=True).agg(
    n=('residual_hf','size'),
    mean_residual=('residual_hf','mean'),
    median_residual=('residual_hf','median'),
    mean_abs_res=('abs_res','mean'),
    out05=('abs_res', lambda v: (v>0.05).mean()*100),
    out10=('abs_res', lambda v: (v>0.10).mean()*100),
).round(5).to_string())

# Multivariate regression: residual ~ n_hetmask + n_premissing + n_cactus
X = HI[['n_hetmask_pg','n_premissing_pg','n_cactus_missing']].values.astype(float)
y = HI['residual_hf'].values.astype(float)
X1 = np.hstack([X, np.ones((len(X),1))])
b, *_ = np.linalg.lstsq(X1, y, rcond=None)
yhat = X1 @ b
r2 = 1 - ((y - yhat)**2).sum() / ((y - y.mean())**2).sum()
print(f'\nMultivariate OLS  residual = a*n_hetmask + b*n_premissing + c*n_cactus + d')
print(f'  coef n_hetmask_pg    = {b[0]:+.6f}  (per cell)')
print(f'  coef n_premissing_pg = {b[1]:+.6f}')
print(f'  coef n_cactus_missing= {b[2]:+.6f}')
print(f'  intercept            = {b[3]:+.6f}')
print(f'  R^2                  = {r2:.4f}')

# Decompose: at the mean per-record, how much of residual is each bucket?
mh = HI['n_hetmask_pg'].mean()
mp = HI['n_premissing_pg'].mean()
mc = HI['n_cactus_missing'].mean()
print(f'\nAt mean per-record cell counts (n_h={mh:.1f}, n_p={mp:.1f}, n_c={mc:.1f}):')
print(f'  contribution het-mask    : {b[0]*mh:+.5f}')
print(f'  contribution pre-missing : {b[1]*mp:+.5f}')
print(f'  contribution cactus      : {b[2]*mc:+.5f}')
print(f'  intercept                : {b[3]:+.5f}')
print(f'  total predicted residual : {b[0]*mh+b[1]*mp+b[2]*mc+b[3]:+.5f}')
print(f'  actual mean residual     : {HI["residual_hf"].mean():+.5f}')

# Save the high-missing decomposition
out = ROOT / 'scratch/hetmask_diag/highmissing_decomposed.tsv'
HI[['chrom','pos','ref','alt','alt_freq','hapfire_af','residual_hf',
    'F_MISSING_post','AN_post',
    'n_hetmask_pg','n_premissing_pg','n_cactus_missing','n_missing_total',
    'share_hetmask','share_premissing','share_cactus']].to_csv(out, sep='\t', index=False)
print(f'\nWrote {out}')

# --- Plots ---
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for ax, share_col, title, color in [
    (axes[0], 'share_hetmask',    'het-mask share of missing',    'firebrick'),
    (axes[1], 'share_premissing', 'pre-existing share of missing','steelblue'),
    (axes[2], 'share_cactus',     'cactus-haploid share of missing','seagreen'),
]:
    # binned mean residual
    bins = np.arange(0, 1.05, 0.1)
    HI['bin'] = pd.cut(HI[share_col], bins=bins, include_lowest=True)
    g = HI.groupby('bin', observed=True)['residual_hf'].agg(['mean','count']).reset_index()
    centers = bins[:-1] + 0.05
    centers = centers[:len(g)]
    ax.bar(centers, g['mean'], width=0.08, color=color, alpha=0.85)
    for i,(m,n) in enumerate(zip(g['mean'], g['count'])):
        ax.text(centers[i], m + (0.001 if m>=0 else -0.001),
                f'n={n:,}', ha='center', va='bottom' if m>=0 else 'top', fontsize=7,
                rotation=90)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xlabel(title)
    ax.set_ylabel('mean (hapFIRE − cactus_em)')
    ax.set_title(f'Residual vs {title}\n(F_MISSING_post >= 0.10, n={len(HI):,})')
plt.tight_layout()
plt.savefig(ROOT / 'plots/MISSING_SOURCE_DECOMP_v3qc_v3_chr1.png', dpi=140, bbox_inches='tight')
print('Wrote plot.')

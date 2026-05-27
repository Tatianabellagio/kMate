#!/usr/bin/env python3
"""
Extended diagnostic on the joined hapFIRE / cactus_em / pre-mask-PG table.

1. Stratify residual by (het_rate, F_MISSING) jointly.
2. Outlier accounting: of records with |residual|>0.05, what fraction has
   het_rate>0 (potentially het-mask) vs het_rate==0 (other mechanism)?
3. Sign check: among records with het_rate>0, is hapFIRE-cactus_em
   systematically positive (consistent with cactus_em losing ALT signal)?
4. Toy reconstruction: if het-mask is the cause, the SHIFT applied is roughly
   -0.5 * (AC_Het / AN_pre) on the PG-side AF. Compare to observed shift.
"""
import numpy as np
import pandas as pd
from pathlib import Path

J = pd.read_csv('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/hetmask_diag/joined_hetmask_chr1.tsv',
                sep='\t')
print(f'rows: {len(J):,}')

# F_MISSING bins as in PIPELINE_STATE
fbins = [-1e-9, 0.01, 0.05, 0.1, 0.2, 0.5, 1.01]
flabels = ['<0.01','0.01-0.05','0.05-0.1','0.1-0.2','0.2-0.5','>=0.5']
J['fmiss_bin_pre'] = pd.cut(J['F_MISSING_pre'], bins=fbins, labels=flabels)

hbins = [-1e-9, 0.0, 0.01, 0.05, 0.5, 1.0]
hlabels = ['=0','(0,0.01]','(0.01,0.05]','(0.05,0.5]','(0.5,1]']
J['het_bin'] = pd.cut(J['het_rate'], bins=hbins, labels=hlabels)

print('\n=== Joint cross-tab: mean (hapFIRE - cactus_em) ===')
ct_mean = J.pivot_table(index='het_bin', columns='fmiss_bin_pre',
                       values='residual', aggfunc='mean', observed=True).round(4)
print(ct_mean.to_string())

print('\n=== Joint cross-tab: cell count n ===')
ct_n = J.pivot_table(index='het_bin', columns='fmiss_bin_pre',
                    values='residual', aggfunc='size', observed=True)
print(ct_n.to_string())

# Outlier accounting
print('\n=== Outliers (|residual| > 0.05) split by het_rate==0 vs >0 ===')
J['abs_res'] = J['residual'].abs()
out = J[J['abs_res'] > 0.05].copy()
print(f'    n outliers: {len(out):,}  ({100*len(out)/len(J):.2f}% of all rows)')
print(f'    of which het_rate == 0  : {(out["het_rate"]==0).sum():,} '
      f'({100*(out["het_rate"]==0).mean():.1f}%)')
print(f'    of which het_rate >  0  : {(out["het_rate"]>0).sum():,} '
      f'({100*(out["het_rate"]>0).mean():.1f}%)')

print('\n=== Outliers (|residual| > 0.1) — likely PG-genotype-source disagreement ===')
out10 = J[J['abs_res'] > 0.1].copy()
print(f'    n outliers: {len(out10):,}  ({100*len(out10)/len(J):.2f}% of all rows)')
print(f'    of which het_rate == 0  : {(out10["het_rate"]==0).sum():,} '
      f'({100*(out10["het_rate"]==0).mean():.1f}%)')
print(f'    of which pre-mask PG AF == 0 (PG sees no ALT here): '
      f'{(out10["pg_af_pre"]==0).sum():,} ({100*(out10["pg_af_pre"]==0).mean():.1f}%)')
print(f'    of which pre-mask PG AF >= 0.95: '
      f'{(out10["pg_af_pre"]>=0.95).sum():,} ({100*(out10["pg_af_pre"]>=0.95).mean():.1f}%)')

# Toy reconstruction: predicted shift from het-mask alone
# Pre-mask PG AC contributes equally per allele. Het cells contribute 1 ALT allele
# each (out of 2). Masking turns them into ./. so we lose:
#   AC_lost   = AC_Het    (one ALT allele per het)
#   AN_lost   = 2*AC_Het  (the whole genotype goes to missing)
# Post-mask AF (PG-side, if uniform h on PG) =
#   (AC_pre - AC_Het) / (AN_pre - 2*AC_Het)
# Pre-mask AF (PG-side) = AC_pre / AN_pre
# Δ_pg_af = pre - post  (positive = mask reduces PG-side AF) -- this is the
# predicted "missing ALT mass" cactus_em loses on the PG side ONLY.
J['AC_pre'] = J['pg_af_pre'] * J['AN_pre']  # reconstructed
J['post_pg_af'] = np.where(
    J['AN_pre'] - 2*J['AC_Het'] > 0,
    (J['AC_pre'] - J['AC_Het']) / (J['AN_pre'] - 2*J['AC_Het']),
    np.nan
)
J['delta_pg_af_pred'] = J['pg_af_pre'] - J['post_pg_af']   # >=0

# If the residual were entirely driven by the het-mask PG-side shift, we'd expect
# residual ~ (151/231) * delta_pg_af_pred  (PG founders are 151 of 231 mass under uniform h)
PG_FRACTION = 151.0/231.0
J['expected_residual_from_hetmask'] = PG_FRACTION * J['delta_pg_af_pred']

# Compare on records with het_rate > 0
hh = J[(J['het_rate'] > 0) & (J['het_rate'] <= 0.1)].copy()
print(f'\n=== Het-mask reconstruction (rows with 0 < het_rate <= 0.1, n={len(hh):,}) ===')
print('    mean observed residual                   :',
      f'{hh["residual"].mean():+.5f}')
print('    mean predicted residual (PG-share shift) :',
      f'{hh["expected_residual_from_hetmask"].mean():+.5f}')
print('    ratio observed/predicted                 :',
      f'{hh["residual"].mean() / hh["expected_residual_from_hetmask"].mean():.3f}')

# Per-bin observed vs predicted
hh['hbin'] = pd.cut(hh['het_rate'], bins=[0, 0.005, 0.01, 0.02, 0.05, 0.1],
                    labels=['(0,0.005]','(0.005,0.01]','(0.01,0.02]','(0.02,0.05]','(0.05,0.1]'])
recon = hh.groupby('hbin', observed=True).agg(
    n=('residual','size'),
    obs_mean=('residual','mean'),
    pred_mean=('expected_residual_from_hetmask','mean'),
    obs_med=('residual','median'),
    pred_med=('expected_residual_from_hetmask','median'),
).round(5)
recon['obs_over_pred'] = (recon['obs_mean'] / recon['pred_mean']).round(3)
print('\n    Per-het-bin observed-vs-predicted residual:')
print(recon.to_string())

# Variance explained by het-mask under the simple model
# Residual variance: total vs residual-of-(observed - predicted)
ss_total = (J['residual']**2).sum()
ss_after = ((J['residual'] - J['expected_residual_from_hetmask'].fillna(0))**2).sum()
print(f'\n=== Variance explained by het-mask model ===')
print(f'    sum(residual^2)                        = {ss_total:.1f}')
print(f'    sum((residual - predicted_hetmask)^2)  = {ss_after:.1f}')
print(f'    R^2 explained by het-mask model        = {1 - ss_after/ss_total:.4f}')

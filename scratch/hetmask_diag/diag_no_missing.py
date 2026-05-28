#!/usr/bin/env python3
"""
Restrict the hapFIRE vs cactus_em comparison to records where the v3qc_v3 panel
has NO missing data (AN_post == 231, F_MISSING_post == 0).

Compares overall agreement metrics on:
  (a) all biallelic Chr1 SNPs (the 517k baseline)
  (b) subset with F_MISSING_post == 0       (truly complete panel)
  (c) subset with F_MISSING_post < 0.01    (~near-complete; matches PIPELINE_STATE §5 bin)

Then re-runs the het_rate residual breakdown within (b).
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')

J = pd.read_csv(ROOT / 'scratch/hetmask_diag/joined_hetmask_chr1.tsv', sep='\t')

# Load F_MISSING from merged 231 haploid VCF
fm = pd.read_csv(
    ROOT / 'scratch/hetmask_diag/panel231_chr1_snp_fmissing.tsv',
    sep='\t', header=None,
    names=['chrom','pos','ref','alt','AC_post','AN_post','F_MISSING_post']
)
for c in ('AC_post','AN_post','F_MISSING_post'):
    fm[c] = pd.to_numeric(fm[c], errors='coerce')

before = len(J)
J = J.merge(fm, on=['chrom','pos','ref','alt'], how='inner')
print(f'Joined (before F_MISSING) : {before:,}')
print(f'Joined (after  F_MISSING) : {len(J):,}')

J['residual'] = J['hapfire_af'] - J['alt_freq']
J['abs_res']  = J['residual'].abs()

def metrics(df, label):
    r2 = pearsonr(df['hapfire_af'], df['alt_freq']).statistic**2
    rsp = spearmanr(df['hapfire_af'], df['alt_freq']).statistic
    mae = df['abs_res'].mean()
    rmse = np.sqrt((df['residual']**2).mean())
    mean_bias = df['residual'].mean()
    outliers_05 = (df['abs_res'] > 0.05).mean() * 100
    outliers_10 = (df['abs_res'] > 0.10).mean() * 100
    return {
        'n': len(df), 'R2': r2, 'spearman': rsp,
        'MAE': mae, 'RMSE': rmse, 'mean_bias': mean_bias,
        '%|d|>0.05': outliers_05, '%|d|>0.10': outliers_10,
    }

scenarios = {
    'all (joined)'           : J,
    'F_MISSING_post == 0'    : J[J['F_MISSING_post'] == 0],
    'F_MISSING_post < 0.01'  : J[J['F_MISSING_post'] < 0.01],
    'F_MISSING_post < 0.05'  : J[J['F_MISSING_post'] < 0.05],
    'F_MISSING_post >= 0.1'  : J[J['F_MISSING_post'] >= 0.1],
}

rows = []
for k,d in scenarios.items():
    m = metrics(d, k)
    m['scenario'] = k
    rows.append(m)
out = pd.DataFrame(rows).set_index('scenario')
print('\n=== hapFIRE vs cactus_em agreement under panel-missingness filters ===')
cols = ['n','R2','spearman','MAE','RMSE','mean_bias','%|d|>0.05','%|d|>0.10']
print(out[cols].round(5).to_string())

# Re-do the het_rate breakdown within the no-missing subset
sub = J[J['F_MISSING_post'] == 0].copy()
bins   = [-1e-9, 0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
labels = ['=0', '(0,0.005]', '(0.005,0.01]', '(0.01,0.02]', '(0.02,0.05]',
          '(0.05,0.1]', '(0.1,0.2]', '(0.2,0.5]', '(0.5,1.0]']
sub['het_bin'] = pd.cut(sub['het_rate'], bins=bins, labels=labels)
summary = sub.groupby('het_bin', observed=True).agg(
    n=('residual','size'),
    mean_residual=('residual','mean'),
    median_residual=('residual','median'),
    mean_abs_res=('abs_res','mean'),
    mean_hapfire=('hapfire_af','mean'),
    mean_cactus_em=('alt_freq','mean'),
    mean_pg_af_pre=('pg_af_pre','mean'),
).round(5)
print('\n=== Residual vs pre-mask PG het_rate, restricted to F_MISSING_post==0 ===')
print(summary.to_string())

# Pearson + Spearman within the clean subset
r_p = pearsonr(sub['het_rate'], sub['residual'])
r_s = spearmanr(sub['het_rate'], sub['residual'])
print(f'\nPearson  residual vs het_rate (no-missing subset): r={r_p.statistic:+.4f}  p={r_p.pvalue:.2e}')
print(f'Spearman residual vs het_rate (no-missing subset): r={r_s.statistic:+.4f}  p={r_s.pvalue:.2e}')

x = sub['het_rate'].values
y = sub['residual'].values
A = np.vstack([x, np.ones_like(x)]).T
slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
print(f'\nresidual = {slope:+.4f} * het_rate + {intercept:+.5f}  (within F_MISSING_post==0)')

# Save the no-missing subset for plotting
out_no_missing = ROOT / 'scratch/hetmask_diag/joined_no_missing.tsv'
sub.to_csv(out_no_missing, sep='\t', index=False)
print(f'\nWrote {out_no_missing}  ({len(sub):,} rows)')

# Also produce a (residual vs het_rate) bin-mean for plotting
plot_df = sub.groupby('het_bin', observed=True).agg(
    n=('residual','size'),
    mean_residual=('residual','mean'),
    sem_residual=('residual', lambda v: v.std(ddof=1)/np.sqrt(len(v))),
    mean_abs_res=('abs_res','mean'),
).reset_index()
plot_df.to_csv(ROOT / 'scratch/hetmask_diag/bin_means_no_missing.tsv', sep='\t', index=False)
print('Wrote bin_means_no_missing.tsv')
print('\nBin means for plotting (no-missing subset):')
print(plot_df.to_string())

# Outlier accounting on no-missing subset
print('\n=== Outliers (|residual| > 0.05/0.10) in no-missing subset ===')
n = len(sub)
n05 = (sub['abs_res'] > 0.05).sum()
n10 = (sub['abs_res'] > 0.10).sum()
print(f'    |residual| > 0.05: {n05:,} ({100*n05/n:.2f}%)')
print(f'    |residual| > 0.10: {n10:,} ({100*n10/n:.2f}%)')
# Split by het_rate
o5 = sub[sub['abs_res'] > 0.05]
o10 = sub[sub['abs_res'] > 0.10]
print(f'    Among |d|>0.05: het_rate==0 fraction = {(o5["het_rate"]==0).mean()*100:.1f}%')
print(f'    Among |d|>0.10: het_rate==0 fraction = {(o10["het_rate"]==0).mean()*100:.1f}%')

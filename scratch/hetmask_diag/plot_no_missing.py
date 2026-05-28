#!/usr/bin/env python3
"""Plot hapFIRE vs cactus_em agreement, all-records vs no-missing subset."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

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
J['residual'] = J['hapfire_af'] - J['alt_freq']
J['abs_res']  = J['residual'].abs()

sub_no  = J[J['F_MISSING_post'] == 0].copy()
sub_lo  = J[(J['F_MISSING_post'] > 0) & (J['F_MISSING_post'] < 0.05)].copy()
sub_hi  = J[J['F_MISSING_post'] >= 0.1].copy()

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

def scatter_panel(ax, df, title):
    n = len(df)
    if n > 30000:
        d = df.sample(30000, random_state=0)
    else:
        d = df
    ax.scatter(d['alt_freq'], d['hapfire_af'], s=2, alpha=0.18, edgecolor='none')
    ax.plot([0,1],[0,1], color='red', lw=1, ls='--')
    mae = df['abs_res'].mean()
    r2  = np.corrcoef(df['hapfire_af'], df['alt_freq'])[0,1]**2
    out5 = (df['abs_res']>0.05).mean()*100
    ax.set_xlabel('cactus_em alt_freq (v3qc_v3 mixed-loose)')
    ax.set_ylabel('hapFIRE alt_freq')
    ax.set_title(f'{title}\n n={n:,}  MAE={mae:.4f}  R²={r2:.4f}  %|d|>0.05={out5:.2f}')
    ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.set_aspect('equal')

scatter_panel(axes[0,0], J,       'All biallelic Chr1 SNPs (joined)')
scatter_panel(axes[0,1], sub_no,  'F_MISSING_post = 0 (no missing in panel)')
scatter_panel(axes[0,2], sub_hi,  'F_MISSING_post ≥ 0.10 (high-missing tail)')

# Residual histograms
def res_hist(ax, df, title, color='steelblue'):
    r = df['residual'].values
    rng = (-0.15, 0.15)
    ax.hist(r, bins=120, range=rng, color=color, alpha=0.85, density=True)
    ax.axvline(0, color='k', lw=0.5)
    ax.set_xlabel('hapFIRE - cactus_em')
    ax.set_ylabel('density')
    mae = np.abs(r).mean(); med = np.median(r); mean = r.mean()
    ax.set_title(f'{title}\n n={len(df):,}  mean={mean:+.4f}  median={med:+.4f}  MAE={mae:.4f}')
    ax.set_xlim(rng)

res_hist(axes[1,0], J,      'All joined')
res_hist(axes[1,1], sub_no, 'F_MISSING_post = 0', color='seagreen')
res_hist(axes[1,2], sub_hi, 'F_MISSING_post >= 0.10', color='firebrick')

plt.suptitle('hapFIRE vs cactus_em agreement on SEEDMIX_S1 Chr1 biallelic SNPs\n'
             '(v3qc_v3 mixed-loose cactus_em)', y=1.00)
plt.tight_layout()
fig_path = ROOT / 'scratch/hetmask_diag/agreement_no_missing.png'
plt.savefig(fig_path, dpi=140, bbox_inches='tight')
print(f'Wrote {fig_path}')

# Also a clean F_MISSING_post-binned summary plot
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 4))
bins   = [-1e-9, 0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.01]
labels = ['=0', '(0,0.01]','(0.01,0.05]','(0.05,0.1]','(0.1,0.2]','(0.2,0.5]','>=0.5']
J['fmiss_bin_post'] = pd.cut(J['F_MISSING_post'], bins=bins, labels=labels)
binsum = J.groupby('fmiss_bin_post', observed=True).agg(
    n=('residual','size'),
    mae=('abs_res','mean'),
    mean_bias=('residual','mean'),
    out05=('abs_res', lambda v: (v>0.05).mean()*100),
).reset_index()
print(binsum.to_string())

x = np.arange(len(binsum))
axes2[0].bar(x, binsum['mae'], color='steelblue')
axes2[0].set_xticks(x); axes2[0].set_xticklabels(binsum['fmiss_bin_post'], rotation=30)
axes2[0].set_ylabel('MAE'); axes2[0].set_title('MAE by F_MISSING_post bin')
for i,(v,n) in enumerate(zip(binsum['mae'], binsum['n'])):
    axes2[0].text(i, v, f' n={n:,}', rotation=90, va='bottom', ha='center', fontsize=8)

axes2[1].bar(x, binsum['out05'], color='firebrick')
axes2[1].set_xticks(x); axes2[1].set_xticklabels(binsum['fmiss_bin_post'], rotation=30)
axes2[1].set_ylabel('% records with |d| > 0.05'); axes2[1].set_title('Outlier rate by F_MISSING_post bin')
plt.tight_layout()
fig2_path = ROOT / 'scratch/hetmask_diag/agreement_by_fmiss.png'
plt.savefig(fig2_path, dpi=140, bbox_inches='tight')
print(f'Wrote {fig2_path}')

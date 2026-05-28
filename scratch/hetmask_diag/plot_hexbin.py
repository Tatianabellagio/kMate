#!/usr/bin/env python3
"""Hexbin version of hapFIRE vs cactus_em agreement plot."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
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

sub_no = J[J['F_MISSING_post'] == 0].copy()
sub_hi = J[J['F_MISSING_post'] >= 0.1].copy()

fig, axes = plt.subplots(2, 3, figsize=(16, 11))

def hex_panel(ax, df, title, gridsize=80):
    n = len(df)
    hb = ax.hexbin(df['alt_freq'], df['hapfire_af'],
                   gridsize=gridsize, extent=(0,1,0,1),
                   norm=LogNorm(vmin=1), mincnt=1, cmap='viridis')
    ax.plot([0,1],[0,1], color='red', lw=1.2, ls='--')
    mae = df['abs_res'].mean()
    r2  = np.corrcoef(df['hapfire_af'], df['alt_freq'])[0,1]**2
    out5 = (df['abs_res']>0.05).mean()*100
    out10 = (df['abs_res']>0.10).mean()*100
    ax.set_xlabel('cactus_em alt_freq (v3qc_v3 mixed-loose)')
    ax.set_ylabel('hapFIRE alt_freq')
    ax.set_title(f'{title}\n n={n:,}  MAE={mae:.4f}  R²={r2:.4f}  '
                 f'%|d|>0.05={out5:.2f}  %|d|>0.10={out10:.2f}')
    ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.set_aspect('equal')
    cb = fig.colorbar(hb, ax=ax, shrink=0.85)
    cb.set_label('records / hex (log)')

hex_panel(axes[0,0], J,       'All biallelic Chr1 SNPs (joined)')
hex_panel(axes[0,1], sub_no,  'F_MISSING_post = 0 (no missing in panel)')
hex_panel(axes[0,2], sub_hi,  'F_MISSING_post >= 0.10 (high-missing tail)')

def res_hist(ax, df, title, color='steelblue'):
    r = df['residual'].values
    rng = (-0.15, 0.15)
    ax.hist(r, bins=120, range=rng, color=color, alpha=0.85, density=True)
    ax.axvline(0, color='k', lw=0.5)
    ax.set_xlabel('hapFIRE - cactus_em')
    ax.set_ylabel('density')
    mae = np.abs(r).mean(); med = np.median(r); mean = r.mean()
    ax.set_title(f'{title}\n n={len(df):,}  mean={mean:+.4f}  '
                 f'median={med:+.4f}  MAE={mae:.4f}')
    ax.set_xlim(rng)

res_hist(axes[1,0], J,      'All joined')
res_hist(axes[1,1], sub_no, 'F_MISSING_post = 0', color='seagreen')
res_hist(axes[1,2], sub_hi, 'F_MISSING_post >= 0.10', color='firebrick')

plt.suptitle('hapFIRE vs cactus_em agreement on SEEDMIX_S1 Chr1 biallelic SNPs '
             '(v3qc_v3 mixed-loose)', y=1.00)
plt.tight_layout()
fig_path = ROOT / 'plots/HAPFIRE_VS_CACTUSEM_NO_MISSING_v3qc_v3.png'
plt.savefig(fig_path, dpi=140, bbox_inches='tight')
print(f'Wrote {fig_path}')

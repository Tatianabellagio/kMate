#!/usr/bin/env python3
"""Generate hexbin scatter plots of SAFER coord-aware results vs hapFIRE."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')
M = pd.read_csv(ROOT / 'scratch/f2_honest/coord_aware_safer.tsv', sep='\t')
M['cell'] = (M['high_fmiss'].astype(int)*2 + M['is_mixed_bubble']).map({
    0:'[a] low_F pure_SNP',
    1:'[b] low_F mixed_bub',
    2:'[c] high_F pure_SNP',
    3:'[d] high_F mixed_bub',
})
print(f'rows: {len(M):,}')

def hex_panel(ax, df, x_col, y_col, title, gridsize=80):
    n = len(df)
    if n == 0:
        ax.set_title(f'{title}\nempty')
        return
    hb = ax.hexbin(df[x_col], df[y_col], gridsize=gridsize, extent=(0,1,0,1),
                   norm=LogNorm(vmin=1), mincnt=1, cmap='viridis')
    ax.plot([0,1],[0,1], color='red', lw=1.2, ls='--')
    mae = (df[y_col] - df[x_col]).abs().mean()
    r2 = np.corrcoef(df[y_col], df[x_col])[0,1]**2 if len(df) > 1 else float('nan')
    out5 = ((df[y_col] - df[x_col]).abs() > 0.05).mean() * 100
    out10 = ((df[y_col] - df[x_col]).abs() > 0.10).mean() * 100
    ax.set_title(f'{title}\nn={n:,}  MAE={mae:.4f}  R²={r2:.4f}\n%|d|>0.05={out5:.2f}  %|d|>0.10={out10:.2f}')
    ax.set_aspect('equal'); ax.set_xlim(0,1); ax.set_ylim(0,1)
    return hb

# ─── Figure 1: overall (all Chr1 SNPs joined) — orig vs SAFER vs hapFIRE ───
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
hb = hex_panel(axes[0], M, 'alt_freq', 'hapfire_af',
               'cactus_em (orig) vs hapFIRE — ALL Chr1 joined SNPs')
axes[0].set_xlabel('cactus_em alt_freq (orig)'); axes[0].set_ylabel('hapFIRE alt_freq')
fig.colorbar(hb, ax=axes[0], shrink=0.85).set_label('records/hex (log)')
hb = hex_panel(axes[1], M, 'cactus_em_safer', 'hapfire_af',
               'cactus_em + SAFER vs hapFIRE — ALL Chr1 joined SNPs')
axes[1].set_xlabel('cactus_em + SAFER alt_freq'); axes[1].set_ylabel('hapFIRE alt_freq')
fig.colorbar(hb, ax=axes[1], shrink=0.85).set_label('records/hex (log)')
plt.tight_layout()
fp = ROOT / 'plots/SAFER_overall_v3qc_v3.png'
plt.savefig(fp, dpi=140, bbox_inches='tight')
print(f'Wrote {fp}')

# ─── Figure 2: per-cell SAFER vs hapFIRE (a, b, c, d) ───
fig, axes = plt.subplots(2, 4, figsize=(18, 9))
for i, (cell_name, sub) in enumerate(M.groupby('cell')):
    # Row 0: orig
    hex_panel(axes[0, i], sub, 'alt_freq', 'hapfire_af',
              f'{cell_name}\nORIGINAL')
    axes[0, i].set_xlabel('cactus_em (orig)'); axes[0, i].set_ylabel('hapFIRE')
    # Row 1: SAFER
    hex_panel(axes[1, i], sub, 'cactus_em_safer', 'hapfire_af',
              f'{cell_name}\nSAFER')
    axes[1, i].set_xlabel('cactus_em + SAFER'); axes[1, i].set_ylabel('hapFIRE')
plt.suptitle('SAFER coord-aware: per-cell comparison (Chr1 SNPs)', y=1.00, fontsize=14)
plt.tight_layout()
fp = ROOT / 'plots/SAFER_per_cell_v3qc_v3.png'
plt.savefig(fp, dpi=140, bbox_inches='tight')
print(f'Wrote {fp}')

# ─── Figure 3: residual histograms ───
fig, axes = plt.subplots(2, 4, figsize=(18, 7))
for i, (cell_name, sub) in enumerate(M.groupby('cell')):
    rng = (-0.2, 0.2)
    bins = 80
    r_o = sub['hapfire_af'] - sub['alt_freq']
    r_s = sub['hapfire_af'] - sub['cactus_em_safer']
    axes[0, i].hist(r_o, bins=bins, range=rng, color='steelblue', alpha=0.85, density=True)
    axes[0, i].axvline(0, color='k', lw=0.5)
    axes[0, i].set_title(f'{cell_name}\norig: mean={r_o.mean():+.4f}, MAE={r_o.abs().mean():.4f}')
    axes[0, i].set_xlim(rng); axes[0, i].set_xlabel('hapFIRE − cactus_em (orig)')
    axes[1, i].hist(r_s, bins=bins, range=rng, color='seagreen', alpha=0.85, density=True)
    axes[1, i].axvline(0, color='k', lw=0.5)
    axes[1, i].set_title(f'{cell_name}\nSAFER: mean={r_s.mean():+.4f}, MAE={r_s.abs().mean():.4f}')
    axes[1, i].set_xlim(rng); axes[1, i].set_xlabel('hapFIRE − (cactus_em+SAFER)')
plt.suptitle('Residual histograms: orig (top) vs SAFER (bottom)', y=1.01, fontsize=14)
plt.tight_layout()
fp = ROOT / 'plots/SAFER_residual_hist_v3qc_v3.png'
plt.savefig(fp, dpi=140, bbox_inches='tight')
print(f'Wrote {fp}')

# ─── Figure 4: focus on cell [d] (the one that matters) ───
D = M[M['cell']=='[d] high_F mixed_bub']
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
hex_panel(axes[0], D, 'alt_freq', 'hapfire_af',
          'cell [d] ORIGINAL', gridsize=60)
axes[0].set_xlabel('cactus_em (orig)'); axes[0].set_ylabel('hapFIRE')
hex_panel(axes[1], D, 'cactus_em_safer', 'hapfire_af',
          'cell [d] SAFER', gridsize=60)
axes[1].set_xlabel('cactus_em + SAFER'); axes[1].set_ylabel('hapFIRE')
if 'cem_unsafe' in D.columns:
    hex_panel(axes[2], D, 'cem_unsafe', 'hapfire_af',
              'cell [d] UNSAFE (for reference)', gridsize=60)
    axes[2].set_xlabel('cactus_em + UNSAFE'); axes[2].set_ylabel('hapFIRE')
plt.suptitle(f'Cell [d] (high F_MISSING + mixed bubble, n={len(D):,})', y=1.02, fontsize=14)
plt.tight_layout()
fp = ROOT / 'plots/SAFER_cell_d_zoom_v3qc_v3.png'
plt.savefig(fp, dpi=140, bbox_inches='tight')
print(f'Wrote {fp}')

print('\nDone. Plots at plots/SAFER_*.png')

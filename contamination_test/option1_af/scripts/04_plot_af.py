#!/usr/bin/env python
"""
Visualize the AF analysis. Plots:
1. Scatter: obs_AF vs panel_AF per sample, with regression line.
2. Histogram: obs_AF distribution per sample (all sites with DP>=8).
3. Combined cross-sample histogram with sequencing-error reference.

If contamination existed at fraction c, obs_AF would track c * panel_AF.
Slope = 0 + intercept = sequencing error => no contamination.
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path('/home/tbellagio/scratch/hapfire_sv/contamination_test/option1_af')
RESULTS = ROOT / 'results'

# --- load data ---
diag = pd.concat([
    pd.read_csv(p, sep='\t', dtype={'chrom': str})
    for p in sorted((ROOT / 'sites').glob('diag_chr*.tsv'))
], ignore_index=True)
diag['chrom_pos'] = diag['chrom'] + ':' + diag['pos'].astype(str)
diag = diag.set_index('chrom_pos')

samples = ['S1', 'S2', 'S3', 'S4']
af = {}
for s in samples:
    df = pd.read_csv(ROOT / 'af' / f'seedmix_{s}.af.tsv', sep='\t', dtype={'chrom': str})
    df['chrom_pos'] = df['chrom'] + ':' + df['pos'].astype(str)
    df = df.set_index('chrom_pos')
    af[s] = df.join(diag[['af_non231']], how='inner')
    print(f'{s}: {len(af[s]):,} sites')

# --- Plot 1: obs_AF vs panel_AF (one panel per sample) ---
fig, axes = plt.subplots(1, 4, figsize=(20, 4.5), sharex=True, sharey=True)
for ax, s in zip(axes, samples):
    df = af[s]
    df = df[df['dp'] >= 8]
    df = df[df['af_non231'] >= 0.05]  # focus on detectable bins
    ax.scatter(df['af_non231'], df['af'], s=4, alpha=0.3, color='#444444')

    # WLS regression
    x, y, w = df['af_non231'].values, df['af'].values, df['dp'].values.astype(float)
    W = np.sqrt(w)
    A = np.column_stack([x*W, W])
    sol, *_ = np.linalg.lstsq(A, y*W, rcond=None)
    c_hat, eps_hat = sol
    xs = np.linspace(0, df['af_non231'].max(), 50)
    ax.plot(xs, c_hat*xs + eps_hat, color='red', lw=2,
            label=f'ĉ = {c_hat:.4f}\nintercept = {eps_hat:.4f}')

    # reference lines for c = 1%, 5%, 10% contamination
    for c in (0.01, 0.05):
        ax.plot(xs, c*xs + 0.005, ls='--', alpha=0.4,
                label=f'c = {c*100:.0f}% (expected)')

    ax.axhline(0.005, color='gray', ls=':', alpha=0.6, label='0.5% seq error')
    ax.set_xlabel('panel AF in non-231')
    ax.set_ylabel('observed AF in BAM')
    ax.set_title(f'SEEDMIX_{s} (n={len(df):,} sites)')
    ax.legend(fontsize=7, loc='upper left')
    ax.set_ylim(-0.005, 0.06)
    ax.set_xlim(0.04, df['af_non231'].max() + 0.01)

fig.suptitle('Direct contamination test: obs_AF vs panel_AF at non-231-private SNPs', y=1.02)
plt.tight_layout()
plt.savefig(RESULTS / 'af_scatter.png', dpi=140, bbox_inches='tight')
plt.show()
print('saved af_scatter.png')

# --- Plot 2: histogram of obs_AF (all sites with DP>=8) ---
fig, ax = plt.subplots(figsize=(10, 5))
all_obs = np.concatenate([af[s][af[s]['dp']>=8]['af'].values for s in samples])
ax.hist(all_obs, bins=np.linspace(0, 0.05, 51), alpha=0.7, color='#444444')
ax.axvline(0.005, color='red', ls='--', label='0.5% seq error')
ax.axvline(np.median(all_obs), color='blue', ls='-', label=f'median = {np.median(all_obs):.4f}')
ax.set_xlabel('observed AF in seedmix BAM (at non-231-private sites)')
ax.set_ylabel('# (sample × site)')
ax.set_title(f'Distribution of obs_AF across {len(samples)} samples × {sum(len(af[s]) for s in samples):,} site-measurements')
ax.legend()
plt.tight_layout()
plt.savefig(RESULTS / 'af_hist.png', dpi=140, bbox_inches='tight')
plt.show()
print('saved af_hist.png')

# --- Per-bin obs_AF ---
fig, ax = plt.subplots(figsize=(9, 5))
bins = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.15]
labels = [f'{b*100:.1f}-{bs*100:.1f}%' for b, bs in zip(bins[:-1], bins[1:])]
for s in samples:
    df = af[s][af[s]['dp']>=8].copy()
    df['bin'] = pd.cut(df['af_non231'], bins=bins, labels=labels, include_lowest=True)
    means = df.groupby('bin', observed=True)['af'].mean()
    ax.plot(range(len(labels)), [means.get(l, np.nan) for l in labels], 'o-', label=s)

# expected lines for hypothetical contamination rates
for c in (0.001, 0.005, 0.01, 0.05, 0.10):
    expected = [c * (b + bn) / 2 + 0.005 for b, bn in zip(bins[:-1], bins[1:])]
    ax.plot(range(len(labels)), expected, ls='--', alpha=0.4, label=f'c = {c*100:.1f}% (expected)')

ax.axhline(0.005, color='gray', ls=':', label='0.5% seq error')
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=30, ha='right')
ax.set_xlabel('panel AF bin (non-231)')
ax.set_ylabel('mean observed AF in BAM')
ax.set_title('Mean obs_AF by panel AF bin — flat = no contamination')
ax.legend(fontsize=8, loc='upper left')
plt.tight_layout()
plt.savefig(RESULTS / 'af_by_bin.png', dpi=140, bbox_inches='tight')
plt.show()
print('saved af_by_bin.png')

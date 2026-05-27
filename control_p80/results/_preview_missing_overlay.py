"""Quick standalone preview of the 'color by missingness' scatter.
Iterate on the visualization here without rebuilding the whole notebook.

Usage:
    python _preview_missing_overlay.py
Outputs preview_missing_overlay.png in this directory.

Cells previewed: a small grid covering the most informative regimes for the
'fan from missingness' signal:
  - n50_g3 (fan-prone, lots of recomb)
  - n231_g0 (clean)
  - n50_g3_dom500 (selection — should show different residual structure)
For each: global and star2 side-by-side.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p80')
RESULTS = ROOT / 'results'
SIMS = ROOT / 'sims'
DATA = ROOT / 'data'

# Order: easiest → hardest, based on the headline-table MAE on the realized-pool truth.
REGIME_SIMDIR = {
    'n231_g0':        'cov10_n231_g0_s42_hotspots_p80_chr1',         # perfect mix, no recomb
    'n50_g0':         'cov10_n50_g0_s42_hotspots_p80_chr1',           # subset, no recomb
    'n231_g1':        'cov10_n231_g1_s42_hotspots_p80_chr1',          # perfect mix + 1 gen recomb
    'n50_g1':         'cov10_n50_g1_s42_hotspots_p80_chr1',           # subset + 1 gen recomb
    'n50_g3':         'cov10_n50_g3_s42_hotspots_p80_chr1',           # subset + 3 gen recomb
    'n50_g3_dom500':  'cov10_n50_g3_s42_hotspots_dom500_p80_chr1',    # recomb + selection (hardest)
}
METHODS = ['global', 'star2']

# Subsample size for plotting — keep responsive
SUBSAMPLE = 200_000

print('Loading meta + cn_var_called ...')
meta = np.load(DATA / 'cn_var_p80.meta.npz', allow_pickle=True)
N = len(meta['pos'])
cn_called = sp.load_npz(DATA / 'cn_var_p80.cn_var_called.npz')
F = cn_called.shape[0]
called_per_rec = np.asarray(cn_called.sum(axis=0)).ravel()
missing_frac = (F - called_per_rec) / F
print(f'  N={N:,}  F={F}  mean missing_frac={missing_frac.mean()*100:.2f}%')


def load_pair(method, regime):
    """Returns (truth, est, missing_frac) for one (method, regime) cell."""
    tsv = RESULTS / f'cactus_em_{method}' / regime / f'p80_{regime}_cov10_s42.tsv'
    truth_gz = SIMS / REGIME_SIMDIR[regime] / 'recomb_truth.tsv.gz'
    est = pd.read_csv(tsv, sep='\t', usecols=['alt_freq'])['alt_freq'].to_numpy()
    truth = pd.read_csv(truth_gz, sep='\t', usecols=['truth_af'])['truth_af'].to_numpy()
    return truth, est


GREY = '#888888'

def scatter_by_missing(ax, x, y, miss, title, s=15,
                       subsample=SUBSAMPLE,
                       alpha=0.5):
    """Single-layer plot: COLOR = local 2D-histogram density (viridis, log)."""
    x = np.asarray(x); y = np.asarray(y)
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]; y = y[finite]
    keep = (x > 0) | (y > 0)
    x = x[keep]; y = y[keep]
    if len(x) == 0:
        return None

    # 2D-histogram density lookup per point
    nbins = 100
    H, xedges, yedges = np.histogram2d(x, y, bins=nbins, range=[[0, 1], [0, 1]])
    ix = np.clip(np.digitize(x, xedges) - 1, 0, nbins - 1)
    iy = np.clip(np.digitize(y, yedges) - 1, 0, nbins - 1)
    d = H[ix, iy].astype(float)

    if subsample is not None and len(x) > subsample:
        np.random.seed(42)
        sel = np.random.choice(len(x), subsample, replace=False)
        x, y, d = x[sel], y[sel], d[sel]

    cmap = plt.get_cmap('viridis')
    norm = matplotlib.colors.LogNorm(vmin=max(d.min(), 1), vmax=max(d.max(), 1))
    # High-density drawn on top
    order = np.argsort(d)
    ax.scatter(x[order], y[order], c=d[order], s=s, cmap=cmap, norm=norm,
               alpha=alpha, edgecolors='none', rasterized=True)

    ax.plot([0, 1], [0, 1], '--', lw=0.8, alpha=0.7, color=GREY)
    err = y - x
    mae = float(np.abs(err).mean())
    n = int(np.isfinite(err).sum())
    # Subtitle ABOVE the plot, non-bold
    ax.set_title(f'{title}\nMAE = {mae:.4f}   n = {n:,}',
                 fontsize=10, color=GREY)

    # Padding so points near 0 / 1 aren't clipped at edge.
    ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
    ax.set_aspect('equal', adjustable='box')
    # Major ticks at 0.0, 0.2, ..., 1.0; grid lines emanate from them.
    import matplotlib.ticker as mticker
    ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    ax.set_xticks(ticks); ax.set_yticks(ticks)
    # Drop all spines + tick marks; keep tick labels; thin grid.
    for s_name in ('top', 'right', 'left', 'bottom'):
        ax.spines[s_name].set_visible(False)
    ax.tick_params(length=0, colors=GREY, which='both')
    ax.grid(True, color='#dddddd', linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)   # grid behind the scatter
    ax.xaxis.label.set_color(GREY)
    ax.yaxis.label.set_color(GREY)

    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    return sm


# Method display names — friendlier than "global" / "star2"
METHOD_LABELS = {
    'global': 'h chromosome',
    'star2':  'h 10kb window',
}

REGIMES = list(REGIME_SIMDIR.keys())   # 6 regimes, ordered easiest → hardest
# Square cells: figsize_w / ncols == figsize_h / nrows. With 4" per cell:
CELL = 4.0
fig, axes = plt.subplots(len(REGIMES), len(METHODS),
                        figsize=(CELL * len(METHODS) + 1, CELL * len(REGIMES) + 0.5),
                        sharex=True, sharey=True)
last_sm = None
for i, regime in enumerate(REGIMES):
    for j, method in enumerate(METHODS):
        ax = axes[i, j]
        truth, est = load_pair(method, regime)
        title = f'{METHOD_LABELS[method]}  ·  {regime}'
        sm = scatter_by_missing(ax, truth, est, missing_frac, title)
        if sm is not None:
            last_sm = sm
        if j == 0: ax.set_ylabel('Estimated AF')
        if i == len(REGIMES) - 1: ax.set_xlabel('True AF')

# Shared horizontal colorbar at the bottom (one for all subplots).
fig.subplots_adjust(bottom=0.08, top=0.96, hspace=0.28, wspace=0.12)
cbar_ax = fig.add_axes([0.22, 0.035, 0.56, 0.011])
cb = fig.colorbar(last_sm, cax=cbar_ax, orientation='horizontal',
                  label='Local density (log count of records)')
cb.outline.set_edgecolor(GREY)
cb.ax.xaxis.set_tick_params(color=GREY, labelcolor=GREY)
cb.ax.xaxis.label.set_color(GREY)
out = RESULTS / 'preview_missing_overlay.png'
plt.savefig(out, dpi=130, bbox_inches='tight')
print(f'saved: {out}')

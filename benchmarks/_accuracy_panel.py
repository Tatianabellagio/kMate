"""Canonical kMate accuracy-panel aesthetic — single source of truth.

The density-colored truth-vs-estimate scatter used across every kMate "final
results" panel (p231 / p80, SNP/SV/atomized, sims benchmarks). Lifted verbatim
from benchmarks/p231/results/FINAL_RESULTS_p231.ipynb so that scripts and
notebooks share ONE definition of the look — import `scatter_density` /
`grid_panel` instead of re-implementing, and every figure stays consistent.

Aesthetic: 100-bin histogram2d local density, 150k subsample, viridis + LogNorm,
points drawn densest-last, dashed grey identity line, "MAE = … n = …" subtitle,
grey (#888888) text/spines, horizontal "Local density (log count of records)"
colorbar, equal aspect on [0,1].
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42   # editable text in vector exports
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

GREY = "#888888"


def load_cell(truth_path, est_path, miss_thr=None, F=None, truth_col="truth_af"):
    """Load (truth, est) for one panel cell from a truth .tsv.gz and a kMate est
    .tsv (cols alt_freq, n_called). If miss_thr is set (with F = #founders), keep
    only records whose missing fraction (F - n_called)/F <= miss_thr, i.e.
    n_called >= (1-miss_thr)*F — drops the low-support (SV / `./.`-heavy) tail."""
    tr = pd.read_csv(truth_path, sep="\t", usecols=[truth_col])[truth_col].values.astype(np.float32)
    es = pd.read_csv(est_path, sep="\t", usecols=["alt_freq", "n_called"])
    est = es["alt_freq"].values.astype(np.float32)
    nc = es["n_called"].values
    if len(tr) != len(est):
        raise RuntimeError(f"len mismatch {len(tr)} vs {len(est)}: {est_path}")
    if miss_thr is not None and F:
        keep = nc >= (1.0 - miss_thr) * F
        tr, est = tr[keep], est[keep]
    return tr, est


def scatter_density(ax, x, y, title, s=15, alpha=0.5, subsample=150_000, seed=42):
    """Density-colored truth(x) vs estimate(y) scatter on `ax`. Returns the
    mappable (for a shared colorbar) or None if empty."""
    x = np.asarray(x, dtype=np.float32); y = np.asarray(y, dtype=np.float32)
    m = np.isfinite(x) & np.isfinite(y); x = x[m]; y = y[m]
    keep = (x > 0) | (y > 0); x = x[keep]; y = y[keep]
    if len(x) == 0:
        ax.set_title(f"{title}\n(no data)", fontsize=10, color=GREY); return None
    nb = 100
    H, xe, ye = np.histogram2d(x, y, bins=nb, range=[[0, 1], [0, 1]])
    ix = np.clip(np.digitize(x, xe) - 1, 0, nb - 1)
    iy = np.clip(np.digitize(y, ye) - 1, 0, nb - 1)
    d = H[ix, iy].astype(np.float32)
    if subsample and len(x) > subsample:
        np.random.seed(seed)
        sel = np.random.choice(len(x), subsample, replace=False)
        x, y, d = x[sel], y[sel], d[sel]
    cmap = plt.get_cmap("viridis")
    norm = LogNorm(vmin=max(d.min(), 1), vmax=max(d.max(), 1))
    order = np.argsort(d)
    sm = ax.scatter(x[order], y[order], c=d[order], s=s, cmap=cmap, norm=norm,
                    alpha=alpha, edgecolors="none", rasterized=True)
    ax.plot([0, 1], [0, 1], "--", lw=0.8, alpha=0.7, color=GREY)
    err = y - x; mae = float(np.abs(err).mean()); n = int(np.isfinite(err).sum())
    ax.set_title(f"{title}\nMAE = {mae:.4f}   n = {n:,}", fontsize=10, color=GREY)
    ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([0, .2, .4, .6, .8, 1.0]); ax.set_yticks([0, .2, .4, .6, .8, 1.0])
    for sp in ("top", "right", "left", "bottom"):
        ax.spines[sp].set_edgecolor(GREY)
    ax.tick_params(colors=GREY, labelcolor=GREY)
    return sm


def grid_panel(cells, suptitle, outfile, ncols=2, cell=4.0, dpi=130):
    """Render a grid of density scatters in the canonical aesthetic.

    cells: list of (title, truth_array, est_array) in row-major order. A `None`
    entry (or trailing shortfall) leaves that grid slot blank — use it to align
    a baseline panel against a comparison row.
    Lays them out in `ncols` columns with a shared bottom colorbar + suptitle.
    """
    n = len(cells); nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(cell * ncols + 1, cell * nrows + 0.5),
                             sharex=True, sharey=True, squeeze=False)
    last_sm = None
    for k, ax in enumerate(axes.flat):
        if k >= n or cells[k] is None:
            ax.axis("off"); continue
        title, truth, est = cells[k]
        sm = scatter_density(ax, truth, est, title=title)
        if sm is not None:
            last_sm = sm
        if k % ncols == 0:
            ax.set_ylabel("Estimated AF", color=GREY)
        if k // ncols == nrows - 1:
            ax.set_xlabel("True AF", color=GREY)
    fig.subplots_adjust(bottom=0.07, top=0.95, hspace=0.30, wspace=0.10)
    if last_sm is not None:
        cax = fig.add_axes([0.25, 0.03, 0.50, 0.01])
        cb = fig.colorbar(last_sm, cax=cax, orientation="horizontal",
                          label="Local density (log count of records)")
        cb.outline.set_edgecolor(GREY)
        cb.ax.xaxis.set_tick_params(color=GREY, labelcolor=GREY)
        cb.ax.xaxis.label.set_color(GREY)
    fig.suptitle(suptitle, y=0.985, fontsize=11, color=GREY)
    fig.savefig(outfile, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {outfile}")

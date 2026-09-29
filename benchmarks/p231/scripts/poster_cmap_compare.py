#!/usr/bin/env python3
"""One-off: compare candidate perceptual colormaps for the poster accuracy scatter.
Same data (n50_g1, miss<=0.90, SNP+indel+SV), small grid, one colormap per panel."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT = Path("/global/scratch/users/tbellg/kmate/benchmarks/p231")
RES, SIMS, PLOTS = ROOT/"results", ROOT/"sims", ROOT/"results"/"plots"
F, REG, MISS_THR = 231, "n50_g1", 0.90

# candidate maps; seaborn's mako/rocket/flare/crest registered if seaborn present
CANDIDATES = ["viridis", "plasma", "magma", "cividis", "mako", "rocket"]
try:
    import seaborn  # noqa: registers mako/rocket/flare/crest
except Exception:
    CANDIDATES = [c for c in CANDIDATES if c not in ("mako", "rocket", "flare", "crest")]

def subdir(reg):
    n, g = reg.split("_g"); return f"cov10_{n}_g{g}_s42_hotspots_p231_chr1"

tr = pd.read_csv(SIMS/subdir(REG)/"recomb_truth_raw.tsv.gz", sep="\t", usecols=["truth_af"])
es = pd.read_csv(RES/"kmate_global_filt2invu_raw"/REG/f"p231_filt2invu_raw_{REG}_cov10_s42.tsv",
                 sep="\t", usecols=["alt_freq", "n_called"])
truth = tr["truth_af"].values.astype(np.float32)
est   = es["alt_freq"].values.astype(np.float32)
miss  = ((F - es["n_called"].values) / F).astype(np.float32)
sel = (miss <= MISS_THR) & np.isfinite(truth) & np.isfinite(est)
x, y = truth[sel], est[sel]
keep = (x > 0) | (y > 0); x, y = x[keep], y[keep]

nb = 100
H, xe, ye = np.histogram2d(x, y, bins=nb, range=[[0, 1], [0, 1]])
ix = np.clip(np.digitize(x, xe)-1, 0, nb-1); iy = np.clip(np.digitize(y, ye)-1, 0, nb-1)
dens = H[ix, iy].astype(np.float32)
SUB = 150_000
rng = np.random.default_rng(42); pick = rng.choice(len(x), SUB, replace=False)
x, y, dens = x[pick], y[pick], dens[pick]
order = np.argsort(dens)
vmax = max(np.percentile(dens, 98), 2)
norm = LogNorm(vmin=max(dens.min(), 1), vmax=vmax)

ncol = 3; nrow = int(np.ceil(len(CANDIDATES)/ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4*ncol, 4*nrow))
for ax, cm in zip(axes.ravel(), CANDIDATES):
    ax.scatter(x[order], y[order], c=dens[order], s=6, cmap=cm, norm=norm,
               alpha=0.95, edgecolors="none", rasterized=True)
    ax.plot([0,1],[0,1],"--",lw=1,color="#444",alpha=.8)
    ax.set_xlim(-.02,1.02); ax.set_ylim(-.02,1.02); ax.set_aspect("equal","box")
    ax.set_title(cm, fontsize=14); ax.set_xticks([0,.5,1]); ax.set_yticks([0,.5,1])
for ax in axes.ravel()[len(CANDIDATES):]:
    ax.axis("off")
fig.suptitle(f"colormap candidates — {REG}, miss<=0.90 (SNP+indel+SV)", fontsize=15)
fig.tight_layout(rect=[0,0,1,0.97])
out = PLOTS/"poster_cmap_compare.png"
fig.savefig(out, dpi=140, bbox_inches="tight"); print("saved:", out, "| maps:", CANDIDATES)

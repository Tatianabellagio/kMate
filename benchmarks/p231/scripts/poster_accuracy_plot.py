#!/usr/bin/env python3
"""One-off: poster accuracy scatter for kMate.

Single panel, all variant classes together (SNP+indel+SV, raw arm), missing_frac<=0.90.
Density-colored truth-vs-estimate, matching the manuscript panel figures' convention.

Usage:
    poster_accuracy_plot.py [REGIME] [CMAP]
        REGIME : n50_g1 (default), n50_g0, n231_g0, ...
        CMAP   : viridis (default) | earth (poster earth-tone ramp)
"""
import sys, math
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
# keep text as EDITABLE TEXT (not outlined paths) in the vector outputs, so the
# labels open as words in Illustrator. fonttype 42 = embedded TrueType.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
try:
    import seaborn  # noqa: registers mako/rocket/flare/crest colormaps
except Exception:
    pass

ROOT = Path("/global/scratch/users/tbellg/kmate/benchmarks/p231")
RES, SIMS, PLOTS = ROOT/"results", ROOT/"sims", ROOT/"results"/"plots"
F = 231
REG  = sys.argv[1] if len(sys.argv) > 1 else "n50_g1"
CMAP = sys.argv[2] if len(sys.argv) > 2 else "mako"
MISS_THR = 0.50
GREY = "#888888"

# poster pie earth-tone palette (eye-sampled; tweak hex here to taste)
BLUE, ORANGE = "#1f3b5c", "#e07e37"   # navy + orange slices of the pie
CMAPS = {
    # full 7-tone earth ramp: dark teal (sparse) -> warm gold (dense)
    "earth": ["#2c4e63", "#6e5068", "#6e4a33", "#a84b34", "#e07b39", "#bfa94e", "#f0e6cf"],
    # two-hue: blue (sparse) -> orange (dense)
    "blueorange": [BLUE, ORANGE],
    # two-hue with a pale midpoint to avoid the muddy blue->orange blend
    "blueorange3": [BLUE, "#e9e2d0", ORANGE],
}
cmap = LinearSegmentedColormap.from_list(CMAP, CMAPS[CMAP]) if CMAP in CMAPS else plt.get_cmap(CMAP)

def subdir(reg):
    if reg == "n50_g3_dom500": return "cov10_n50_g3_s42_hotspots_dom500_p231_chr1"
    n, g = reg.split("_g"); return f"cov10_{n}_g{g}_s42_hotspots_p231_chr1"

tr = pd.read_csv(SIMS/subdir(REG)/"recomb_truth_raw.tsv.gz", sep="\t", usecols=["truth_af"])
es = pd.read_csv(RES/"kmate_global_filt2invu_raw"/REG/f"p231_filt2invu_raw_{REG}_cov10_s42.tsv",
                 sep="\t", usecols=["alt_freq", "n_called"])
assert len(tr) == len(es), "len mismatch"
truth = tr["truth_af"].values.astype(np.float32)
est   = es["alt_freq"].values.astype(np.float32)
miss  = ((F - es["n_called"].values) / F).astype(np.float32)

sel = (miss <= MISS_THR) & np.isfinite(truth) & np.isfinite(est)
x, y = truth[sel], est[sel]
keep = (x > 0) | (y > 0)
x, y = x[keep], y[keep]

d = y - x
mae = float(np.abs(d).mean()); rmse = float(math.sqrt(np.mean(d**2)))
ss = float(np.sum(d**2)); st = float(np.sum((x - x.mean())**2)); r2 = 1 - ss/st
n = len(x)
print(f"{REG} miss<= {MISS_THR}: n={n:,} MAE={mae:.4f} RMSE={rmse:.4f} R2={r2:.4f}")

nb = 100
H, xe, ye = np.histogram2d(x, y, bins=nb, range=[[0, 1], [0, 1]])
ix = np.clip(np.digitize(x, xe) - 1, 0, nb-1); iy = np.clip(np.digitize(y, ye) - 1, 0, nb-1)
dens = H[ix, iy].astype(np.float32)

SUB = 200_000
if n > SUB:
    rng = np.random.default_rng(42); pick = rng.choice(n, SUB, replace=False)
    x, y, dens = x[pick], y[pick], dens[pick]

plt.rcParams.update({"font.size": 15})
fig, ax = plt.subplots(figsize=(6.4, 6.4))
order = np.argsort(dens)
norm = LogNorm(vmin=max(dens.min(), 1), vmax=max(dens.max(), 1))
sc = ax.scatter(x[order], y[order], c=dens[order], s=11, cmap=cmap, norm=norm,
                alpha=0.95, edgecolors="none", rasterized=True, zorder=3)
ax.plot([0, 1], [0, 1], "--", lw=0.9, color=GREY, alpha=0.8, zorder=4)

# --- paper-panel aesthetic: no spines, grey ticks/labels, light background grid ---
ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
ax.set_aspect("equal", adjustable="box")
ax.set_xticks([0, .2, .4, .6, .8, 1.0]); ax.set_yticks([0, .2, .4, .6, .8, 1.0])
for sp in ("top", "right", "left", "bottom"):
    ax.spines[sp].set_visible(False)
ax.tick_params(length=0, colors=GREY)
ax.grid(True, color="#dddddd", linewidth=0.7, zorder=0); ax.set_axisbelow(True)
ax.set_xlabel("True allele frequency", fontsize=17); ax.xaxis.label.set_color(GREY)
ax.set_ylabel("Estimated allele frequency", fontsize=17); ax.yaxis.label.set_color(GREY)
ax.set_title(f"MAE = {mae:.4f}   R² = {r2:.3f}   n = {n:,}", fontsize=14, color=GREY)

cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
cb.set_label("Records (log density)", fontsize=12, color=GREY)
cb.outline.set_edgecolor(GREY)
cb.ax.yaxis.set_tick_params(color=GREY, labelcolor=GREY, labelsize=11)
fig.tight_layout()
tag = REG + ("" if CMAP == "viridis" else f"_{CMAP}")
misslab = f"miss{int(round(MISS_THR*100))}"
for ext in ("png", "pdf"):
    out = PLOTS/f"poster_accuracy_{tag}_{misslab}.{ext}"
    fig.savefig(out, dpi=300, bbox_inches="tight"); print("saved:", out)

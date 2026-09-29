#!/usr/bin/env python
"""Scatter of per-site SV purging vs per-site parallelism (the Q3 confound test), analogous to the
Section-1 climate scatter but with parallelism / selection-intensity on x. Points colored by bio1 so
you can see climate is NOT aligned with the parallelism axis. Env: basic. Reads site_parallelism.csv.
Saves to notebooks/plots/site_parallelism_vs_purging.png . No titles (repo convention)."""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy import stats
from adjustText import adjust_text
G = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r1_sv_negative_selection/results/sv_adaptive"
PLOTS = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/notebooks/plots"
df = pd.read_csv(f"{G}/site_parallelism.csv")

xs = [("par_pearson", "parallelism  (mean pairwise Pearson corr across plots)"),
      ("intensity", "selection intensity  (mean |founder logit-slope|)"),
      ("diversity", "evolved diversity  (eff. # surviving founders, 1/Σh²)\n(low = few ecotypes won = strong selection)")]
norm = mpl.colors.Normalize(vmin=df.bio1.min(), vmax=df.bio1.max()); cmap = mpl.cm.coolwarm
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
for ax, (col, lab) in zip(axes, xs):
    ax.set_axisbelow(True); ax.grid(color="0.9", lw=0.7)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.axhline(0, color="0.6", lw=0.8, ls=":")
    sc = ax.scatter(df[col], df.shift_sv, c=df.bio1, cmap=cmap, norm=norm, s=55, zorder=3)
    ax.margins(0.10)   # room so edge points' labels aren't clipped by the axes
    b, a0 = np.polyfit(df[col], df.shift_sv, 1); x = np.array([df[col].min(), df[col].max()])
    ax.plot(x, a0 + b * x, color="0.35", lw=1.4, ls="--")
    r = stats.spearmanr(df[col], df.shift_sv)
    ax.annotate(f"ρ = {r.statistic:+.2f}   p = {r.pvalue:.3f}", xy=(0.5, 0.98),
                xycoords="axes fraction", ha="center", va="top", fontsize=9, color="0.25",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", alpha=0.85))
    texts = [ax.text(rr[col], rr.shift_sv, str(int(rr.site)), fontsize=5.5, color="0.45")
             for _, rr in df.iterrows()]
    adjust_text(texts, ax=ax, expand=(1.4, 1.6), force_text=(0.4, 0.6),
                arrowprops=dict(arrowstyle="-", color="0.75", lw=0.4))
    ax.set_xlabel(lab, fontsize=8); ax.set_ylabel("SV purging   mean (s − SNP baseline)")
cb = fig.colorbar(sc, ax=axes, fraction=0.025, pad=0.01); cb.set_label("site bio1 (°C)", fontsize=8)
fig.savefig(f"{PLOTS}/plots/site_parallelism_vs_purging.png", dpi=130, bbox_inches="tight")
print(f"saved {PLOTS}/site_parallelism_vs_purging.png")

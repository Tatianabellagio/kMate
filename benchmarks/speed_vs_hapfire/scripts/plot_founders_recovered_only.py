#!/usr/bin/env python3
"""Founders-recovered-only figure (no speed row), all N including 231, plotted
as recovery FRACTION (n_found/N) so all N share one axis. Same visual language
as plot_hapfire_vs_kmate.py's founders_recovered figure but stripped to this one
metric. Reads results/hapfire_vs_kmate_table.tsv (from score_hapfire_vs_kmate.py).

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/speed_vs_hapfire/scripts/plot_founders_recovered_only.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/global/scratch/users/tbellg/kmate"
RES = f"{ROOT}/benchmarks/speed_vs_hapfire/results"
df = pd.read_csv(f"{RES}/hapfire_vs_kmate_table.tsv", sep="\t")
df["kmate_frac_found"] = df["kmate_n_found"] / df["N"]
df["hapfire_frac_found"] = df["hapfire_n_found"] / df["N"]

POOL_SIZES = [2, 5, 20, 50, 150, 231]
DEPTHS = [1, 10]
TOOL_COLOR = {"kmate": "#54a24b", "hapfire": "#e45756"}
TOOL_LABEL = {"kmate": "kMate", "hapfire": "hapFIRE"}

GREY = "#4d4d4d"
plt.rcParams.update({
    "text.color": GREY, "axes.labelcolor": GREY,
    "xtick.color": GREY, "ytick.color": GREY, "font.family": "sans-serif",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.spines.left": False, "axes.spines.bottom": False,
    "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.8,
    "axes.axisbelow": True, "xtick.bottom": False, "ytick.left": False,
})
JITTER_RNG = np.random.RandomState(0)


def box_panel(ax, sub, kmate_col, hapfire_col):
    x = np.arange(len(POOL_SIZES))
    w = 0.32
    for ti, (tool, col) in enumerate([("kmate", kmate_col), ("hapfire", hapfire_col)]):
        vals = [sub[sub.N == n][col].dropna().values for n in POOL_SIZES]
        pos = x + (ti - 0.5) * w
        ax.boxplot(vals, positions=pos, widths=w * 0.9, patch_artist=True,
                   showfliers=False, medianprops=dict(color="0.3", lw=1.3),
                   boxprops=dict(facecolor="none", edgecolor="0.3"),
                   whiskerprops=dict(color="0.3"), capprops=dict(color="0.3"))
        for p, v in zip(pos, vals):
            if len(v) == 0: continue
            jitter = JITTER_RNG.uniform(-w * 0.25, w * 0.25, size=len(v))
            ax.scatter(np.full(len(v), p) + jitter, v, color=TOOL_COLOR[tool],
                       s=22, alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([str(n) for n in POOL_SIZES])


def legend_handles():
    from matplotlib.lines import Line2D
    return [Line2D([0], [0], marker="o", color="none", markerfacecolor=TOOL_COLOR[t],
                  markeredgecolor="white", markersize=7, label=TOOL_LABEL[t])
            for t in ("kmate", "hapfire")]


fig, axes = plt.subplots(1, len(DEPTHS), figsize=(4.6 * len(DEPTHS), 4.0),
                         squeeze=False, sharey=True)
for ci, depth in enumerate(DEPTHS):
    sub = df[df.depth == depth]
    box_panel(axes[0][ci], sub, "kmate_frac_found", "hapfire_frac_found")
    axes[0][ci].text(0.5, 1.04, f"{depth}×", transform=axes[0][ci].transAxes,
                     fontsize=10, color="#999999", ha="center", va="bottom")
    axes[0][ci].set_xlabel("number of pooled founders")
    if ci == 0: axes[0][ci].set_ylabel("founders recovered / N pooled", fontsize=10)

fig.tight_layout()
fig.legend(handles=legend_handles(), loc="lower center", ncol=2, fontsize=9,
           bbox_to_anchor=(0.5, -0.11), frameon=False)
out_path = f"{RES}/hapfire_vs_kmate_founders_recovered_only.png"
fig.savefig(out_path, dpi=140, bbox_inches="tight")
plt.close(fig)
print("saved", out_path)

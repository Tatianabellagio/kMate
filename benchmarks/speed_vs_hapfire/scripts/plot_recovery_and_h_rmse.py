#!/usr/bin/env python3
"""h-accuracy variant: top panel = founders recovered (raw counts of N pooled),
bottom panel = RMSE/sd(truth) for founder proportions (h). Replaces the R2 top of
hapfire_vs_kmate_h_accuracy.png with founder recovery. 2 rows x 2 depths (1x,10x).
Reads results/hapfire_vs_kmate_table.tsv (from score_hapfire_vs_kmate.py).

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/speed_vs_hapfire/scripts/plot_recovery_and_h_rmse.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/benchmarks/speed_vs_hapfire/results"
df = pd.read_csv(f"{OUT}/hapfire_vs_kmate_table.tsv", sep="\t")

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


def box_panel(ax, sub, kmate_col, hapfire_col, log=False):
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
                       s=20, alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([str(n) for n in POOL_SIZES])
    if log:
        ax.set_yscale("log")


def legend_handles():
    from matplotlib.lines import Line2D
    return [Line2D([0], [0], marker="o", color="none", markerfacecolor=TOOL_COLOR[t],
                  markeredgecolor="white", markersize=7, label=TOOL_LABEL[t])
            for t in ("kmate", "hapfire")]


fig, axes = plt.subplots(2, len(DEPTHS), figsize=(4.6 * len(DEPTHS), 6.6),
                         squeeze=False, sharex="col", sharey="row")
for ci, depth in enumerate(DEPTHS):
    sub = df[df.depth == depth]
    # TOP: founders recovered (raw counts)
    box_panel(axes[0][ci], sub, "kmate_n_found", "hapfire_n_found")
    axes[0][ci].text(0.5, 1.04, f"{depth}×", transform=axes[0][ci].transAxes,
                     fontsize=10, color="#999999", ha="center", va="bottom")
    if ci == 0: axes[0][ci].set_ylabel("founders recovered (of N pooled)", fontsize=10)
    # BOTTOM: raw h RMSE (log) -- defined at every N incl. 231 (unlike the
    # sd-normalized version, which is NaN at N=231 where truth variance is 0)
    box_panel(axes[1][ci], sub, "kmate_h_RMSE", "hapfire_h_RMSE", log=True)
    axes[1][ci].set_xlabel("number of pooled founders")
    if ci == 0: axes[1][ci].set_ylabel("RMSE (h, raw, log)\n[defined at every N]", fontsize=10)

fig.tight_layout()
fig.legend(handles=legend_handles(), loc="lower center", ncol=2, fontsize=9,
           bbox_to_anchor=(0.5, -0.03), frameon=False)
out_path = f"{OUT}/hapfire_vs_kmate_recovery_and_h_rmse.png"
fig.savefig(out_path, dpi=140, bbox_inches="tight")
plt.close(fig)
print("saved", out_path)

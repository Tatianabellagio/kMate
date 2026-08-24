#!/usr/bin/env python3
"""h-accuracy plot for the hapFIRE-vs-kMate poolsize x depth comparison --
each tool on its OWN native panel (kMate: arch3, hapFIRE: greneNet), matched
pools. Founder-recovery + speed figures now come from
plot_recovery_and_h_rmse.py / plot_speed_vs_coverage.py.
Reads benchmarks/speed_vs_hapfire/results/hapfire_vs_kmate_table.tsv (see
score_hapfire_vs_kmate.py). Same visual language as
benchmarks/poolsize_depth/scripts/score_and_plot.py (transparent boxes, jittered
colored points, shared y per row) but hue = tool (kMate vs hapFIRE) instead of
depth, faceted by depth instead of by panel.

Depth columns: 1x / 5x / 10x -- the 5x column is the intermediate depth added to
show the N=150 low-coverage weakness closing; both tools run there (full fair
grid), same as 1x/10x.

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/speed_vs_hapfire/scripts/plot_hapfire_vs_kmate.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/benchmarks/speed_vs_hapfire/results"
os.makedirs(f"{OUT}/plots", exist_ok=True)
df = pd.read_csv(f"{OUT}/hapfire_vs_kmate_table.tsv", sep="\t")

POOL_SIZES = [2, 5, 20, 50, 150, 231]
DEPTHS = [1, 5, 10]
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


def make_figure(rows, out_path):
    n_rows = len(rows)
    fig, axes = plt.subplots(n_rows, len(DEPTHS), figsize=(4.2 * len(DEPTHS), 3.3 * n_rows),
                             squeeze=False, sharex="col", sharey="row")
    for ci, depth in enumerate(DEPTHS):
        sub = df[df.depth == depth]
        for ri, (kmate_col, hapfire_col, label, log) in enumerate(rows):
            box_panel(axes[ri][ci], sub, kmate_col, hapfire_col, log=log)
            if ri == 0:
                axes[ri][ci].text(0.5, 1.04, f"{depth}×", transform=axes[ri][ci].transAxes,
                                 fontsize=10, fontweight="normal", color="#999999",
                                 ha="center", va="bottom")
            if ci == 0: axes[ri][ci].set_ylabel(label, fontsize=10)
        axes[n_rows - 1][ci].set_xlabel("number of pooled founders")
    fig.tight_layout()
    fig.legend(handles=legend_handles(), loc="lower center", ncol=2, fontsize=9,
               bbox_to_anchor=(0.5, -0.03), frameon=False)
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved", out_path)


make_figure([
    ("kmate_h_R2", "hapfire_h_R2", "R² (h)", False),
    ("kmate_h_RMSE_norm", "hapfire_h_RMSE_norm", "RMSE / sd(truth) (h)", False),
    ("kmate_h_RMSE", "hapfire_h_RMSE", "RMSE (h)", False),
], f"{OUT}/plots/hapfire_vs_kmate_h_accuracy.png")

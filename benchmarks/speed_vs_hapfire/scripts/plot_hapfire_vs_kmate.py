#!/usr/bin/env python3
"""Plots for the hapFIRE-vs-kMate poolsize x depth comparison (h-accuracy,
founder recovery, speed/compute) -- each tool on its OWN native panel (kMate:
arch3, hapFIRE: greneNet), matched pools. See ../PANEL_MISMATCH_BUG.md for why
there is no AF-accuracy figure: the two tools' native variant catalogs are
different calling lineages, so a single shared "truth AF" isn't well-defined.
Reads benchmarks/speed_vs_hapfire/results/hapfire_vs_kmate_table.tsv (see
score_hapfire_vs_kmate.py). Same visual language as
benchmarks/poolsize_depth/scripts/score_and_plot.py (transparent boxes, jittered
colored points, shared y per row) but hue = tool (kMate vs hapFIRE) instead of
depth, faceted by depth instead of by panel.

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/speed_vs_hapfire/scripts/plot_hapfire_vs_kmate.py
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


def make_figure(row1_kmate, row1_hapfire, row1_label, row2_kmate, row2_hapfire, row2_label,
               out_path, row1_log=False, row2_log=False):
    fig, axes = plt.subplots(2, len(DEPTHS), figsize=(4.6 * len(DEPTHS), 6.6),
                             squeeze=False, sharex="col", sharey="row")
    for ci, depth in enumerate(DEPTHS):
        sub = df[df.depth == depth]
        box_panel(axes[0][ci], sub, row1_kmate, row1_hapfire, log=row1_log)
        axes[0][ci].text(0.5, 1.04, f"{depth}×", transform=axes[0][ci].transAxes,
                         fontsize=10, fontweight="normal", color="#999999",
                         ha="center", va="bottom")
        if ci == 0: axes[0][ci].set_ylabel(row1_label, fontsize=10)
        box_panel(axes[1][ci], sub, row2_kmate, row2_hapfire, log=row2_log)
        axes[1][ci].set_xlabel("number of pooled founders")
        if ci == 0: axes[1][ci].set_ylabel(row2_label, fontsize=10)
    fig.tight_layout()
    fig.legend(handles=legend_handles(), loc="lower center", ncol=2, fontsize=9,
               bbox_to_anchor=(0.5, -0.03), frameon=False)
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved", out_path)


make_figure("kmate_h_R2", "hapfire_h_R2", "R² (h)",
            "kmate_h_RMSE_norm", "hapfire_h_RMSE_norm", "RMSE / sd(truth) (h)",
            f"{OUT}/hapfire_vs_kmate_h_accuracy.png")

make_figure("kmate_n_found", "hapfire_n_found", "founders recovered (of N pooled)",
            "kmate_elapsed_s", "hapfire_total_elapsed_s", "wall-clock (s, log)",
            f"{OUT}/hapfire_vs_kmate_founders_recovered.png", row2_log=True)

make_figure("kmate_elapsed_s", "hapfire_total_elapsed_s", "wall-clock (s, log)",
            "kmate_cpu_s", "hapfire_total_cpu_s", "total CPU-seconds (log)",
            f"{OUT}/hapfire_vs_kmate_speed.png", row1_log=True, row2_log=True)

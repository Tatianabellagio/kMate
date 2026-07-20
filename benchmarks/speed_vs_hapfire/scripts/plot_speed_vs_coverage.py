#!/usr/bin/env python3
"""Speed/compute vs sequencing depth, at FIXED pool size N=50 (isolates the
coverage-scaling question from the pool-size question -- see
hapfire_vs_kmate_speed.png for the N-swept view). Coverage in {1,10,30,50}x.
Each tool on its own native panel (kMate: arch3, hapFIRE: greneNet); see
../PANEL_MISMATCH_BUG.md. Reads
benchmarks/speed_vs_hapfire/results/hapfire_vs_kmate_table.tsv (see
score_hapfire_vs_kmate.py -- DEPTHS extended to include 30/50, N=50 only).

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/speed_vs_hapfire/scripts/plot_speed_vs_coverage.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, ScalarFormatter, FuncFormatter

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/benchmarks/speed_vs_hapfire/results"
df = pd.read_csv(f"{OUT}/hapfire_vs_kmate_table.tsv", sep="\t")
df = df[df.N == 50]

DEPTHS = [1, 10, 30, 50]
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


def box_panel(ax, kmate_col, hapfire_col):
    x = np.arange(len(DEPTHS))
    w = 0.32
    for ti, (tool, col) in enumerate([("kmate", kmate_col), ("hapfire", hapfire_col)]):
        vals = [df[df.depth == d][col].dropna().values for d in DEPTHS]
        pos = x + (ti - 0.5) * w
        ax.boxplot(vals, positions=pos, widths=w * 0.9, patch_artist=True,
                   showfliers=False, medianprops=dict(color="0.3", lw=1.3),
                   boxprops=dict(facecolor="none", edgecolor="0.3"),
                   whiskerprops=dict(color="0.3"), capprops=dict(color="0.3"))
        for p, v in zip(pos, vals):
            if len(v) == 0: continue
            jitter = JITTER_RNG.uniform(-w * 0.25, w * 0.25, size=len(v))
            ax.scatter(np.full(len(v), p) + jitter, v, color=TOOL_COLOR[tool],
                       s=24, alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}×" for d in DEPTHS])
    ax.set_xlabel("sequencing depth (N=50 founders)")
    # log scale, but only 2x/5x per decade as labeled reference lines (not all
    # of 2-9, which gave unlabeled, increasingly-cramped gridlines with no way
    # to read off an approximate value away from a bare power of 10).
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(2, 5)))
    ax.yaxis.set_minor_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.tick_params(axis="y", which="minor", labelsize=8, labelcolor="#999999")
    ax.grid(True, which="minor", axis="y", color="#ececec", linewidth=0.6)
    ax.grid(True, which="major", axis="y", color="#dddddd", linewidth=0.8)


fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.6))
box_panel(axes[0], "kmate_elapsed_s", "hapfire_total_elapsed_s")
axes[0].set_ylabel("wall-clock (s, log)", fontsize=10)
box_panel(axes[1], "kmate_cpu_s", "hapfire_total_cpu_s")
axes[1].set_ylabel("total CPU-seconds (log)", fontsize=10)

from matplotlib.lines import Line2D
handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=TOOL_COLOR[t],
                  markeredgecolor="white", markersize=7, label=TOOL_LABEL[t])
          for t in ("kmate", "hapfire")]
fig.tight_layout()
fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=9,
           bbox_to_anchor=(0.5, -0.06), frameon=False)
out_path = f"{OUT}/hapfire_vs_kmate_speed_vs_coverage.png"
fig.savefig(out_path, dpi=140, bbox_inches="tight")
plt.close(fig)
print("saved", out_path)

#!/usr/bin/env python3
"""Figure for the native-panel ecotype-resolution comparison (kMate vs hapFIRE).

Reads ecotype_nativepanel_table.tsv (score_ecotype_nativepanel.py) and renders the
resolution metrics that speed_vs_hapfire's h-accuracy figure does NOT cover.
Matches the house style of hapfire_vs_kmate_af_accuracy.png: box + strip over
seeds, one column per depth, tools as colour.

Project convention: no chart titles, no subplot titles (CLAUDE.md). Panel identity
goes in an in-panel corner annotation; depth goes in the column annotation.

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/ecotype_count/scripts/plot_ecotype_nativepanel.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = "/global/scratch/users/tbellg/kmate"
RES = f"{ROOT}/benchmarks/ecotype_count/results"
PLOTS = f"{RES}/plots"

# Same palette/role mapping as the sibling native-panel figures.
TOOLS = [("kMate", "#4c9a52"), ("hapFIRE", "#d1495b")]
METRICS = [("recall", "recall (true founders detected)"),
           ("mass_on_true", "mass on true founders"),
           ("top_n_jaccard", "Jaccard(top-N, truth)")]
DEPTHS = [1, 10]          # the depths every N is run at; N=50 also has 30/50


def main():
    df = pd.read_csv(f"{RES}/ecotype_nativepanel_table.tsv", sep="\t")
    df = df[df.depth.isin(DEPTHS)]
    Ns = sorted(df.N.unique())
    x = np.arange(len(Ns))

    fig, axes = plt.subplots(len(METRICS), len(DEPTHS),
                             figsize=(5.4 * len(DEPTHS), 3.0 * len(METRICS)),
                             squeeze=False, sharex=True, sharey="row")

    for ri, (metric, ylab) in enumerate(METRICS):
        for ci, depth in enumerate(DEPTHS):
            ax = axes[ri][ci]
            sub = df[df.depth == depth]
            for ti, (tool, col) in enumerate(TOOLS):
                off = (ti - (len(TOOLS) - 1) / 2) * 0.3
                vals = [sub[(sub.tool == tool) & (sub.N == n)][metric].dropna().values
                        for n in Ns]
                ax.boxplot(vals, positions=x + off, widths=0.24, showfliers=False,
                           medianprops=dict(color="0.25"),
                           boxprops=dict(color="0.45"),
                           whiskerprops=dict(color="0.45"),
                           capprops=dict(color="0.45"))
                for xi, v in zip(x + off, vals):
                    if len(v):
                        ax.scatter(np.full(len(v), xi) + np.linspace(-.05, .05, len(v)),
                                   v, s=14, color=col, zorder=3, alpha=0.85,
                                   edgecolors="none")
            ax.grid(axis="y", ls=":", alpha=0.4)
            if ci == 0:
                ax.set_ylabel(ylab)
            if ri == 0:
                ax.text(0.5, 1.02, f"{depth}×", transform=ax.transAxes,
                        ha="center", fontsize=11, color="0.4")
    for ax in axes[-1]:
        ax.set_xticks(x)
        ax.set_xticklabels(Ns)
        ax.set_xlabel("number of pooled founders")

    handles = [Patch(fc=c, label=t) for t, c in TOOLS]
    fig.legend(handles=handles, loc="lower center", ncol=len(TOOLS), fontsize=10,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.03, 1, 1])

    os.makedirs(PLOTS, exist_ok=True)
    out = f"{PLOTS}/ecotype_nativepanel_resolution.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"saved -> {out}  ({len(df)} rows, N={Ns}, depths={DEPTHS})")


if __name__ == "__main__":
    main()

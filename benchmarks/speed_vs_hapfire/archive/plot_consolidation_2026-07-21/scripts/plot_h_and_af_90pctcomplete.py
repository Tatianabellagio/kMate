#!/usr/bin/env python3
"""Combined kMate-vs-hapFIRE figure: founder-proportion accuracy (h) + per-SNP
allele-frequency accuracy (AF), with kMate's AF restricted to >=90%-called sites
(the fair basis vs hapFIRE's zero-missingness panel). 5 rows x 2 depths:
  1. R2 (h)
  2. RMSE/sd(truth) (h)         [NaN at N=231: truth has zero variance there]
  3. RMSE (h, raw, log)         [defined at every N, incl. 231]
  4. R2 (AF, >=90% called)
  5. RMSE (AF, >=90% called)

Reads (all under results/): hapfire_vs_kmate_table.tsv (h), hapfire_vs_kmate_af_table.tsv
(hapFIRE AF), kmate_af_missingness_percondition.tsv (kMate >=90%-called AF, from
score_af_missingness_thresholds.py -- run that first). Same visual language as
plot_hapfire_vs_kmate.py.

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/speed_vs_hapfire/scripts/plot_h_and_af_90pctcomplete.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/global/scratch/users/tbellg/kmate"
RES = f"{ROOT}/benchmarks/speed_vs_hapfire/results"

h_df = pd.read_csv(f"{RES}/hapfire_vs_kmate_table.tsv", sep="\t")
af_hf = pd.read_csv(f"{RES}/hapfire_vs_kmate_af_table.tsv", sep="\t")[
    ["N", "depth", "seed", "hapfire_af_R2", "hapfire_af_RMSE"]]
p90 = pd.read_csv(f"{RES}/kmate_af_missingness_percondition.tsv", sep="\t")[
    ["N", "depth", "seed", "p90_R2", "p90_RMSE"]]
af_df = af_hf.merge(p90, on=["N", "depth", "seed"], how="inner")
print(f"AF merged rows: {len(af_df)} (expect 60)")

POOL_SIZES = [2, 5, 20, 50, 150, 231]
DEPTHS = [1, 10]
TOOL_COLOR = {"kmate": "#54a24b", "hapfire": "#e45756"}
TOOL_LABEL = {"kmate": "kMate (>=90% called)", "hapfire": "hapFIRE"}

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


fig, axes = plt.subplots(5, len(DEPTHS), figsize=(4.6 * len(DEPTHS), 15.6),
                         squeeze=False, sharex="col")
for ci, depth in enumerate(DEPTHS):
    hsub = h_df[h_df.depth == depth]
    box_panel(axes[0][ci], hsub, "kmate_h_R2", "hapfire_h_R2")
    axes[0][ci].text(0.5, 1.04, f"{depth}×", transform=axes[0][ci].transAxes,
                     fontsize=10, color="#999999", ha="center", va="bottom")
    if ci == 0: axes[0][ci].set_ylabel("R² (h)", fontsize=10)

    box_panel(axes[1][ci], hsub, "kmate_h_RMSE_norm", "hapfire_h_RMSE_norm")
    if ci == 0: axes[1][ci].set_ylabel("RMSE / sd(truth) (h)\n[NaN at N=231: truth has\nzero variance there]", fontsize=9)

    box_panel(axes[2][ci], hsub, "kmate_h_RMSE", "hapfire_h_RMSE", log=True)
    if ci == 0: axes[2][ci].set_ylabel("RMSE (h, raw, log)\n[defined at every N]", fontsize=9)

    asub = af_df[af_df.depth == depth]
    box_panel(axes[3][ci], asub, "p90_R2", "hapfire_af_R2")
    if ci == 0: axes[3][ci].set_ylabel("R² (AF, ≥90% called)", fontsize=10)

    box_panel(axes[4][ci], asub, "p90_RMSE", "hapfire_af_RMSE")
    axes[4][ci].set_xlabel("number of pooled founders")
    if ci == 0: axes[4][ci].set_ylabel("RMSE (AF, ≥90% called)", fontsize=10)

for ax_row in axes:
    for ax in ax_row:
        ax.margins(y=0.08)

fig.tight_layout()
fig.legend(handles=legend_handles(), loc="lower center", ncol=2, fontsize=9,
           bbox_to_anchor=(0.5, -0.015), frameon=False)
out_path = f"{RES}/hapfire_vs_kmate_h_and_af_90pctcomplete.png"
fig.savefig(out_path, dpi=140, bbox_inches="tight")
plt.close(fig)
print("saved", out_path)

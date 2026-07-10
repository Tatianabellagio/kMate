#!/usr/bin/env python3
"""h-imbalance grid: does the private-k-mer filter still matter under the current
production estimator (--unit chrom, per-founder M-step normalization)?

Two views of the same data:
  1. scatter grid  - per-founder h estimate vs truth, long-read/short-read colored,
                     2 rows (raw / filt2inv) x 3 columns (regimes)
  2. error boxplot - distribution of per-founder signed error (est - truth), split
                     by founder class (long-read/short-read), one panel per regime,
                     boxes grouped by arm

Data: benchmarks/p231/results/kmate_chrom_himbalance_{raw,filt2inv}/<regime>/*.h_per_chrom.npz
Produced by: benchmarks/p231/scripts/09_run_h_imbalance_chrom_p231.sh

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/h_imbalance/scripts/plot_h_imbalance.py
"""
import os, json, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT  = "/global/scratch/users/tbellg/kmate"
P231  = f"{ROOT}/benchmarks/p231/results"
SIMS  = f"{ROOT}/benchmarks/p231/sims"
OUT   = f"{ROOT}/benchmarks/h_imbalance/results"
os.makedirs(OUT, exist_ok=True)

_split   = json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))
LONGREAD = set(map(str, _split["cactus"]))          # 80 long-read founders

C_LONG, C_SHORT, C_TRUTH = "#006d2c", "#a1d76a", "0.45"
LAB_LONG, LAB_SHORT = "long-read founders", "short-read founders"

APPROACHES = [("raw\n(no filter)", "raw"), ("private k-mer filter\n(production)", "filt2inv")]
COLUMNS = [
    ("n231_g0\n(uniform, all 231)",          "n231_g0",       False),
    ("n50 long-read-heavy\n(40 LR + 10 SR)", "n50_cactheavy", True),
    ("n50 short-read-heavy\n(45 SR + 5 LR)", "n50_pgheavy",   True),
]

GREY = "#4d4d4d"
plt.rcParams.update({
    "text.color": GREY, "axes.labelcolor": GREY, "axes.titlecolor": GREY,
    "xtick.color": GREY, "ytick.color": GREY, "font.family": "sans-serif",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.spines.left": False, "axes.spines.bottom": False,
    "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.8,
    "axes.axisbelow": True, "xtick.bottom": False, "ytick.left": False,
})


def sim_dir(regime):
    return f"{SIMS}/cov10_{regime}_s42_hotspots_p231_chr1"


def load_h(arm_tag, regime):
    f = f"{P231}/kmate_chrom_himbalance_{arm_tag}/{regime}/p231_chrom_himbalance_{arm_tag}_{regime}_cov10_s42.h_per_chrom.npz"
    if not os.path.exists(f):
        return None, None
    d = np.load(f, allow_pickle=True)
    return d["founders"].astype(str), d["Chr1"].astype(float)


def load_truth(founders, regime):
    th = {}
    for r in csv.DictReader(open(f"{sim_dir(regime)}/pool_weights.tsv"), delimiter="\t"):
        th[str(r["founder"])] = float(r["weight"])
    return np.array([th.get(f, 0.0) for f in founders])


def longmask(founders):
    return np.array([f in LONGREAD for f in founders])


def plot_scatter_grid():
    fig, axes = plt.subplots(len(APPROACHES), len(COLUMNS),
                             figsize=(4.8 * len(COLUMNS), 3.7 * len(APPROACHES)),
                             squeeze=False, sharey=True)
    for ci, (clab, regime, by_truth) in enumerate(COLUMNS):
        for ri, (rlab, arm) in enumerate(APPROACHES):
            ax = axes[ri, ci]
            fo, h = load_h(arm, regime)
            if fo is None:
                ax.text(0.5, 0.5, "pending", ha="center", va="center",
                        transform=ax.transAxes, color="0.5")
                ax.set_xticks([]); ax.set_yticks([])
                if ci == 0: ax.set_ylabel(rlab + "\n\n$h$", fontsize=9, color=GREY)
                if ri == 0: ax.annotate(clab, xy=(0.5, 1.16), xycoords="axes fraction",
                                        ha="center", va="bottom", fontsize=11, color=GREY)
                continue
            truth = load_truth(fo, regime)
            lm = longmask(fo)
            order = np.argsort(truth, kind="stable") if by_truth else np.arange(len(fo))
            x = np.arange(len(fo))
            t_o, h_o, lm_o = truth[order], h[order], lm[order]
            ax.scatter(x, t_o, s=6, color=C_TRUTH, zorder=1)
            ax.scatter(x[~lm_o], h_o[~lm_o], marker="x", s=15, lw=.8, color=C_SHORT, zorder=2)
            ax.scatter(x[lm_o],  h_o[lm_o],  marker="x", s=15, lw=.9, color=C_LONG,  zorder=3)
            rmse = float(np.sqrt(np.mean((h - truth) ** 2)))
            ax.set_title(f"RMSE {rmse:.4f}", fontsize=10, color=GREY)
            ax.tick_params(colors=GREY)
            if ci == 0: ax.set_ylabel(rlab + "\n\n$h$", fontsize=9, color=GREY)
            if ri == len(APPROACHES) - 1:
                ax.set_xlabel("founder (ordered by true freq)" if by_truth else "founder index")
            if ri == 0: ax.annotate(clab, xy=(0.5, 1.16), xycoords="axes fraction",
                                    ha="center", va="bottom", fontsize=11, color=GREY)
    handles = [Line2D([], [], marker='x', ls='', color=C_LONG,  label=LAB_LONG),
               Line2D([], [], marker='x', ls='', color=C_SHORT, label=LAB_SHORT),
               Line2D([], [], marker='o', ls='', color=C_TRUTH, label="truth")]
    fig.tight_layout()
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9,
               bbox_to_anchor=(0.5, -0.02), frameon=False)
    out = f"{OUT}/h_imbalance_grid_p231_chrom.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("saved", out)


def plot_error_boxplot():
    """Distribution of per-founder signed error (est - truth), split by founder
    class, boxes grouped by arm -- one panel per regime (hapFIRE-paper style
    box-and-whisker, distribution over founders rather than seeds)."""
    fig, axes = plt.subplots(1, len(COLUMNS), figsize=(4.4 * len(COLUMNS), 4.2), sharey=False)
    if len(COLUMNS) == 1:
        axes = [axes]
    for ci, (clab, regime, _) in enumerate(COLUMNS):
        ax = axes[ci]
        box_data, positions, colors = [], [], []
        pos = 0
        for arm_i, (rlab, arm) in enumerate(APPROACHES):
            fo, h = load_h(arm, regime)
            if fo is None:
                pos += 3
                continue
            truth = load_truth(fo, regime)
            err = h - truth
            lm = longmask(fo)
            for cls_mask, color in ((lm, C_LONG), (~lm, C_SHORT)):
                box_data.append(err[cls_mask])
                positions.append(pos)
                colors.append(color)
                pos += 1
            pos += 1  # gap between arms
        bp = ax.boxplot(box_data, positions=positions, widths=0.8, patch_artist=True,
                         showfliers=True, flierprops=dict(marker='.', ms=3, alpha=0.5,
                                                           markeredgecolor="0.5"),
                         medianprops=dict(color="0.2", lw=1.3))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.75); patch.set_edgecolor("0.3")
        # x-tick per arm, centered between its two boxes
        arm_centers = [1.0 + 3 * i for i in range(len(APPROACHES))]
        ax.set_xticks(arm_centers)
        ax.set_xticklabels([rlab.replace("\n", " ") for rlab, _ in APPROACHES], fontsize=8.5)
        ax.axhline(0, color="0.6", lw=0.8, ls=":")
        ax.set_title(clab.replace("\n", " "), fontsize=10, color=GREY)
        if ci == 0: ax.set_ylabel("h error (est − truth)", fontsize=9, color=GREY)
        ax.tick_params(colors=GREY)
    handles = [Patch(fc=C_LONG, alpha=0.75, label=LAB_LONG),
               Patch(fc=C_SHORT, alpha=0.75, label=LAB_SHORT)]
    fig.tight_layout()
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=9,
               bbox_to_anchor=(0.5, -0.06), frameon=False)
    out = f"{OUT}/h_imbalance_error_boxplot_p231_chrom.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("saved", out)


if __name__ == "__main__":
    plot_scatter_grid()
    plot_error_boxplot()

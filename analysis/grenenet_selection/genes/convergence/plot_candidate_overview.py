#!/usr/bin/env python
"""Overview figures for the candidate shortlist: what was found, how, and in what theme.

Three panels, each answering one question, no chart or subplot titles (house rule --
context comes from the markdown around them and from in-panel corner annotations):

  A  the shortlist itself: one row per gene, ordered by evidence then mechanism.
     x = evidence strength (GEA climate clusters, and gardens for the GWAS side),
     marker = which scan found it, colour = theme, size = variant size.
  B  theme composition of the full candidate pool vs the shortlist -- does gating on
     evidence enrich for the biology of interest, or just shrink the list?
  C  where the shortlist sits on the genome, so the LD pileups are visible rather
     than implied. Loci contributing many genes are the ones to distrust.

env: kmate (renders headless).  Figures -> results/plots/
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                              # noqa: E402
from matplotlib.lines import Line2D                          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = f"{HERE}/results"
PLOTS = f"{OUT}/plots"

THEME_COL = {"stress": "#b2182b", "climate_temp": "#ef8a62",
             "circadian_light": "#5aae61", "flowering": "#762a83",
             "multiple": "#1a1a1a", "none": "#bdbdbd"}
SCAN_MARK = {"GEA+GWAS": ("*", 320), "GEA": ("o", 70), "GWAS": ("s", 70)}
CHROM_LEN = {"Chr1": 30427671, "Chr2": 19698289, "Chr3": 23459830,
             "Chr4": 18585056, "Chr5": 26975502}


def theme_of(s: str) -> str:
    if not isinstance(s, str) or not s:
        return "none"
    parts = [p for p in s.split(",") if p]
    if len(parts) > 1:
        return "multiple"
    return parts[0]


def panel_shortlist(ax, S):
    S = S.iloc[::-1].reset_index(drop=True)          # best at the top
    y = np.arange(len(S))
    # evidence strength: climate clusters for the GEA, gardens for the GWAS
    ax.barh(y, S.gea_n_clusters, color="#c6dbef", height=0.62, zorder=1,
            label="GEA climate clusters")
    ax.barh(y, -S.gwas_n_gardens, color="#fdd0a2", height=0.62, zorder=1,
            label="GWAS gardens")
    for r, yy in zip(S.itertuples(), y):
        mk, sz = SCAN_MARK[r.found_by]
        col = THEME_COL[theme_of(r.themes)]
        x = r.gea_n_clusters if r.gea_n_clusters >= r.gwas_n_gardens \
            else -r.gwas_n_gardens
        ax.scatter([x], [yy], marker=mk, s=sz, color=col, edgecolor="black",
                   linewidth=0.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([(r.symbol if isinstance(r.symbol, str) and r.symbol
                         else r.gene) for r in S.itertuples()], fontsize=7)
    ax.axvline(0, color="black", lw=0.8, zorder=2)
    ax.set_xlabel("← gardens (per-garden GWAS)      climate clusters (GEA) →")
    ax.tick_params(axis="x", labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.annotate("shortlist: ≥1 gated line of evidence", xy=(0.98, 0.01),
                xycoords="axes fraction", ha="right", va="bottom", fontsize=8,
                style="italic", color="#555555")


def panel_themes(ax, M, S):
    order = ["stress", "climate_temp", "circadian_light", "flowering"]
    pool = [M.themes.fillna("").str.contains(t).mean() * 100 for t in order]
    short = [S.themes.fillna("").str.contains(t).mean() * 100 for t in order]
    x = np.arange(len(order))
    ax.bar(x - 0.19, pool, 0.38, color="#bdbdbd", label=f"all candidates (n={len(M)})")
    ax.bar(x + 0.19, short, 0.38, color="#4292c6", label=f"shortlist (n={len(S)})")
    for i, (p, s) in enumerate(zip(pool, short)):
        ax.text(i - 0.19, p + 0.6, f"{p:.0f}%", ha="center", fontsize=7)
        ax.text(i + 0.19, s + 0.6, f"{s:.0f}%", ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(["stress", "temperature", "circadian/\nlight", "flowering"],
                       fontsize=8)
    ax.set_ylabel("% of genes carrying the theme")
    ax.legend(fontsize=7, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)


def panel_genome(ax, M, S):
    chroms = list(CHROM_LEN)
    for i, ch in enumerate(chroms):
        ax.plot([0, CHROM_LEN[ch] / 1e6], [i, i], color="#d9d9d9", lw=5,
                solid_capstyle="round", zorder=1)
    # every candidate faint, shortlist emphasised, marker size = genes in the locus
    for r in M.itertuples():
        if r.chrom not in CHROM_LEN or pd.isna(r.gene_start):
            continue
        ax.scatter([r.gene_start / 1e6], [chroms.index(r.chrom)], s=4,
                   color="#9ecae1", alpha=0.35, zorder=2, linewidth=0)
    for r in S.itertuples():
        if r.chrom not in CHROM_LEN or pd.isna(r.gene_start):
            continue
        ax.scatter([r.gene_start / 1e6], [chroms.index(r.chrom)],
                   s=12 + 3.0 * float(r.locus_n_genes),
                   color=THEME_COL[theme_of(r.themes)], edgecolor="black",
                   linewidth=0.4, alpha=0.9, zorder=3)
    ax.set_yticks(range(len(chroms)))
    ax.set_yticklabels(chroms, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("position (Mb)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.annotate("marker size = candidate genes in the locus (LD pileup)",
                xy=(0.99, 0.04), xycoords="axes fraction", ha="right",
                va="bottom", fontsize=8, style="italic", color="#555555")


def main():
    os.makedirs(PLOTS, exist_ok=True)
    M = pd.read_csv(f"{OUT}/master_candidate_genes.csv")
    S = M[M.n_independent >= 1].sort_values(
        ["n_independent", "gea_n_clusters", "gwas_n_gardens"], ascending=False)

    fig = plt.figure(figsize=(15, 13))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[3, 1.15],
                          hspace=0.16, wspace=0.28)
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])

    panel_shortlist(ax1, S)
    panel_themes(ax2, M, S)
    panel_genome(ax3, M, S)

    # legend goes INSIDE panel A's empty left half -- ax2 already carries the bar
    # legend, and a second ax2.legend() would silently replace it.
    handles = [Line2D([], [], marker=m, color="w", markerfacecolor="#777777",
                      markeredgecolor="black", markersize=min(np.sqrt(s) * 0.8, 13),
                      label=f"found by {k}", linestyle="")
               for k, (m, s) in SCAN_MARK.items()]
    handles += [Line2D([], [], marker="o", color="w", markerfacecolor=c,
                       markeredgecolor="black", markersize=7,
                       label=k.replace("climate_temp", "temperature")
                              .replace("circadian_light", "circadian/light")
                              .replace("multiple", "multiple themes")
                              .replace("none", "no theme"), linestyle="")
                for k, c in THEME_COL.items()]
    ax1.legend(handles=handles, fontsize=7.5, frameon=True, framealpha=0.95,
               edgecolor="#cccccc", ncol=2, loc="lower left",
               bbox_to_anchor=(0.005, 0.045))

    out = f"{PLOTS}/candidate_overview"
    fig.savefig(out + ".png", dpi=170, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print(f"wrote {out}.png / .pdf   ({len(S)} shortlist of {len(M)} candidates)")


if __name__ == "__main__":
    main()

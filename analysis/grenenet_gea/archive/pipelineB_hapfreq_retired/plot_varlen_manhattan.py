#!/usr/bin/env python
"""Two-panel Manhattan for the 31-site variable-length (varlen) hap-block GEA.

Top  = raw per-block test (block_gea.csv)      -> inflated (lambda ~ 8.6)
Bottom = permutation-calibrated block-WZA (block_wza.csv) -> lambda ~ 1.0

Mirrors the 19-site pipelineB/manhattan_raw_vs_wza.png so the two N's are visually
comparable. Run in kmate env; instant (16k blocks). Writes manhattan_raw_vs_wza.png
into the chosen PB_DIR.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

H = "results/grenenet_gea/hapfreq"
B = os.environ.get("PB_DIR", f"{H}/pipelineB_varlen")
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
COLORS = ["#3b4cc0", "#7aa0c4"]  # alternating chrom shades


def lam(p):
    p = np.clip(p, 1e-12, 1)
    return np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)


def cumpos(df):
    """Genome-cumulative x + chrom tick centers."""
    df = df.copy()
    df["mid"] = (df.unit_start + df.unit_end) / 2
    off, centers, ticks = 0.0, [], []
    df["x"] = 0.0
    for ch in CHROMS:
        m = df.chrom == ch
        if not m.any():
            continue
        cmax = df.loc[m, "mid"].max()
        df.loc[m, "x"] = df.loc[m, "mid"] + off
        centers.append(off + cmax / 2)
        off += cmax * 1.02
        ticks.append(ch)
    return df, centers, ticks


def panel(ax, df, title, draw_thresholds=True):
    df, centers, ticks = cumpos(df)
    df["nlp"] = -np.log10(np.clip(df.block_p, 1e-12, 1))
    for i, ch in enumerate(CHROMS):
        m = df.chrom == ch
        ax.scatter(df.loc[m, "x"], df.loc[m, "nlp"], s=7,
                   c=COLORS[i % 2], alpha=0.55, edgecolors="none", rasterized=True)
    m = len(df)
    # permutation floor: smallest achievable p = 1/(N_perm+1); here data min
    pfloor = df.block_p[df.block_p > 0].min()
    ax.axhline(-np.log10(pfloor), color="grey", lw=0.8, ls=":",
               label=f"perm floor (p={pfloor:.1e})")
    if draw_thresholds:
        bonf = 0.05 / m
        ax.axhline(-np.log10(bonf), color="firebrick", lw=0.9, ls="--",
                   label=f"Bonferroni 0.05/{m}")
        # BH p-threshold at q<0.05 (largest p with q<0.05), if any
        ps = np.sort(df.block_p.values)
        bh = ps[ps <= 0.05 * (np.arange(1, m + 1)) / m]
        if bh.size:
            ax.axhline(-np.log10(bh.max()), color="darkorange", lw=0.9, ls="-.",
                       label=f"BH q<0.05 ({(df.q<0.05).sum()} blocks)")
    ax.set_xticks(centers); ax.set_xticklabels(ticks)
    ax.set_ylabel(r"$-\log_{10}\,p$")
    ax.set_title(title, fontsize=11, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=7, loc="upper right", frameon=False)
    return df


def main():
    raw = pd.read_csv(f"{B}/block_gea.csv")
    wza = pd.read_csv(f"{B}/block_wza.csv")
    lr, lw = lam(raw.block_p.values), lam(wza.block_p.values)
    nfdr_r, nfdr_w = int((raw.q < 0.05).sum()), int((wza.q < 0.05).sum())

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=False)
    panel(axes[0], raw,
          f"Raw per-block test — 31 sites (varlen)   |   $\\lambda$={lr:.2f},  "
          f"FDR q<0.05: {nfdr_r} blocks")
    dfw = panel(axes[1], wza,
                f"Permutation-calibrated block-WZA — 31 sites (varlen)   |   "
                f"$\\lambda$={lw:.2f},  FDR q<0.05: {nfdr_w} blocks")
    axes[1].set_xlabel("genome position")
    # annotate the most extreme WZA block
    top = dfw.loc[dfw.nlp.idxmax()]
    axes[1].annotate(f"top block {top.unit}\np={top.block_p:.2g}, q={top.q:.2f}",
                     xy=(top.x, top.nlp), fontsize=7, color="black",
                     xytext=(8, 6), textcoords="offset points")
    fig.suptitle("Variable-length trajectories (all 31 GrENE-net sites) — "
                 "haploblock climate-GEA (bio1)", fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = f"{B}/manhattan_raw_vs_wza.png"
    fig.savefig(out, dpi=150)
    print(f"[done] wrote {out}")
    print(f"  raw : lambda={lr:.2f}  FDR<0.05={nfdr_r}")
    print(f"  WZA : lambda={lw:.2f}  FDR<0.05={nfdr_w}  "
          f"min block_p={wza.block_p.min():.4g} (q={wza.q.min():.3f})")


if __name__ == "__main__":
    main()

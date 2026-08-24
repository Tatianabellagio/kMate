#!/usr/bin/env python
"""3x2 (model x class) WZA Manhattan grid for the clq0.9 replication (deg-2 primary).

Reads clq90/wza/wza_{model}_{snp|nonsnp}_gen9_bio1_deg2.csv (block id in col 0 /
'gene', block genomic coord in retained chrom/pos, empirical block p Z_pVal), plots
-log10 p vs genome position per (model x class), marks the BH q<0.05 threshold and
the CAM5 clq0.9 blocks. Runs in the `basic` env (matplotlib present; the `plotting`
env hangs on import per project memory).

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
  $PY plot_clq90_manhattan.py --regime deg2
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

W = f"{lib.GEA}/phase1_replication/results/clq90/wza"
CHROM_LEN = {"Chr1": 30427671, "Chr2": 19698289, "Chr3": 23459830,
             "Chr4": 18585056, "Chr5": 26975502}
MODELS = ["kendall", "lfmm", "binomial"]
CLASSES = ["snp", "sv", "smallindel"]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = (p[o] * n) / (np.arange(n) + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]; return np.clip(q, 0, 1)


def offsets():
    off, cum = {}, 0
    for c in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
        off[c] = cum; cum += CHROM_LEN[c]
    return off, cum


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regime", default="isotonic")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    off, total = offsets()
    cam = pd.read_csv(f"{lib.CLQ_BLOCKS_DIR}/chr2_clq0.9_blocks_clq0.9.tsv", sep="\t")
    fig, axes = plt.subplots(len(MODELS), len(CLASSES), figsize=(15, 9), sharex=True)
    for i, model in enumerate(MODELS):
        for j, cls in enumerate(CLASSES):
            ax = axes[i, j]
            f = f"{W}/wza_{model}_{cls}_gen9_bio1_{args.regime}.csv"
            if not os.path.exists(f):
                ax.text(0.5, 0.5, f"missing\n{model} {cls}", ha="center", va="center",
                        transform=ax.transAxes); continue
            w = pd.read_csv(f)
            gc = "gene" if "gene" in w.columns else w.columns[0]
            w = w[w["Z_pVal"].notna()].copy()
            if "chrom" not in w.columns or "pos" not in w.columns:
                ax.text(0.5, 0.5, "no coord", ha="center", va="center",
                        transform=ax.transAxes); continue
            w = w[w["chrom"].isin(off)]
            w["gx"] = w["chrom"].map(off) + w["pos"]
            w["mlp"] = -np.log10(w["Z_pVal"].clip(lower=1e-300))
            w["q"] = bh(w["Z_pVal"].to_numpy())
            qthr = w.loc[w.q < 0.05, "Z_pVal"].max() if (w.q < 0.05).any() else None
            for k, c in enumerate(["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]):
                s = w[w.chrom == c]
                ax.scatter(s.gx, s.mlp, s=3, color=("#4c78a8", "#9ecae1")[k % 2],
                           rasterized=True, linewidths=0)
            # CAM5 blocks (Chr2 gene span) as a red marker
            camx = off["Chr2"] + 11532004
            ax.axvline(camx, color="crimson", lw=0.7, alpha=0.6)
            if qthr is not None:
                ax.axhline(-np.log10(qthr), color="red", ls="--", lw=0.8,
                           label=f"BH q<.05 ({int((w.q<0.05).sum())})")
                ax.legend(fontsize=7, loc="upper right")
            nsig = int((w.q < 0.05).sum())
            # No titles (repo convention): panel identity via an in-panel corner
            # annotation; model/class/regime context lives in the notebook markdown.
            ax.annotate(f"{model} · {cls}\n{len(w):,} blk · {nsig} BH-sig",
                        xy=(0.015, 0.97), xycoords="axes fraction",
                        ha="left", va="top", fontsize=8,
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
            ax.set_ylabel("-log10 p", fontsize=8)
    for j, cls in enumerate(CLASSES):
        axes[-1, j].set_xticks([off[c] + CHROM_LEN[c] / 2 for c in off])
        axes[-1, j].set_xticklabels(list(off), fontsize=8)
    fig.tight_layout()
    out = args.out or f"{lib.GEA}/phase1_replication/results/clq90/manhattan_clq90_{args.regime}.png"
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main()

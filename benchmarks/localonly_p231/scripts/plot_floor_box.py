#!/usr/bin/env python3
"""'How many k-mers does a local block need?' — distribution version.

Same data as plot_floor_curve.py but shows the FULL per-block RMSE distribution
per nnz bin as boxplots (box=IQR, whiskers=5-95 pct, line=median) instead of
median+p90 dots. One panel per unit."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

FD = Path("benchmarks/localonly_p231/floor_diag"); SIMS = Path("benchmarks/p231/sims")
RES = Path("benchmarks/localonly_p231/results")
KEYS = ["chrom", "pos", "ref_len", "alt_len"]
POOLS = ["cov10_n50_g0_s42_hotspots_p231_chr1",
         "cov10_n231_g1_s42_self97_hotspots_p231_chr1",
         "cov10_n50_g3_s42_hotspots_dom500nr_p231_chr1"]
EDGES = [1, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 100000]
LBL = ["1-7", "8-15", "16-31", "32-63", "64-127", "128-255", "256-511",
       "512-1k", "1k-2k", "2k-4k", "4k+"]


def per_block(pool, unit, min_rec=20):
    diag = pd.read_csv(FD / f"{pool}_{unit}.blockdiag.tsv", sep="\t")
    est = pd.read_csv(FD / f"{pool}_{unit}.recest.tsv", sep="\t")
    tr = pd.read_csv(SIMS / pool / "recomb_truth_raw.tsv.gz", sep="\t").dropna(subset=["truth_af"])
    tr["occ"] = tr.groupby(KEYS).cumcount(); est["occ"] = est.groupby(KEYS).cumcount()
    m = tr.merge(est[KEYS + ["occ", "alt_freq", "block"]], on=KEYS + ["occ"], how="inner")
    m = m.merge(diag[["block", "nnz"]], on="block", how="left")
    m = m[np.isfinite(m.alt_freq.values)].copy()
    m["e2"] = (m.alt_freq - m.truth_af) ** 2
    g = pd.DataFrame({"rmse": np.sqrt(m.groupby("block").e2.mean()),
                      "n": m.groupby("block").alt_freq.size(),
                      "nnz": m.groupby("block").nnz.first()})
    return g[g.n >= min_rec]


def boxes(ax, df, color):
    df = df.copy(); df["bi"] = pd.cut(df.nnz, EDGES, right=False, labels=False)
    data, pos, labs = [], [], []
    for i in range(len(EDGES) - 1):
        v = df.loc[df.bi == i, "rmse"].values
        if len(v) >= 10:
            data.append(v); pos.append(i); labs.append(f"{LBL[i]}\n(n={len(v)})")
    bp = ax.boxplot(data, positions=pos, widths=0.6, whis=(5, 95),
                    showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", lw=1.4))
    for b in bp["boxes"]:
        b.set(facecolor=color, alpha=0.55, edgecolor=color)
    for w in bp["whiskers"] + bp["caps"]:
        w.set(color=color)
    ax.set_xticks(pos); ax.set_xticklabels(labs, fontsize=7)
    return pos


fig, axes = plt.subplots(2, 1, figsize=(11, 8.5), sharey=True)
for ax, unit, col in [(axes[0], "w10kb", "C0"), (axes[1], "dynldK500", "C1")]:
    allb = pd.concat([per_block(p, unit) for p in POOLS], ignore_index=True)
    boxes(ax, allb, col)
    ax.axhline(0.048, color="gray", ls=":", lw=1)
    ax.text(0.1, 0.052, "accuracy floor ≈ 0.048", fontsize=8, color="gray")
    # bins 8..10 are 1k-2k, 2k-4k, 4k+ -> plateau region
    ax.axvspan(7.5, 10.5, color="green", alpha=0.07)
    # 50 falls in bin 3 (32-63)->4(64-127); 200 in bin 6 (128-255)->... mark boundaries
    ax.axvline(3.5, color="k", ls=":", lw=1); ax.text(3.55, 0.46, "bench 50", fontsize=8, rotation=90, va="top")
    ax.axvline(5.5, color="k", ls=":", lw=1); ax.text(5.55, 0.46, "default 200", fontsize=8, rotation=90, va="top")
    ax.set_ylim(0, 0.5); ax.set_ylabel("per-block AF RMSE")
    ax.set_title(f"{unit}", fontsize=10, loc="left")
axes[1].set_xlabel("nonzero k-mers in block (nnz)")
fig.suptitle("How many k-mers does a local block need? (per-block RMSE distribution)\n"
             "box = IQR, whiskers = 5–95 pct, line = median; green = plateau region",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(RES / "block_floor_box.png", dpi=140)
print(f"-> {RES/'block_floor_box.png'}")

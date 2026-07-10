#!/usr/bin/env python3
"""Distribution (strip/swarm) version of the panel & private-filter comparison.

One panel per series (p231 filt2inv / p80 filt2 / p80 unfiltered); per nnz bin,
jittered raw per-block RMSE points with the median bar. Centromere blocks in red
so the panel-quality effect on the centromere is visible. Pools the 3 scenarios."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RES = Path("benchmarks/localonly_p231/results")
KEYS = ["chrom", "pos", "ref_len", "alt_len"]
SCEN = ["cov10_n50_g0_s42_hotspots", "cov10_n231_g1_s42_self97_hotspots",
        "cov10_n50_g3_s42_hotspots_dom500nr"]
EDGES = [1, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 100000]
LBL = ["1-7","8-15","16-31","32-63","64-127","128-255","256-511","512-1k","1k-2k","2k-4k","4k+"]
RNG = np.random.RandomState(0)
# (label, arm_color, floor_diag dir, pool suffix, variant, sims dir, truth filename)
SERIES = [
    ("p231 filt2inv  (231 mixed-tech, private dropped)", "C0", "floor_diag", "_p231_chr1", "",
     "benchmarks/p231/sims", "recomb_truth_raw.tsv.gz"),
    ("p80 filt2  (80 long-read, private dropped)", "C1", "floor_diag_p80", "_p80_chr1", "filt2",
     "benchmarks/p80/sims", "recomb_truth.tsv.gz"),
    ("p80 unfiltered  (80 long-read, private KEPT)", "C2", "floor_diag_p80", "_p80_chr1", "unfiltered",
     "benchmarks/p80/sims", "recomb_truth.tsv.gz"),
]


def per_block(fd, pool, variant, sims, truthfn, min_rec=20):
    base = Path("benchmarks/localonly_p231") / fd
    tag = f"{pool}_w10kb" + (f"_{variant}" if variant else "")
    if not (base / f"{tag}.blockdiag.tsv").exists():
        return None
    diag = pd.read_csv(base / f"{tag}.blockdiag.tsv", sep="\t")
    est = pd.read_csv(base / f"{tag}.recest.tsv", sep="\t")
    tr = pd.read_csv(Path(sims) / pool / truthfn, sep="\t").dropna(subset=["truth_af"])
    tr["occ"] = tr.groupby(KEYS).cumcount(); est["occ"] = est.groupby(KEYS).cumcount()
    m = tr.merge(est[KEYS + ["occ", "alt_freq", "block"]], on=KEYS + ["occ"], how="inner")
    m = m.merge(diag[["block", "nnz", "start"]], on="block", how="left")
    m = m[np.isfinite(m.alt_freq.values)].copy(); m["e2"] = (m.alt_freq - m.truth_af) ** 2
    g = pd.DataFrame({"rmse": np.sqrt(m.groupby("block").e2.mean()),
                      "n": m.groupby("block").alt_freq.size(),
                      "nnz": m.groupby("block").nnz.first(),
                      "start": m.groupby("block").start.first()})
    g = g[g.n >= min_rec].copy()
    g["cen"] = (g.start >= 12_500_000) & (g.start <= 17_500_000)
    return g


def strip(ax, df, armcol, maxpts=500):
    df = df.copy(); df["bi"] = pd.cut(df.nnz, EDGES, right=False, labels=False)
    for i in range(len(EDGES) - 1):
        sub = df[df.bi == i]
        if len(sub) < 8:
            continue
        for mask, pcol in [(~sub.cen, armcol), (sub.cen, "crimson")]:
            vv = sub.loc[mask, "rmse"].values
            if not len(vv):
                continue
            show = vv if len(vv) <= maxpts else RNG.choice(vv, maxpts, replace=False)
            alpha = (0.5 if len(sub) <= maxpts else 0.3) if pcol != "crimson" else 0.7
            x = i + RNG.uniform(-0.28, 0.28, len(show))
            ax.scatter(x, show, s=7, color=pcol, alpha=alpha, edgecolors="none",
                       zorder=4 if pcol == "crimson" else 3)
        ax.plot([i - 0.34, i + 0.34], [sub.rmse.median()] * 2, color="black", lw=2, zorder=5)
    pos = list(range(len(EDGES) - 1))
    ax.set_xticks(pos)
    ax.set_xticklabels([f"{LBL[i]}\n(n={int((df.bi==i).sum())})" for i in pos], fontsize=6.5)
    ax.axhline(0.048, color="gray", ls=":", lw=1)
    ax.axvspan(7.5, 10.5, color="green", alpha=0.06)
    ax.axvline(5.5, color="k", ls=":", lw=1); ax.text(5.55, 0.47, "200", fontsize=7, rotation=90, va="top")
    ax.set_ylim(0, 0.5)


fig, axes = plt.subplots(3, 1, figsize=(11, 11), sharey=True)
for ax, (lab, col, fd, suf, var, sims, tfn) in zip(axes, SERIES):
    parts = [per_block(fd, s + suf, var, sims, tfn) for s in SCEN]
    df = pd.concat([p for p in parts if p is not None], ignore_index=True)
    strip(ax, df, col)
    ax.set_ylabel("per-block AF RMSE"); ax.set_title(lab, fontsize=10, loc="left")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=col, label="chromosome arm"),
                       Line2D([], [], marker="o", ls="", color="crimson", label="centromere 12.5–17.5 Mb")],
              fontsize=7.5, loc="upper right")
axes[2].set_xlabel("nonzero k-mers in block (nnz)")
fig.suptitle("Per-block AF RMSE distribution by panel & private-k-mer filter (w10kb)\n"
             "each dot = one block; black bar = median; big bins subsampled to 500",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(RES / "panel_compare_strip.png", dpi=140)
print(f"-> {RES/'panel_compare_strip.png'}")

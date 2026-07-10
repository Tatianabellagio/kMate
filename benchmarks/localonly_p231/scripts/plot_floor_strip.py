#!/usr/bin/env python3
"""Raw-data (swarm/strip) version of the floor plot + diagnose the w10kb 16-31 bump.

Per nnz bin, jittered raw per-block RMSE points (big bins subsampled for display)
with the median overlaid. Then investigate WHY the low-nnz w10kb bins are
non-monotonic: are they noise (small n) or a real region (centromere/repeat)?"""
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
LBL = ["1-7","8-15","16-31","32-63","64-127","128-255","256-511","512-1k","1k-2k","2k-4k","4k+"]
RNG = np.random.RandomState(0)


def per_block(pool, unit, min_rec=20):
    diag = pd.read_csv(FD / f"{pool}_{unit}.blockdiag.tsv", sep="\t")
    est = pd.read_csv(FD / f"{pool}_{unit}.recest.tsv", sep="\t")
    tr = pd.read_csv(SIMS / pool / "recomb_truth_raw.tsv.gz", sep="\t").dropna(subset=["truth_af"])
    tr["occ"] = tr.groupby(KEYS).cumcount(); est["occ"] = est.groupby(KEYS).cumcount()
    m = tr.merge(est[KEYS + ["occ", "alt_freq", "block"]], on=KEYS + ["occ"], how="inner")
    m = m.merge(diag[["block", "nnz", "start", "end"]], on="block", how="left")
    m = m[np.isfinite(m.alt_freq.values)].copy()
    m["e2"] = (m.alt_freq - m.truth_af) ** 2
    g = pd.DataFrame({"rmse": np.sqrt(m.groupby("block").e2.mean()),
                      "n": m.groupby("block").alt_freq.size(),
                      "nnz": m.groupby("block").nnz.first(),
                      "start": m.groupby("block").start.first(),
                      "truth_mean": m.groupby("block").truth_af.mean()})
    g["pool"] = pool
    g["cen"] = (g.start >= 12_500_000) & (g.start <= 17_500_000)
    return g[g.n >= min_rec]


# ---------- diagnosis: w10kb low-nnz bins ----------
w = pd.concat([per_block(p, "w10kb") for p in POOLS], ignore_index=True)
w["bi"] = pd.cut(w.nnz, EDGES, right=False, labels=False)
# Arabidopsis Chr1 (peri)centromere ~ 12.5-17.5 Mb
w["cen"] = (w.start >= 12_500_000) & (w.start <= 17_500_000)
print("=== w10kb low-nnz bins: is the 16-31 bump real? ===")
print(f"{'bin':>8} {'n':>5} {'med':>6} {'p90':>6} {'max':>6} {'%cen':>6} {'medAFtruth':>10}")
for i in range(6):
    s = w[w.bi == i]
    if len(s):
        print(f"{LBL[i]:>8} {len(s):>5} {s.rmse.median():>6.3f} {s.rmse.quantile(.9):>6.3f} "
              f"{s.rmse.max():>6.3f} {100*s.cen.mean():>5.0f}% {s.truth_mean.median():>10.3f}")
print("\nworst 8 windows in the 16-31 bin (w10kb):")
b2 = w[w.bi == 2].sort_values("rmse", ascending=False).head(8)
for _, r in b2.iterrows():
    print(f"  {r.pool[:22]:22} start={int(r.start):>9,}  nnz={int(r.nnz):>3}  "
          f"n={int(r.n):>4}  rmse={r.rmse:.3f}  truthAF_mean={r.truth_mean:.3f}  "
          f"{'CENTROMERE' if r.cen else ''}")


# ---------- strip/swarm figure ----------
def strip(ax, df, col, maxpts=500):
    df = df.copy(); df["bi"] = pd.cut(df.nnz, EDGES, right=False, labels=False)
    for i in range(len(EDGES) - 1):
        sub = df[df.bi == i]
        if len(sub) < 10:
            continue
        # plot arm and centromere points separately so the centromere shows up red
        for mask, pcol, a in [(~sub.cen, col, None), (sub.cen, "crimson", None)]:
            vv = sub.loc[mask, "rmse"].values
            if not len(vv):
                continue
            show = vv if len(vv) <= maxpts else RNG.choice(vv, maxpts, replace=False)
            alpha = (0.5 if len(sub) <= maxpts else 0.3) if pcol != "crimson" else 0.7
            x = i + RNG.uniform(-0.28, 0.28, len(show))
            ax.scatter(x, show, s=8, color=pcol, alpha=alpha, edgecolors="none",
                       zorder=4 if pcol == "crimson" else 3)
        ax.plot([i - 0.34, i + 0.34], [sub.rmse.median()] * 2, color="black", lw=2, zorder=5)
    pos = list(range(len(EDGES) - 1))
    ax.set_xticks(pos)
    ax.set_xticklabels([f"{LBL[i]}\n(n={int((df.bi==i).sum())})" for i in pos], fontsize=7)
    ax.axhline(0.048, color="gray", ls=":", lw=1)
    ax.axvspan(7.5, 10.5, color="green", alpha=0.07)
    ax.axvline(3.5, color="k", ls=":", lw=1); ax.text(3.55, 0.47, "bench 50", fontsize=8, rotation=90, va="top")
    ax.axvline(5.5, color="k", ls=":", lw=1); ax.text(5.55, 0.47, "default 200", fontsize=8, rotation=90, va="top")
    ax.set_ylim(0, 0.5)


from matplotlib.lines import Line2D
fig, axes = plt.subplots(2, 1, figsize=(11, 8.5), sharey=True)
for ax, unit, col in [(axes[0], "w10kb", "C0"), (axes[1], "dynldK500", "C1")]:
    allb = pd.concat([per_block(p, unit) for p in POOLS], ignore_index=True)
    strip(ax, allb, col)
    ax.set_ylabel("per-block AF RMSE"); ax.set_title(unit, fontsize=10, loc="left")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=col, label="chromosome arm"),
                       Line2D([], [], marker="o", ls="", color="crimson",
                              label="Chr1 (peri)centromere 12.5–17.5 Mb")],
              fontsize=8, loc="upper right")
axes[1].set_xlabel("nonzero k-mers in block (nnz)")
fig.suptitle("Raw per-block RMSE — low-nnz blocks are the centromere (red)\n"
             "each dot = one block; black bar = median; big bins subsampled to 500",
             fontsize=10.5)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(RES / "block_floor_strip.png", dpi=140)
print(f"\n-> {RES/'block_floor_strip.png'}")

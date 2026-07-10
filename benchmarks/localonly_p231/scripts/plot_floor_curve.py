#!/usr/bin/env python3
"""Final figure: per-block AF error vs k-mer supply — where does it plateau?

The decision-relevant curve. Pools the 3 representative sims; for each unit bins
blocks by nonzero-k-mer count and plots median + p90 block RMSE. Marks the
current floors (50, 200) and the accuracy floor. Answers "how many k-mers does
a local block need" directly.
"""
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


def per_block(pool, unit):
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
    return g[g.n >= 20]


fig, ax = plt.subplots(figsize=(8.2, 5.6))
ctr = [np.sqrt(EDGES[i] * EDGES[i + 1]) for i in range(len(EDGES) - 1)]
ctr[-1] = 6000
for unit, col in [("w10kb", "C0"), ("dynldK500", "C1")]:
    allb = pd.concat([per_block(p, unit) for p in POOLS], ignore_index=True)
    allb["bin"] = pd.cut(allb.nnz, EDGES, right=False)
    med = allb.groupby("bin", observed=True).rmse.median()
    p90 = allb.groupby("bin", observed=True).rmse.quantile(.9)
    x = ctr[:len(med)]
    ax.plot(x, med.values, "-o", color=col, label=f"{unit} median")
    ax.plot(x, p90.values, "--^", color=col, alpha=.55, label=f"{unit} p90")
floor = 0.048
ax.axhline(floor, color="gray", ls=":", lw=1)
ax.text(1.2, floor + 0.003, f"accuracy floor ≈ {floor:.03f}", fontsize=8, color="gray")
for v, lab in [(50, "bench 50"), (200, "default 200")]:
    ax.axvline(v, color="k", ls=":", lw=1)
    ax.text(v * 1.05, 0.40, lab, rotation=90, fontsize=8, va="top")
ax.axvspan(1024, 4096, color="green", alpha=0.06)
ax.text(1500, 0.36, "error plateaus\n(~1–2k k-mers)", fontsize=8, color="green")
ax.set_xscale("log"); ax.set_xlabel("nonzero k-mers in block (nnz)")
ax.set_ylabel("per-block AF RMSE (local-only, no smoothing)")
ax.set_title("How many k-mers does a local block need?\n"
             "Per-block AF error vs k-mer supply (pooled 3 p231 sims, Chr1)")
ax.legend(fontsize=8, ncol=2)
ax.set_ylim(0, 0.42)
fig.tight_layout()
fig.savefig(RES / "block_floor_curve.png", dpi=140)
print(f"-> {RES/'block_floor_curve.png'}")

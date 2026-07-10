#!/usr/bin/env python3
"""Aggregate per-block k-mer-floor rows -> accuracy-vs-(k-mers-in-block) curves.

Reads benchmarks/ecotype_count/results/kmer_floor/*.tsv, bins blocks by the
number of nonzero-count k-mers in the block, and plots per-block ecotype-mixture
recovery (cosine to truth, top-N Jaccard, mass-on-true) as a function of that
k-mer count, stratified by mixture size (true_n) and coverage. Also reports the
"floor": the smallest k-mer count at which median accuracy reaches 90% of the
high-k-mer asymptote, per mixture size.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/global/scratch/users/tbellg/kmate")
DDIR = ROOT / "benchmarks/ecotype_count/results/kmer_floor"
OUT = ROOT / "benchmarks/ecotype_count/results/kmer_floor_summary"

# log-spaced k-mer bins (extended to 100k for the window-size sweep)
EDGES = np.array([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000,
                  10000, 20000, 50000, 1e9])
CENTERS = np.sqrt(EDGES[:-1] * np.minimum(EDGES[1:], EDGES[-2] * 2))


def load():
    fs = sorted(DDIR.glob("*.tsv"))
    if not fs:
        sys.exit(f"no per-pool tsvs in {DDIR}")
    df = pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True)
    df["kbin"] = pd.cut(df["n_nz_kmers"], bins=EDGES, right=False, labels=False)
    print(f"loaded {len(df):,} local-fit blocks from {len(fs)} pools "
          f"({sorted(df['pool'].unique())})", file=sys.stderr)
    return df


def floor_of(sub, metric, frac=0.9):
    """Smallest bin-center where median(metric) >= frac * asymptote.
    Asymptote = median over the highest-k-mer decile actually present."""
    g = sub.groupby("kbin")[metric].median()
    if g.empty:
        return np.nan
    hi = sub["n_nz_kmers"].quantile(0.9)
    asym = sub[sub["n_nz_kmers"] >= hi][metric].median()
    if not np.isfinite(asym):
        asym = g.iloc[-1]
    target = frac * asym
    for kb in sorted(g.index):
        if g[kb] >= target:
            return CENTERS[int(kb)]
    return np.nan


def main():
    df = load()
    metrics = [("cosine", "cosine sim to true mixture"),
               ("jaccard_topN", "top-N Jaccard (identity)"),
               ("mass_on_true", "mass on true ecotypes")]

    # main figure: cov10 mixture sweep
    cov10 = df[df["cov"] == 10]
    ns = sorted(cov10["true_n"].unique())
    cmap = plt.cm.viridis(np.linspace(0, 0.9, len(ns)))

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, (mt, lab) in zip(axes, metrics):
        for c, n in zip(cmap, ns):
            sub = cov10[cov10["true_n"] == n]
            g = sub.groupby("kbin")[mt].agg(["median", "count"])
            x = [CENTERS[int(k)] for k in g.index]
            ax.plot(x, g["median"], "o-", color=c, label=f"n={n}")
        ax.set_xscale("log")
        ax.set_xlabel("nonzero k-mers in block")
        ax.set_ylabel(f"median {lab}")
        ax.set_title(lab)
        ax.grid(alpha=0.3)
        ax.axvline(200, ls="--", color="grey", lw=1)  # old fallback threshold
    axes[0].legend(title="mixture size", fontsize=8)
    fig.suptitle("Block-mode ecotype resolution vs k-mers per block "
                 "(fallback OFF, anchor OFF; cov10, g0)  — grey dashes = old min_kmers=200",
                 y=1.03)
    fig.tight_layout()
    fig.savefig(str(OUT) + "_cov10.png", dpi=130, bbox_inches="tight")

    # coverage contrast for n80
    n80 = df[df["true_n"] == 80]
    if n80["cov"].nunique() > 1:
        fig2, axes2 = plt.subplots(1, 3, figsize=(16, 4.8))
        for ax, (mt, lab) in zip(axes2, metrics):
            for cov in sorted(n80["cov"].unique()):
                sub = n80[n80["cov"] == cov]
                g = sub.groupby("kbin")[mt].median()
                x = [CENTERS[int(k)] for k in g.index]
                ax.plot(x, g.values, "o-", label=f"cov{cov}")
            ax.set_xscale("log"); ax.set_xlabel("nonzero k-mers in block")
            ax.set_ylabel(f"median {lab}"); ax.set_title(lab); ax.grid(alpha=0.3)
        axes2[0].legend()
        fig2.suptitle("Coverage shifts the k-mer supply, not the per-k-mer floor (n80, g0)", y=1.03)
        fig2.tight_layout()
        fig2.savefig(str(OUT) + "_n80_coverage.png", dpi=130, bbox_inches="tight")

    # floor table
    print("\n=== k-mer floor (>=90% of high-k-mer asymptote) ===")
    rows = []
    for cov in sorted(df["cov"].unique()):
        for n in sorted(df[df["cov"] == cov]["true_n"].unique()):
            sub = df[(df["cov"] == cov) & (df["true_n"] == n)]
            rec = dict(cov=cov, true_n=n, n_blocks=len(sub),
                       median_nz=int(sub["n_nz_kmers"].median()))
            hi = sub["n_nz_kmers"].quantile(0.9)
            for mt, _ in metrics:
                rec[f"floor_{mt}"] = round(float(floor_of(sub, mt)), 1)
                rec[f"asym_{mt}"] = round(float(sub[sub["n_nz_kmers"] >= hi][mt].median()), 3)
            rec["max_nz"] = int(sub["n_nz_kmers"].max())
            rows.append(rec)
    ft = pd.DataFrame(rows)
    ft.to_csv(str(OUT) + "_floor_table.tsv", sep="\t", index=False)
    print(ft.to_string(index=False))
    print(f"\n[fig] {OUT}_cov10.png")
    print(f"[fig] {OUT}_n80_coverage.png")
    print(f"[tsv] {OUT}_floor_table.tsv")


if __name__ == "__main__":
    main()

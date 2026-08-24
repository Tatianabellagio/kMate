#!/usr/bin/env python3
"""Sequencing-depth distribution across the 2,415 evolved GrENE-Net pool-seq samples
(Table_S5 weighted-mean coverage, via qc_coverage_audit.csv). Marks 1x/10x since those
are the two depths kMate is benchmarked at (benchmarks/{p80,p231}).

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        analysis/grenenet_gea/qc/plot_coverage_distribution.py
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/global/scratch/users/tbellg/kmate"
CSV = f"{ROOT}/analysis/grenenet_gea/results/qc_coverage_audit.csv"
OUT = f"{ROOT}/analysis/grenenet_gea/notebooks/plots/coverage_distribution.png"

GREY = "#4d4d4d"
plt.rcParams.update({
    "text.color": GREY, "axes.labelcolor": GREY, "axes.titlecolor": GREY,
    "xtick.color": GREY, "ytick.color": GREY, "font.family": "sans-serif",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.spines.left": False, "axes.spines.bottom": False,
    "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.8,
    "axes.axisbelow": True, "xtick.bottom": False, "ytick.left": False,
})

df = pd.read_csv(CSV)
d = df["seq_depth"].dropna()
med = d.median()

fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.hist(d, bins=40, color="#4c78a8", edgecolor="white", alpha=0.9)
ax.axvline(med, color="crimson", ls="--", lw=1.4, label=f"median = {med:.1f}×  (n={len(d)})")
ax.axvline(1, color=GREY, ls=":", lw=1.3, label="1× (benchmarked)")
ax.axvline(10, color=GREY, ls="-.", lw=1.3, label="10× (benchmarked)")
ax.plot(d, np.full_like(d, -d.max() * 0.02), "|", color="k", alpha=0.15, ms=6, clip_on=False)
ax.set_xlabel("sequencing depth (Table_S5 weighted-mean coverage, ×)")
ax.set_ylabel("number of samples")
ax.margins(x=0.01)
ax.legend(fontsize=9, frameon=False)
fig.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print(f"saved {OUT}")
print(f"seq_depth: min={d.min():.2f}x  median={med:.2f}x  mean={d.mean():.2f}x  max={d.max():.2f}x  n={len(d)}")

#!/usr/bin/env python
"""Aggregate per-sample window-vs-global AF comparison CSVs into a cohort summary + plot."""
import os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "results/grenenet_gea/window_vs_global"
rows = [pd.read_csv(f) for f in sorted(glob.glob(f"{BASE}/per_sample/*.csv"))]
df = pd.concat(rows, ignore_index=True)
print(f"samples aggregated: {len(df)}")

hist_cols = [c for c in df.columns if c.startswith("hist_")]
edges = np.array([0, 1e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0])
htot = df[hist_cols].sum(axis=0).to_numpy()
hfrac = htot / htot.sum()

# cohort summary (variant-weighted means where it matters)
nv = df["n_variants"].to_numpy()
def wmean(col):
    return float(np.average(df[col], weights=nv))
summary = {
    "n_samples": len(df),
    "total_variant_comparisons": int(nv.sum()),
    "mean_abs_delta (var-wtd)": wmean("mean_abs_delta"),
    "rms_delta (var-wtd)": wmean("rms_delta"),
    "median pearson_r (across samples)": float(df["pearson_r"].median()),
    "min pearson_r": float(df["pearson_r"].min()),
    "frac_changed (var-wtd)": wmean("frac_changed"),
    "frac |d|>0.001 (var-wtd)": wmean("frac_delta_gt_1e-3"),
    "frac |d|>0.01  (var-wtd)": wmean("frac_delta_gt_1e-2"),
    "frac |d|>0.05  (var-wtd)": wmean("frac_delta_gt_5e-2"),
    "frac |d|>0.10  (var-wtd)": wmean("frac_delta_gt_1e-1"),
}
with open(f"{BASE}/SUMMARY.txt", "w") as fh:
    for k, v in summary.items():
        line = f"{k:38s}: {v:,.6g}" if isinstance(v, float) else f"{k:38s}: {v:,}"
        print(line); fh.write(line + "\n")
    fh.write("\n|delta| histogram (fraction of all variant-comparisons):\n")
    for i in range(len(edges) - 1):
        fh.write(f"  [{edges[i]:.4g}, {edges[i+1]:.4g}) : {hfrac[i]:.5f}\n")
df.to_csv(f"{BASE}/per_sample_summary.csv", index=False)

# ---- plot ----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
labels = [f"{edges[i]:g}–{edges[i+1]:g}" for i in range(len(edges) - 1)]
ax[0].bar(range(len(hfrac)), hfrac, color="#3b6ea5")
ax[0].set_xticks(range(len(hfrac))); ax[0].set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
ax[0].set_ylabel("fraction of variant-comparisons"); ax[0].set_xlabel("|AF_window − AF_global|")
ax[0].set_title("AF-change magnitude distribution"); ax[0].set_yscale("log")

ax[1].hist(df["pearson_r"], bins=40, color="#3b6ea5")
ax[1].set_xlabel("per-sample Pearson r (window vs global)"); ax[1].set_ylabel("samples")
ax[1].set_title(f"correlation per sample (median {df['pearson_r'].median():.4f})")

for col, lab in [("frac_delta_gt_1e-2", "|d|>0.01"), ("frac_delta_gt_5e-2", "|d|>0.05"),
                 ("frac_delta_gt_1e-1", "|d|>0.10")]:
    ax[2].hist(df[col], bins=40, alpha=0.6, label=lab)
ax[2].set_xlabel("fraction of variants exceeding threshold"); ax[2].set_ylabel("samples")
ax[2].set_title("per-sample fraction shifted"); ax[2].legend()

fig.tight_layout()
fig.savefig(f"{BASE}/window_vs_global_af.png", dpi=130)
print(f"\nwrote {BASE}/SUMMARY.txt, per_sample_summary.csv, window_vs_global_af.png")

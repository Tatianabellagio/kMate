#!/usr/bin/env python3
"""Normalized views: does block AF error track raw k-mer COUNT, or k-mer DENSITY
(per variant / per kb)? Tests whether the nnz-collapse of dynld onto w10kb is
because total supply is the sufficient statistic (density irrelevant)."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

FD = Path("benchmarks/localonly_p231/floor_diag"); SIMS = Path("benchmarks/p231/sims")
RES = Path("benchmarks/localonly_p231/results")
KEYS = ["chrom", "pos", "ref_len", "alt_len"]
POOLS = ["cov10_n50_g0_s42_hotspots_p231_chr1",
         "cov10_n231_g1_s42_self97_hotspots_p231_chr1",
         "cov10_n50_g3_s42_hotspots_dom500nr_p231_chr1"]
SHORT = {POOLS[0]: "n50 g0 outcross", POOLS[1]: "n231 g1 SELFING", POOLS[2]: "n50 g3 dom-sel"}


def per_block(pool, unit, min_rec=20):
    diag = pd.read_csv(FD / f"{pool}_{unit}.blockdiag.tsv", sep="\t")
    est = pd.read_csv(FD / f"{pool}_{unit}.recest.tsv", sep="\t")
    tr = pd.read_csv(SIMS / pool / "recomb_truth_raw.tsv.gz", sep="\t").dropna(subset=["truth_af"])
    tr["occ"] = tr.groupby(KEYS).cumcount(); est["occ"] = est.groupby(KEYS).cumcount()
    m = tr.merge(est[KEYS + ["occ", "alt_freq", "block"]], on=KEYS + ["occ"], how="inner")
    m = m[np.isfinite(m.alt_freq.values)].copy()
    m["e2"] = (m.alt_freq - m.truth_af) ** 2
    g = pd.DataFrame({"rmse": np.sqrt(m.groupby("block").e2.mean()),
                      "n_var": m.groupby("block").alt_freq.size()})
    g = g.merge(diag.set_index("block")[["nnz", "start", "end", "n_present"]],
                left_index=True, right_index=True)
    g = g[g.n_var >= min_rec].copy()
    g["len_kb"] = (g.end - g.start + 1) / 1000.0
    g["kmers_per_var"] = g.nnz / g.n_var
    g["kmers_per_kb"] = g.nnz / g.len_kb
    g["unit"] = unit; g["pool"] = SHORT[pool]
    return g


allb = pd.concat([per_block(p, u) for p in POOLS for u in ["w10kb", "dynldK500"]],
                 ignore_index=True)

# ---- correlation comparison: which predictor tracks RMSE best? ----
METRICS = ["nnz", "kmers_per_var", "kmers_per_kb", "n_var", "len_kb", "n_present"]
print("=== Spearman ρ(RMSE, metric)  (more negative = predicts lower error) ===")
rows = []
for (pool, unit), sub in allb.groupby(["pool", "unit"]):
    rows.append(dict(pool=pool, unit=unit,
                     **{k: round(spearmanr(sub[k], sub.rmse).correlation, 3) for k in METRICS}))
corr = pd.DataFrame(rows)
print(corr.to_string(index=False))
corr.to_csv(RES / "block_density_corr.tsv", sep="\t", index=False)

# ---- partial check: within fixed nnz bands, does density still explain error? ----
print("\n=== within fixed nnz bands (w10kb+dynld pooled): ρ(RMSE, kmers_per_var) ===")
nbands = [(50, 100), (100, 200), (200, 400), (400, 800), (800, 1600)]
for lo, hi in nbands:
    s = allb[(allb.nnz >= lo) & (allb.nnz < hi)]
    if len(s) > 30:
        r_kv = spearmanr(s.kmers_per_var, s.rmse).correlation
        r_nz = spearmanr(s.nnz, s.rmse).correlation
        print(f"  nnz[{lo:>4},{hi:>4}) n={len(s):>5}  ρ(per_var)={r_kv:+.3f}  "
              f"ρ(nnz, residual supply)={r_nz:+.3f}")


def binned(df, xcol, q=10):
    df = df[df[xcol] > 0]
    edges = np.unique(np.quantile(df[xcol], np.linspace(0, 1, q + 1)))
    b = pd.cut(df[xcol], edges, include_lowest=True)
    x = df.groupby(b, observed=True)[xcol].median()
    med = df.groupby(b, observed=True).rmse.median()
    p90 = df.groupby(b, observed=True).rmse.quantile(.9)
    return x.values, med.values, p90.values


# ---- figure: RMSE vs the three views (count / per-variant / per-kb) ----
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
views = [("nnz", "nonzero k-mers (raw count)"),
         ("kmers_per_var", "k-mers per variant"),
         ("kmers_per_kb", "k-mers per kb (density)")]
for ax, (col, lab) in zip(axes, views):
    for unit, c in [("w10kb", "C0"), ("dynldK500", "C1")]:
        sub = allb[allb.unit == unit]
        x, med, p90 = binned(sub, col)
        ax.plot(x, med, "-o", color=c, label=f"{unit} median")
        ax.plot(x, p90, "--^", color=c, alpha=.55, label=f"{unit} p90")
    ax.axhline(0.048, color="gray", ls=":", lw=1)
    ax.set_xscale("log"); ax.set_xlabel(col + f"\n({lab})")
    ax.set_ylim(0, 0.30)
axes[0].set_ylabel("per-block AF RMSE (local-only, no smoothing)")
axes[0].legend(fontsize=8, ncol=2)
fig.suptitle("Raw count vs density: which collapses w10kb & dynld onto one curve?\n"
             "(if RAW COUNT is the sufficient statistic, only the left panel collapses)",
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(RES / "block_floor_normalized.png", dpi=140)
print(f"\n-> {RES/'block_floor_normalized.png'}")

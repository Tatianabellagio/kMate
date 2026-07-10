#!/usr/bin/env python3
"""True-vs-estimated AF hexbin panels for local-only on p231, one figure per unit.

For each of the 11 sim setups, join recomb_truth_raw <-> the local-only AF TSV
(same KEYS+occ join as the scorer), keep finite-both records, and draw a 2D hexbin
(log count) of truth_af (x) vs est alt_freq (y) with the y=x line. Title carries
R2/MAE/call-rate. Selfing (self97) panels — the GrENE-Net-relevant regime — get a
highlighted title. One PNG per unit {dynldK500, w10kb}.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/global/scratch/users/tbellg/kmate")
WORK = ROOT / "benchmarks/localonly_p231/work"
SIMS = ROOT / "benchmarks/p231/sims"
OUTD = ROOT / "benchmarks/localonly_p231/results"
KEYS = ["chrom", "pos", "ref_len", "alt_len"]

POOLS = [
    "cov10_n50_g0_s42_hotspots_p231_chr1",
    "cov10_n231_g0_s42_hotspots_p231_chr1",
    "cov10_n50_g1_s42_hotspots_p231_chr1",
    "cov10_n231_g1_s42_hotspots_p231_chr1",
    "cov10_n50_g1_s42_self97_hotspots_p231_chr1",
    "cov10_n231_g1_s42_self97_hotspots_p231_chr1",
    "cov10_n50_g3_s42_hotspots_p231_chr1",
    "cov10_n50_g3_s42_self97_hotspots_p231_chr1",
    "cov10_n50_g3_s42_hotspots_dom500_p231_chr1",
    "cov10_n50_g3_s42_hotspots_dom500nr_p231_chr1",
    "cov10_n50_g3_s42_self97_hotspots_dom500_p231_chr1",
]


def label(pool):
    n = "n231" if "_n231_" in pool else "n50"
    g = "g0" if "_g0_" in pool else ("g1" if "_g1_" in pool else "g3")
    mate = "SELFING" if "self97" in pool else "outcross"
    sel = "+domNR" if "dom500nr" in pool else ("+dom" if "dom500" in pool else "")
    return f"{n} {g} {mate}{sel}", ("self97" in pool)


def load_join(pool, unit):
    est = pd.read_csv(WORK / f"{pool}_localonly_{unit}.tsv", sep="\t")
    est["occ"] = est.groupby(KEYS).cumcount()
    tr = pd.read_csv(SIMS / pool / "recomb_truth_raw.tsv.gz", sep="\t").dropna(subset=["truth_af"])
    tr["occ"] = tr.groupby(KEYS).cumcount()
    m = tr.merge(est[KEYS + ["occ", "alt_freq"]], on=KEYS + ["occ"], how="inner")
    t = m["truth_af"].values.astype(float)
    e = m["alt_freq"].values.astype(float)
    fin = np.isfinite(e) & np.isfinite(t)
    call = fin.sum() / len(m) if len(m) else np.nan
    return t[fin], e[fin], call


def panel(ax, t, e, call, title, is_self):
    hb = ax.hexbin(t, e, gridsize=70, bins="log", cmap="viridis", mincnt=1, extent=(0, 1, 0, 1))
    ax.plot([0, 1], [0, 1], color="red", lw=0.8, ls="--", alpha=0.8)
    d = e - t
    ss = np.sum((t - t.mean()) ** 2)
    r2 = 1 - np.sum(d ** 2) / ss if ss > 0 else np.nan
    mae = np.abs(d).mean()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title(f"{title}\ncall {call*100:.0f}%  R²={r2:.3f}  MAE={mae:.3f}",
                 fontsize=9, color=("darkgreen" if is_self else "black"),
                 fontweight=("bold" if is_self else "normal"))
    ax.tick_params(labelsize=7)
    return hb


def make_fig(unit, unit_label):
    fig, axes = plt.subplots(3, 4, figsize=(15, 11))
    axes = axes.ravel()
    last = None
    for i, pool in enumerate(POOLS):
        t, e, call = load_join(pool, unit)
        title, is_self = label(pool)
        last = panel(axes[i], t, e, call, title, is_self)
        if i % 4 == 0:
            axes[i].set_ylabel("estimated AF", fontsize=8)
        if i >= 7:
            axes[i].set_xlabel("true AF", fontsize=8)
        print(f"  {pool}: n={len(t):,} call={call*100:.0f}%", file=sys.stderr)
    for j in range(len(POOLS), len(axes)):
        axes[j].axis("off")
    fig.suptitle(f"local-only · {unit_label} — true vs estimated AF (p231, 10×, all sims; "
                 f"green/bold = selfing = GrENE-Net regime)", fontsize=13)
    if last is not None:
        cax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
        fig.colorbar(last, cax=cax, label="log10(record count)")
    fig.tight_layout(rect=[0, 0, 0.92, 0.97])
    out = OUTD / f"true_vs_est_localonly_{unit}.png"
    fig.savefig(out, dpi=130)
    print(f"wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    make_fig("dynldK500", "dynld_K500 blocks")
    make_fig("w10kb", "10 kb windows")

#!/usr/bin/env python3
"""Classic truth-vs-estimated AF scatter for block-mode kMate WITHOUT the global
fallback, compared against the global-mode run. All per-record TSVs share the
var_pa meta row order (full Chr1 panel), so they align to recomb_truth by row.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/global/scratch/users/tbellg/kmate")
POOL = "cov10_n80_g0_s42_hotspots_p80_chr1"
W = ROOT / "benchmarks/accuracy_vs_competitors/work"
TRUTH = ROOT / f"benchmarks/p80/sims/{POOL}/recomb_truth.tsv.gz"
GLOBAL = ROOT / f"benchmarks/benchmark_runs/tsv/{POOL}_global.tsv"
OUT = ROOT / f"benchmarks/accuracy_vs_competitors/results/block_nofallback_af_{POOL}"

PANELS = [
    ("global (production)", GLOBAL),
    ("block, NO fallback — pure local\n(anchor 0, no smoothing)", W / f"{POOL}__blocknofb_purelocal.tsv"),
    ("block, NO fallback — anchored\n(anchor 0.3, smooth 5)", W / f"{POOL}__blocknofb_anchored.tsv"),
]


def metrics(t, e):
    m = np.isfinite(t) & np.isfinite(e)
    t, e = t[m], e[m]
    err = e - t
    ss_res = float((err ** 2).sum()); ss_tot = float(((t - t.mean()) ** 2).sum())
    return dict(n=int(m.sum()), MAE=float(np.abs(err).mean()),
                R2=(1 - ss_res / ss_tot if ss_tot > 0 else np.nan),
                r=float(np.corrcoef(t, e)[0, 1]))


tr = pd.read_csv(TRUTH, sep="\t")
truth = tr["truth_af"].values
svmask = (np.abs(tr["alt_len"] - tr["ref_len"]) >= 50).values
print(f"truth: {len(tr):,} records ({svmask.sum():,} SV)", file=sys.stderr)

TAGS = ["global", "purelocal", "anchored"]


def draw(ax, lab, path):
    if not path.exists():
        ax.text(0.5, 0.5, f"missing:\n{path.name}", ha="center", va="center")
        ax.set_title(lab); return
    d = pd.read_csv(path, sep="\t")
    if len(d) != len(tr):
        ax.text(0.5, 0.5, f"len {len(d):,} != truth {len(tr):,}", ha="center")
        ax.set_title(lab); return
    est = d["alt_freq"].values
    hb = ax.hexbin(truth, est, gridsize=80, bins="log", cmap="viridis", mincnt=1)
    ax.plot([0, 1], [0, 1], "r--", lw=1.2)
    mt = metrics(truth, est); mtsv = metrics(truth[svmask], est[svmask])
    ax.set_title(f"{lab}\nALL  r={mt['r']:.3f}  R²={mt['R2']:.3f}  MAE={mt['MAE']:.4f}\n"
                 f"SV   r={mtsv['r']:.3f}  R²={mtsv['R2']:.3f}  MAE={mtsv['MAE']:.4f}",
                 fontsize=12)
    ax.set_xlabel("truth alt-AF"); ax.set_ylabel("kMate estimated alt-AF")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    print(f"{lab.splitlines()[0]:40s} ALL r={mt['r']:.3f} R2={mt['R2']:.3f} "
          f"MAE={mt['MAE']:.4f} | SV r={mtsv['r']:.3f} R2={mtsv['R2']:.3f}", file=sys.stderr)
    return hb


# combined (larger) + one readable single-panel PNG per approach
fig, axes = plt.subplots(1, len(PANELS), figsize=(7.2 * len(PANELS), 7))
for ax, (lab, path) in zip(axes, PANELS):
    draw(ax, lab, path)
fig.suptitle(f"kMate AF: block mode WITHOUT global fallback vs global  —  {POOL}", y=1.02, fontsize=14)
fig.tight_layout()
fig.savefig(str(OUT) + ".png", dpi=130, bbox_inches="tight")
print(f"[fig] {OUT}.png")

for tag, (lab, path) in zip(TAGS, PANELS):
    f1, a1 = plt.subplots(figsize=(7, 7))
    hb = draw(a1, lab, path)
    if hb is not None:
        f1.colorbar(hb, ax=a1, label="log10(count)")
    f1.tight_layout()
    f1.savefig(f"{OUT}__{tag}.png", dpi=140, bbox_inches="tight")
    print(f"[fig] {OUT}__{tag}.png")

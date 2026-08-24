#!/usr/bin/env python
"""Manhattan of WITHIN-SITE selection coefficients s_g for ONE site (default site 4).

For the chosen site, reconstruct the Step-3 weighted-least-squares slope s_g (and its
model SE) of logit(pooled hap freq) over generation, for every testable haplotype, and
draw it along the genome. NO climate, NO permutation null here — this is the per-site
temporal-change signal itself (descriptive): which haplotypes rose/fell at this garden.

Two panels:
  top    = signed s_g (effect size; + = rose over gen 0->3, - = fell)
  bottom = -log10(Wald p) from s_g/SE(s_g), coloured by direction
Wald p is MODEL-BASED (4 points, known V, 0 residual df) → read as a descriptive
z-score, not a calibrated test.

Env: kmate. SITE / PB_DIR via env. Writes site<ID>_scoef_manhattan.png into PB_DIR.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

H = "results/grenenet_gea/hapfreq"
B = os.environ.get("PB_DIR", f"{H}/pipelineB")
SITE = int(os.environ.get("SITE", 4))
EPS = 1e-3
Tg = np.array([0.0, 1.0, 2.0, 3.0])
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
COLORS = ["#cccccc", "#dddddd"]


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def main():
    P = np.load(f"{B}/traj_p.npy").astype(np.float64)        # (S,4,nhap)
    V = np.load(f"{B}/traj_v.npy").astype(np.float64)
    sites = [int(s) for s in np.load(f"{B}/traj_sites.npy")]
    gi = sites.index(SITE)
    reg = pd.read_csv(f"{H}/hapfreq_registry.csv")

    # Step-3 WLS slope for THIS site only
    y = logit(P[gi]); w = 1.0 / np.clip(V[gi], 1e-12, None)   # (4,nhap)
    sw = w.sum(0); tbar = (w * Tg[:, None]).sum(0) / sw
    dt = Tg[:, None] - tbar[None, :]
    Sxx = (w * dt**2).sum(0)
    s_g = (w * dt * y).sum(0) / np.clip(Sxx, 1e-12, None)
    se_g = np.sqrt(1.0 / np.clip(Sxx, 1e-12, None))
    z = s_g / np.clip(se_g, 1e-12, None)
    p = 2 * stats.norm.sf(np.abs(z))

    d = reg[["chrom", "unit_start", "unit_end", "covered", "panel_freq"]].copy()
    d["s"] = s_g; d["p"] = p; d["z"] = z
    d = d[d.covered & d.panel_freq.between(0.05, 0.95)].copy()     # testable haplotypes
    d["mid"] = (d.unit_start + d.unit_end) / 2
    d["nlp"] = -np.log10(np.clip(d.p, 1e-300, 1))

    # genome-cumulative x
    off, centers, x = 0.0, [], np.zeros(len(d))
    d = d.reset_index(drop=True)
    for ch in CHROMS:
        m = (d.chrom == ch).to_numpy()
        if not m.any():
            continue
        cmax = d.loc[m, "mid"].max()
        x[m] = d.loc[m, "mid"] + off
        centers.append(off + cmax / 2); off += cmax * 1.02
    d["x"] = x
    up = d.s > 0

    fig, (a0, a1) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    # chrom shading
    for ax in (a0, a1):
        for i, ch in enumerate(CHROMS):
            m = (d.chrom == ch).to_numpy()
            ax.scatter(d.x[m], np.zeros(m.sum()), s=0)  # keep extent
    a0.scatter(d.x[up], d.s[up], s=5, c="#c0392b", alpha=.5, edgecolors="none",
               rasterized=True, label="rose (s>0)")
    a0.scatter(d.x[~up], d.s[~up], s=5, c="#2471a3", alpha=.5, edgecolors="none",
               rasterized=True, label="fell (s<0)")
    a0.axhline(0, color="k", lw=.6)
    a0.set_ylabel(r"within-site slope $s_g$  (log-odds / gen)")
    a0.set_title(f"Site {SITE} — per-haplotype selection coefficient $s_g$ "
                 f"(gen 0$\\to$3, {len(d):,} testable haplotypes)", loc="left", fontsize=11)
    a0.legend(fontsize=8, frameon=False, loc="upper right")

    a1.scatter(d.x[up], d.nlp[up], s=5, c="#c0392b", alpha=.5, edgecolors="none", rasterized=True)
    a1.scatter(d.x[~up], d.nlp[~up], s=5, c="#2471a3", alpha=.5, edgecolors="none", rasterized=True)
    a1.set_ylabel(r"$-\log_{10}\,p$  (Wald $s_g/\mathrm{SE}$)")
    a1.set_title("Within-site temporal-change significance (model-based Wald; "
                 "descriptive, no permutation null)", loc="left", fontsize=10)
    a1.set_xticks(centers); a1.set_xticklabels(CHROMS)
    a1.set_xlabel("genome position")
    for ax in (a0, a1):
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    out = f"{B}/site{SITE}_scoef_manhattan.png"
    fig.savefig(out, dpi=150)
    print(f"[done] {out}")
    print(f"  site {SITE}: {len(d):,} testable haplotypes | "
          f"s_g median {d.s.median():+.3f}, IQR [{d.s.quantile(.25):+.3f},{d.s.quantile(.75):+.3f}]")
    print(f"  |s_g|>0.5: {(d.s.abs()>0.5).sum():,}  ({(d.s>0.5).sum()} up / {(d.s<-0.5).sum()} down)")
    top = d.reindex(d.s.abs().sort_values(ascending=False).index).head(8)
    print("  strongest movers:")
    print(top[["chrom", "unit_start", "unit_end", "panel_freq", "s", "z", "p"]].to_string(index=False))


if __name__ == "__main__":
    main()

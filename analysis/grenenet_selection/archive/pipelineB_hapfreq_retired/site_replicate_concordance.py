#!/usr/bin/env python
"""Drift-vs-selection at ONE site via REPLICATE CONCORDANCE (default site 4).

Reconstructs the per-PLOT haplotype trajectories the pooled-mean pipeline averages
away, then asks the experimental-evolution question directly: for each haplotype, did
the replicate plots move the SAME direction (gen0 anchor -> last gen)?  Under drift the
plots scatter (sign-concordance ~ binomial(n,1/2)); under selection they agree.

Shows whether the pooled-mean "strong movers" (|s_g|>0.5) survive a concordance
requirement -- i.e. how many were drift that the mean let through.

Env: kmate. SITE / MODE via env.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

H = "results/grenenet_gea/hapfreq"
SITE = int(os.environ.get("SITE", 4))
EPS = 1e-3


def logit(p):
    return np.log(np.clip(p, EPS, 1 - EPS) / (1 - np.clip(p, EPS, 1 - EPS)))


def main():
    mat = np.load(f"{H}/hapfreq_matrix.npy")
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    p0 = np.load(f"{H}/hapfreq_p0_seedmix.npy").astype(np.float64)   # shared gen-0 anchor
    reg = pd.read_csv(f"{H}/hapfreq_registry.csv")

    pt = lib.pool_table()
    pt = pt[(pt.site == SITE) & pt.sampleid.astype(str).isin(sidx)].copy()
    pt["row"] = pt.sampleid.astype(str).map(sidx)

    # per (gen,plot): flower-weighted hapfreq across timepoints (mirror trajectory builder)
    cell = {}
    for (gen, plot), g in pt.groupby(["generation", "plot"]):
        w = g.flowerscollected.to_numpy(float)
        w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        cell[(int(gen), int(plot))] = (mat[g.row.to_numpy()] * w[:, None]).sum(0) / w.sum()

    # persistent plots: present at all of gen 1,2,3
    present = {}
    for (gen, plot) in cell:
        present.setdefault(plot, set()).add(gen)
    plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
    n = len(plots)
    print(f"site {SITE}: {n} persistent plots -> {plots}")

    # per-plot net change on the log-odds scale, gen0(anchor) -> gen3
    g3 = np.vstack([cell[(3, pl)] for pl in plots])            # (n, nhap)
    dlast = logit(g3) - logit(p0)[None, :]                     # per-plot Δlogit
    pooled = (g3.mean(0))                                      # equal-wt mean endpoint (illustrative)
    s_pool = logit(pooled) - logit(p0)                         # pooled-mean net change

    # sign concordance across plots
    nup = (dlast > 0).sum(0)
    conc = np.maximum(nup, n - nup) / n                        # fraction agreeing on direction
    allsame = (nup == 0) | (nup == n)
    # binomial null prob of n-of-n same direction by chance = 2*(1/2)^n
    p_allsame_null = 2 * 0.5**n

    testable = (reg.covered & reg.panel_freq.between(0.05, 0.95)).to_numpy()
    d = pd.DataFrame({
        "chrom": reg.chrom, "start": reg.unit_start, "end": reg.unit_end,
        "panel_freq": reg.panel_freq, "s_pool": s_pool, "conc": conc,
        "nup": nup, "allsame": allsame, "abs_dmean": np.abs(dlast).mean(0),
        "spread": dlast.std(0),
    })[testable].copy()

    strong = d.s_pool.abs() > 0.5
    print(f"\ntestable haplotypes: {len(d):,} | {n} plots, "
          f"chance of all-{n}-agree under drift = {p_allsame_null:.3f}")
    print(f"pooled |Δ|>0.5 ('strong movers'): {int(strong.sum()):,}")
    print(f"  of those, ALL {n} plots same direction: "
          f"{int((strong & d.allsame).sum()):,} "
          f"({100*(strong & d.allsame).mean()/max(strong.mean(),1e-9):.0f}% of strong movers)")
    print(f"  expected by chance if all strong movers were drift: "
          f"{p_allsame_null*100:.0f}%")
    print(f"\nconcordance of strong movers vs all testable (mean fraction agreeing):")
    print(f"  all testable : {d.conc.mean():.3f}")
    print(f"  strong movers: {d.conc[strong].mean():.3f}")

    print(f"\n--- the pooled-mean TOP movers: are the plots actually concordant? ---")
    top = d.reindex(d.s_pool.abs().sort_values(ascending=False).index).head(10)
    show = top[["chrom", "start", "end", "panel_freq", "s_pool", "nup", "conc"]].copy()
    show["plots_dir"] = [f"{int(u)}/{n} up" for u in top.nup]
    print(show.to_string(index=False))

    out = f"{H}/pipelineB/site{SITE}_replicate_concordance.csv"
    d.to_csv(out, index=False)
    print(f"\n[done] {out}")


if __name__ == "__main__":
    main()

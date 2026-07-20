#!/usr/bin/env python
"""Per-founder, per-site SELECTION-COEFFICIENT trait for the variance-partition analysis.

The trait for the ecotype-selection GWAS / variance partition is the per-founder logit-slope
of its genome-wide founder frequency h over generations 0..3 -- a per-generation selection
coefficient s_{f,site}. This is numerically stable where raw h / Delta-h are NOT: h sums to 1
with a long tail of tiny values (skew ~+12, ~10 founders carry half the Delta-h variance);
logit-slope s is ~symmetric (skew ~+1) and spreads variance across founders. See the trait
decision in the session notes.

Design (matches ecotype_selection_site.py, extended to all sites, GLOBAL-mode h):
  * h source = GLOBAL-mode chrom-averaged founder h, reused from the sample_global_h.npz cache
    built by ecotype_fitness.py (no re-read of the ~10k per-chrom npz).
  * founding reference p0 = mean over the 8 SEEDMIX reps (chrom-averaged). Using the ESTIMATED
    seedmix p0 (not forced uniform 1/231) is deliberate: twin-absorbed founders read ~0 at both
    founding and every generation, so the slope cancels that identifiability bias to first order.
  * per site: per-plot trajectory h_{f,plot,gen} (timepoints collapsed to the site_gen_plot pool,
    flower-weighted), gen0 := p0; OLS slope of logit(h) over gens present; s_{f,site} = mean over
    plots. Requires the gen1 anchor. Plot replicates averaged (not pooled) so a plot with few
    flowers does not dominate.
  * ALL 231 founders kept analyzable (decision 2026-07-07): the ESTIMATED seed-mix p0 is FLOORED
    at P0_FLOOR (a small additive number, NOT forced uniform 1/231) so no ecotype disappears — the
    single hardest founder 9977 (p0~4e-6 under uniform ω) gets a defined near-zero slope instead of
    being dropped. Every other founder has p0 >> P0_FLOOR (floor is a no-op). Supersedes the old
    p0>1e-3 drop-floor that discarded ~19 legitimate rare-start founders (was 212 analyzable).
    NO reliability/cross-chrom weighting (chrom-averaging regularizes h).

Output -> analysis/grenenet_gea/varexp/selection_s_matrix.npz
  S[n_site x 231] logit-slope, sites, founders, bio1, p0, analyzable(bool 231), n_plots[n_site],
  freq_last[n_site x 231] (QC). Env: kmate. Light (runs off the cache in seconds).
"""
from __future__ import annotations
import os, sys, glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

OUT = f"{lib.GEA}/varexp"
HCACHE = f"{lib.GEA}/ecotype_fitness/sample_global_h.npz"
CHROMS = [f"Chr{i}" for i in range(1, 6)]
EPS = 1e-4            # logit clip
# p0 FLOOR (decision 2026-07-07): rather than DROP near-zero founders, floor the
# estimated seed-mix p0 at P0_FLOOR so no ecotype disappears — a small additive
# floor, NOT forcing uniform 1/231. Under the per_founder fix + uniform ω, the
# single hardest founder 9977 reads p0~4e-6 (uniform ω sacrifices it; ω=1/m_b kept
# it ~1.9e-4). Flooring p0 keeps ALL 231 founders analyzable: 9977 gets a defined
# (near-zero, ~flat) selection slope instead of being excluded. Every other founder
# has p0 >> P0_FLOOR so the floor is a no-op for them. (The old 1e-3 drop-floor
# discarded ~19 legitimate rare-start founders → only 212; we keep all 231.)
P0_FLOOR = EPS       # = 1e-4
TRAIT_GENS = (1, 2, 3)   # anchor gen1 required; use gens 1..3 (+ gen0=p0)


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def seedmix_p0():
    """Founding p0 = mean over the 8 SEEDMIX reps of the chrom-averaged global h."""
    ids = sorted({os.path.basename(p).split("_Chr")[0]
                  for p in glob.glob(f"{lib.SEEDMIX}/*_Chr1.h_per_chrom.npz")})
    founders, reps = None, []
    for s in ids:
        chs = []
        for c in CHROMS:
            z = np.load(f"{lib.SEEDMIX}/{s}_{c}.h_per_chrom.npz", allow_pickle=True)
            if founders is None:
                founders = z["founders"].astype(str)
            chs.append(z[c].astype(float))
        reps.append(np.vstack(chs).mean(0))
    return founders, np.vstack(reps).mean(0)


def main():
    os.makedirs(OUT, exist_ok=True)
    founders, p0 = seedmix_p0()
    p0 = np.maximum(p0, P0_FLOOR)   # floor so near-zero founders don't disappear (see P0_FLOOR)
    z = np.load(HCACHE, allow_pickle=True)
    assert (z["founders"].astype(str) == founders).all(), "founder order mismatch"
    Hs = z["H"].astype(float)                                  # n_samp x 231 global-mode h
    smap = {s: i for i, s in enumerate(z["samples"].astype(str))}
    nF = len(founders)

    pt = lib.pool_table()
    pt = pt[pt.sampleid.astype(str).isin(smap)].copy()
    clim = lib.load_climate()

    # timepoint samples -> pool (site_gen_plot), flower-weighted mean h
    pool_h, pool_meta = {}, {}
    for pool, g in pt.groupby("pool"):
        idx = [smap[str(s)] for s in g.sampleid]
        w = g.flowerscollected.to_numpy(float)
        w = np.where(np.isfinite(w) & (w > 0), w, 1.0); w = w / w.sum()
        pool_h[pool] = (Hs[idx] * w[:, None]).sum(0)
        r = g.iloc[0]
        pool_meta[pool] = (int(r["site"]), int(r["generation"]), int(r["plot"]))

    Tg_full = np.array([0.0, 1.0, 2.0, 3.0])
    sites, Srows, nplots, freq_last = [], [], [], []
    for site in sorted({m[0] for m in pool_meta.values()}):
        if site not in clim.index or not np.isfinite(clim.loc[site, "bio1"]):
            continue
        # per plot: trajectory over gens present (gen0 = p0)
        by_plot = {}
        for pool, (st, gen, plot) in pool_meta.items():
            if st == site and gen in TRAIT_GENS:
                by_plot.setdefault(plot, {})[gen] = pool_h[pool]
        plot_slopes, plot_last = [], []
        for plot, cells in by_plot.items():
            gens = sorted(cells)
            if 1 not in gens:                     # need the gen1 anchor
                continue
            t = np.array([0.0] + [float(gp) for gp in gens])
            Y = np.vstack([p0] + [cells[gp] for gp in gens])   # (T x nF) freq
            tc = t - t.mean()
            sl = (tc[:, None] * logit(Y)).sum(0) / (tc @ tc)   # per-founder logit slope
            plot_slopes.append(sl); plot_last.append(cells[max(gens)])
        if not plot_slopes:
            continue
        sites.append(site)
        Srows.append(np.mean(plot_slopes, 0))
        nplots.append(len(plot_slopes))
        freq_last.append(np.mean(plot_last, 0))

    sites = np.array(sites, int)
    S = np.vstack(Srows)                                        # n_site x nF
    bio1 = clim.loc[sites, "bio1"].to_numpy(float)
    analyzable = np.ones(nF, dtype=bool)   # keep ALL founders (p0 floored; none dropped)

    np.savez(f"{OUT}/selection_s_matrix.npz", S=S, sites=sites, founders=founders,
             bio1=bio1, p0=p0, analyzable=analyzable,
             n_plots=np.array(nplots), freq_last=np.vstack(freq_last))

    from scipy import stats
    sa = S[:, analyzable].ravel()
    print(f"sites={len(sites)}  founders total={nF}  analyzable(all, p0 floored@{P0_FLOOR:.0e})="
          f"{int(analyzable.sum())}  n_p0_floored={int((p0<=P0_FLOOR).sum())}")
    print(f"plots/site: min={min(nplots)} median={int(np.median(nplots))} max={max(nplots)}")
    print(f"bio1 {bio1.min():.1f}..{bio1.max():.1f} C")
    print(f"s (analyzable pooled): mean={sa.mean():+.3f} sd={sa.std():.3f} "
          f"skew={stats.skew(sa):+.2f} kurt={stats.kurtosis(sa):+.2f}")
    print(f"wrote {OUT}/selection_s_matrix.npz")


if __name__ == "__main__":
    main()

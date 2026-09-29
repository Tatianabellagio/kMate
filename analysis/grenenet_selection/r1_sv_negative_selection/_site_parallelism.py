#!/usr/bin/env python
"""Per-SITE parallelism of FOUNDER (ecotype) frequency change, to test whether the SV-purging
climate gradient (Section 1 climate scatter) is really a SELECTION-INTENSITY effect.

Hypothesis (user 2026-07-20): higher selection makes replicate plots more parallel, which inflates
the per-site selection coefficient, which is why we see SV purging — so the apparent climate gradient
could be a parallelism gradient in disguise.

Per site:
  * reuse build_selection_trait's mapping: sample_global_h.npz -> per (site,gen,plot) founder h
    (flower-weighted over timepoints); gen0 = seedmix p0; per plot, per founder logit-slope over gens.
  * PARALLELISM (same metric as _compute_parallelism.py, applied to founders): for founder f, over the
    n_plots per-plot slopes, rho_f = mean(slope)^2 / mean(slope^2) in [0,1] (1 = all plots agree,
    0 = cancel/drift). Per-site parallelism P = mean_f rho_f over the analyzable founders (p0>1e-3, so
    floored near-absent founders don't add 0/0 noise). Also a robust cross-check: mean pairwise
    Spearman corr of the per-plot founder-slope vectors, and the leading-eigenvalue fraction.
  * selection-intensity proxy: per-site mean |founder logit-slope| (how much founders move overall).

Then merge with the per-site SV purging (shift_sv from s_dist_by_stratum_sitemeta.csv) and climate,
and report: corr(purging, parallelism); corr(parallelism, bio1/bio18); and PARTIAL correlations
(purging~climate controlling parallelism, and purging~parallelism controlling climate) to separate
the two explanations.

Env: kmate. Writes analysis/grenenet_selection/r1_sv_negative_selection/results/sv_adaptive/site_parallelism.csv (+ prints stats).
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

os.chdir("/global/scratch/users/tbellg/kmate")
HCACHE = f"{lib.GEA}/r3_persite_gwas/results/ecotype_fitness/sample_global_h.npz"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
EPS, P0_FLOOR, TRAIT_GENS = 1e-3, 1e-5, (1, 2, 3)
REAL_P0 = 1e-3   # founders with seedmix p0 above this count toward the parallelism mean


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def seedmix_p0():
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


def partial_spearman(x, y, z):
    """Spearman corr of x,y after linearly regressing rank(z) out of both (rank-based partial)."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    rz1 = np.c_[np.ones_like(rz), rz]
    ex = rx - rz1 @ np.linalg.lstsq(rz1, rx, rcond=None)[0]
    ey = ry - rz1 @ np.linalg.lstsq(rz1, ry, rcond=None)[0]
    r, p = stats.pearsonr(ex, ey)
    return r, p


def main():
    founders, p0 = seedmix_p0()
    p0f = np.maximum(p0, P0_FLOOR)
    z = np.load(HCACHE, allow_pickle=True)
    assert (z["founders"].astype(str) == founders).all(), "founder order mismatch"
    Hs = z["H"].astype(float)
    smap = {s: i for i, s in enumerate(z["samples"].astype(str))}
    real = p0 > REAL_P0
    print(f"{len(founders)} founders; {real.sum()} with seedmix p0>{REAL_P0} used for parallelism mean")

    pt = lib.pool_table(); pt = pt[pt.sampleid.astype(str).isin(smap)].copy()
    clim = lib.load_climate()

    pool_h, pool_meta = {}, {}
    for pool, g in pt.groupby("pool"):
        idx = [smap[str(s)] for s in g.sampleid]
        w = g.flowerscollected.to_numpy(float)
        w = np.where(np.isfinite(w) & (w > 0), w, 1.0); w = w / w.sum()
        pool_h[pool] = (Hs[idx] * w[:, None]).sum(0)
        r = g.iloc[0]; pool_meta[pool] = (int(r["site"]), int(r["generation"]), int(r["plot"]))

    rows = []
    for site in sorted({m[0] for m in pool_meta.values()}):
        if site not in clim.index or not np.isfinite(clim.loc[site, "bio1"]):
            continue
        by_plot = {}
        for pool, (st, gen, plot) in pool_meta.items():
            if st == site and gen in TRAIT_GENS:
                by_plot.setdefault(plot, {})[gen] = pool_h[pool]
        slopes, effn = [], []
        for plot, cells in by_plot.items():
            gens = sorted(cells)
            if 1 not in gens:
                continue
            t = np.array([0.0] + [float(gp) for gp in gens]); tc = t - t.mean()
            Y = np.vstack([p0f] + [cells[gp] for gp in gens])
            slopes.append((tc[:, None] * logit(Y)).sum(0) / (tc @ tc))   # per-founder logit slope
            h_last = cells[max(gens)]                                    # evolved (last-gen) founder h
            effn.append(1.0 / np.sum(h_last ** 2))                       # eff. # surviving founders (inv-Simpson)
        if len(slopes) < 3:
            continue
        S = np.vstack(slopes)                      # n_plots x nF
        # per-founder parallelism rho_f = mean^2 / mean(sq); per-site P = mean over real founders
        m1 = S.mean(0); m2 = (S ** 2).mean(0)
        rho = np.divide(m1 ** 2, m2, out=np.zeros_like(m1), where=m2 > 0)
        P = float(rho[real].mean())
        # pairwise across-plot correlation of the founder-slope vectors (real founders only): the
        # simple, standard parallelism metric -- do plots agree on which founders rise/fall.
        Sr = S[:, real]
        pear = [stats.pearsonr(Sr[i], Sr[j])[0]
                for i in range(len(Sr)) for j in range(i + 1, len(Sr))]
        spear = [stats.spearmanr(Sr[i], Sr[j]).statistic
                 for i in range(len(Sr)) for j in range(i + 1, len(Sr))]
        par_pearson = float(np.mean(pear)); par_spear = float(np.mean(spear))
        C = np.corrcoef(Sr)                        # n_plots x n_plots
        eig = np.linalg.eigvalsh(C); eig1 = float(eig[-1] / eig.sum())
        intensity = float(np.abs(m1[real]).mean())  # selection-strength proxy
        diversity = float(np.mean(effn))            # evolved diversity: mean over plots of eff # founders
        rows.append(dict(site=site, bio1=round(float(clim.loc[site, "bio1"]), 2),
                         bio18=round(float(clim.loc[site, "bio18"]), 1),
                         n_plots=len(slopes), par_pearson=round(par_pearson, 4),
                         par_rho=round(P, 4), par_spearman=round(par_spear, 4),
                         par_eig1=round(eig1, 4), intensity=round(intensity, 4),
                         diversity=round(diversity, 2)))
    df = pd.DataFrame(rows)

    # merge per-site SV purging (shift_sv = mean over strata of SV median - SNP baseline)
    site_meta = pd.read_csv(f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive/s_dist_by_stratum_sitemeta.csv")
    df = df.merge(site_meta[["site", "shift_sv", "shift_ind"]], on="site", how="inner")
    df.to_csv(f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive/site_parallelism.csv", index=False)
    print(f"\n[wrote] site_parallelism.csv ({len(df)} sites)")

    def sp(a, b):
        r = stats.spearmanr(df[a], df[b]); return f"rho={r.statistic:+.3f} p={r.pvalue:.4f}"
    METR = ["par_pearson", "par_spearman", "par_rho", "par_eig1", "intensity", "diversity"]
    print("\n=== is SV purging (shift_sv; more negative = more purged) explained by these? ===")
    for par in METR:
        print(f"  shift_sv vs {par:12}: {sp('shift_sv', par)}")
    print("\n=== are these themselves climate-graded (would make climate a proxy)? ===")
    for par in METR:
        print(f"  {par:12} vs bio1 : {sp(par,'bio1')}   vs bio18: {sp(par,'bio18')}")
    print("\n=== SV purging vs climate, RAW and PARTIAL (controlling each selection proxy) ===")
    for ctrl in ["par_pearson", "intensity", "diversity"]:
        for cv in ["bio1", "bio18"]:
            raw = stats.spearmanr(df["shift_sv"], df[cv])
            pr, pp = partial_spearman(df["shift_sv"].to_numpy(), df[cv].to_numpy(), df[ctrl].to_numpy())
            print(f"  shift_sv~{cv} | {ctrl:11}: RAW rho={raw.statistic:+.3f} -> PARTIAL r={pr:+.3f} p={pp:.4f}")
    print("\nREAD: if purging tracks climate but NOT parallelism, and survives controlling for "
          "parallelism, the gradient is climatic, not a selection-intensity artifact. If purging "
          "tracks parallelism and climate~parallelism, and the climate partial dies, it's intensity.")


if __name__ == "__main__":
    main()

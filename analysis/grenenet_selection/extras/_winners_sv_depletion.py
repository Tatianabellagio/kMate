#!/usr/bin/env python
"""Q2 (user 2026-07-20): at the SV-purging sites, are the WINNING ecotypes SV-DEPLETED in the panel?
If the founders that rise at hot sites are the short-read (PanGenie) ecotypes, and short-read founders
carry fewer called SVs than the long-read (cactus) ones, then SV-alt alleles would decline as those
founders rise -> apparent "SV purging" that is a panel-completeness artifact, not selection.

Test:
  * founder SV-content: per founder, # common SV-alt alleles carried (var_pa, SV = |alt_len-ref_len|>50,
    same MAC 12-219 / call>=0.9 filter as the temporal analysis). Also indel/SNP counts for context.
  * founder modality: cactus (long-read, 80) vs PG (short-read, 151) from founder_split_cactus_pg.json.
  * per-site winning-ness: per-founder selection coefficient s (mean over plots of the logit-slope of
    founder h over gens 0..3), reusing build_selection_trait's mapping. s>0 = rising (winner).
  * at each purging site (and cold sites for contrast): Spearman(s, SV-count) across 231 founders;
    winners (top-decile s) vs losers (bottom-decile) mean SV-count and %cactus; plus the enabling
    facts (does SV-count track modality? do winners track modality?).

Interpretation: the artifact needs BOTH (a) winners are SV-depleted / PG-enriched AND (b) it is
specific to purging sites. If winners are not SV-depleted, or the pattern is the same at cold
(non-purging) sites, the SV purging is not a modality/panel-completeness artifact.

Env: kmate. Writes sv_adaptive/winners_sv_depletion.csv + founder_sv_content.csv (+ prints).
"""
import os, sys, glob, json
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

os.chdir("/global/scratch/users/tbellg/kmate")
HCACHE = f"{lib.GEA}/ecotype_fitness/sample_global_h.npz"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
EPS, P0_FLOOR, TRAIT_GENS = 1e-3, 1e-5, (1, 2, 3)
SV_BP, MIN_MAC, CALL_MIN = 50, 12, 0.9
PURGE_SITES = [4, 32, 43, 60, 26]


def logit(p):
    p = np.clip(p, EPS, 1 - EPS); return np.log(p / (1 - p))


def seedmix_p0():
    ids = sorted({os.path.basename(p).split("_Chr")[0]
                  for p in glob.glob(f"{lib.SEEDMIX}/*_Chr1.h_per_chrom.npz")})
    founders, reps = None, []
    for s in ids:
        chs = []
        for c in CHROMS:
            z = np.load(f"{lib.SEEDMIX}/{s}_{c}.h_per_chrom.npz", allow_pickle=True)
            if founders is None: founders = z["founders"].astype(str)
            chs.append(z[c].astype(float))
        reps.append(np.vstack(chs).mean(0))
    return founders, np.vstack(reps).mean(0)


def founder_class_counts(founders):
    """Per-founder count of common (MAC 12-219, call>=0.9) SV / indel / SNP alt alleles from var_pa."""
    F = len(founders); sv = np.zeros(F); ind = np.zeros(F); snp = np.zeros(F)
    for c in CHROMS:
        base = f"panel/arch3/{c.lower()}/var_pa_231_arch3_{c.lower()}"
        m = np.load(f"{base}.meta.npz", allow_pickle=True)
        assert (m["founders"].astype(str) == founders).all(), f"{c} founder order mismatch"
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
        vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        n_alt = np.asarray(vp.sum(0)).ravel(); n_cal = np.asarray(vc.sum(0)).ravel()
        keep = (n_alt >= MIN_MAC) & (n_alt <= F - MIN_MAC) & (n_cal >= CALL_MIN * F)
        dl = np.abs(m["alt_len"].astype(int) - m["ref_len"].astype(int))
        for arr, msk in ((sv, keep & (dl > SV_BP)),
                         (ind, keep & (dl >= 1) & (dl <= SV_BP)),
                         (snp, keep & (dl == 0))):
            cols = np.where(msk)[0]
            arr += np.asarray(vp[:, cols].sum(1)).ravel()
    return sv, ind, snp


def main():
    founders, p0 = seedmix_p0(); p0f = np.maximum(p0, P0_FLOOR)
    z = np.load(HCACHE, allow_pickle=True)
    assert (z["founders"].astype(str) == founders).all()
    Hs = z["H"].astype(float); smap = {s: i for i, s in enumerate(z["samples"].astype(str))}

    split = json.load(open("data/founder_split_cactus_pg.json"))
    cactus = set(map(str, split["cactus"]))
    is_cactus = np.array([f in cactus for f in founders])
    print(f"modality: {is_cactus.sum()} cactus (long-read), {(~is_cactus).sum()} PG (short-read)")

    sv, ind, snp = founder_class_counts(founders)
    print(f"\n=== founder SV-content by modality (common SVs, MAC {MIN_MAC}-{231-MIN_MAC}, call>={CALL_MIN}) ===")
    print(f"  SV-alt per founder:  cactus median={np.median(sv[is_cactus]):.0f}  "
          f"PG median={np.median(sv[~is_cactus]):.0f}  "
          f"MWU p={stats.mannwhitneyu(sv[is_cactus], sv[~is_cactus]).pvalue:.1e}")
    print(f"  (indel per founder:  cactus={np.median(ind[is_cactus]):.0f} PG={np.median(ind[~is_cactus]):.0f}; "
          f"SNP: cactus={np.median(snp[is_cactus]):.0f} PG={np.median(snp[~is_cactus]):.0f})")
    pd.DataFrame(dict(founder=founders, is_cactus=is_cactus, sv=sv, indel=ind, snp=snp, p0=p0)).to_csv(
        f"{lib.GEA}/sv_adaptive/founder_sv_content.csv", index=False)

    # per-site founder s (mean over plots of logit-slope over gens 0..3)
    pt = lib.pool_table(); pt = pt[pt.sampleid.astype(str).isin(smap)].copy()
    clim = lib.load_climate()
    pool_h, pool_meta = {}, {}
    for pool, g in pt.groupby("pool"):
        idx = [smap[str(s)] for s in g.sampleid]
        w = g.flowerscollected.to_numpy(float); w = np.where(np.isfinite(w) & (w > 0), w, 1.0); w /= w.sum()
        pool_h[pool] = (Hs[idx] * w[:, None]).sum(0)
        r = g.iloc[0]; pool_meta[pool] = (int(r["site"]), int(r["generation"]), int(r["plot"]))

    def site_s(site):
        by_plot = {}
        for pool, (st, gen, plot) in pool_meta.items():
            if st == site and gen in TRAIT_GENS: by_plot.setdefault(plot, {})[gen] = pool_h[pool]
        sl = []
        for plot, cells in by_plot.items():
            gens = sorted(cells)
            if 1 not in gens: continue
            t = np.array([0.0] + [float(gp) for gp in gens]); tc = t - t.mean()
            Y = np.vstack([p0f] + [cells[gp] for gp in gens])
            sl.append((tc[:, None] * logit(Y)).sum(0) / (tc @ tc))
        return np.mean(sl, 0) if sl else None

    all_sites = sorted({m[0] for m in pool_meta.values() if m[0] in clim.index})
    rows = []
    for site in all_sites:
        s = site_s(site)
        if s is None: continue
        # winners = top-decile s, losers = bottom-decile s
        hi = s >= np.quantile(s, 0.9); lo = s <= np.quantile(s, 0.1)
        rows.append(dict(
            site=site, bio1=round(float(clim.loc[site, "bio1"]), 1), purge=site in PURGE_SITES,
            rho_s_sv=round(stats.spearmanr(s, sv).statistic, 3),
            rho_s_svp=round(stats.spearmanr(s, sv).pvalue, 4),
            win_sv=round(float(sv[hi].mean()), 0), lose_sv=round(float(sv[lo].mean()), 0),
            allf_sv=round(float(sv.mean()), 0),
            win_pct_cactus=round(float(is_cactus[hi].mean()), 2),
            lose_pct_cactus=round(float(is_cactus[lo].mean()), 2),
            rho_s_cactus=round(stats.spearmanr(s, is_cactus.astype(float)).statistic, 3)))
    df = pd.DataFrame(rows)
    df.to_csv(f"{lib.GEA}/sv_adaptive/winners_sv_depletion.csv", index=False)

    print(f"\n=== does winning-ness (s) correlate with founder SV-content? (Spearman across 231 founders) ===")
    print("  negative rho_s_sv = winners are SV-depleted (the artifact direction)")
    print(df[["site", "bio1", "purge", "rho_s_sv", "rho_s_svp", "win_sv", "lose_sv",
              "win_pct_cactus", "lose_pct_cactus", "rho_s_cactus"]].to_string(index=False))
    pg = df[df.purge]; cold = df[~df.purge]
    print(f"\n=== PURGE sites {PURGE_SITES} vs the rest ===")
    print(f"  rho(s, SV-count):   purge mean={pg.rho_s_sv.mean():+.3f}   rest mean={cold.rho_s_sv.mean():+.3f}")
    print(f"  winners %cactus:    purge mean={pg.win_pct_cactus.mean():.2f}   rest mean={cold.win_pct_cactus.mean():.2f}  "
          f"(all-founder baseline={is_cactus.mean():.2f})")
    print(f"  rho(s, is_cactus):  purge mean={pg.rho_s_cactus.mean():+.3f}   rest mean={cold.rho_s_cactus.mean():+.3f}")
    print("\nREAD: artifact requires winners SV-depleted (rho_s_sv<0) AND PG-enriched (win_pct_cactus low, "
          "rho_s_cactus<0), specifically at purge sites. If purge and non-purge look alike, or winners "
          "aren't SV-poor/PG, the purging is not a short-read-vs-long-read panel artifact.")


if __name__ == "__main__":
    main()

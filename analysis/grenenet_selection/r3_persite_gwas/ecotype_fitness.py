"""Tier 1 — the ecotype-sorting map + per-founder fitness axes (the identifiable signal).

In global mode + ~97% selfing, per-variant selection is not identifiable (a SV and a
SNP on the same founders share one projected trajectory). What IS identifiable is the
per-founder frequency `h` — the ecotype frequencies. This script turns the per-pool `h`
vectors into the per-FOUNDER realized-fitness phenotypes that feed a 231-ecotype fitness
GWAS. Two FITNESS FLAVOURS x several AXES (user decision 2026-07-02):

  FITNESS FLAVOUR
   - relative : Δh from founding, plots equal-weighted within a site (what selection
                acts on within the pool — compositional).
   - census   : Δh from founding, plots weighted by their lifetime flowering census
                (cum_flowers, survival.csv) — a founder that ends up in demographically
                PRODUCTIVE plots scores higher. Demographic fitness. The relative-vs-census
                gap is itself a result (stress-selection decoupling, cf Price work).

  AXES (each computed for both flavours)
   - w_global : mean Δh across all sites (generalist success).
   - w_<zone> : mean Δh within a bio1 climate ZONE (cold/mid/hot tercile). The substrate
                for CONDITIONAL NEUTRALITY — a founder favoured in one zone, ~neutral in
                others (NOT the antagonistic win-hot/lose-cold trade-off).
   - c_climate: OLS slope of Δh on standardized bio1. Kept ONLY as the antagonistic-
                pleiotropy FOIL (the phase-1 monotonic signature we do NOT expect to
                dominate under conditional neutrality).

Per-founder CONDITIONAL-NEUTRALITY class from the zone means + a plot/site bootstrap:
  global (all zones same sign, sig) | antagonistic (hot & cold opposite sign, both sig) |
  cond_hot / cond_cold (sig in one zone, ~neutral in the other) | neutral.

Also a RELIABILITY proxy: cross-chromosome SD of h (identifiability wobble; founders whose
h disagrees across the 5 chromosomes are the non-identifiable ones — h-certainty work).

Outputs -> analysis/grenenet_selection/r3_persite_gwas/results/ecotype_fitness/
  ecotype_fitness.csv     (founder, h0, + w_/c_ for both flavours, cn_class, xchrom_sd)
  founder_site_dh.npz     (founders, sites, bio1, zone, DH_rel, DH_cen, H0)
  sample_global_h.npz     (cache: samples, H[n_samp x 231], SD, founders) — built once
Light once cached (231 x 745 pools); runs in the `kmate` env.
"""
from __future__ import annotations
import os, glob, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import lib

OUT = f"{lib.GEA}/r3_persite_gwas/results/ecotype_fitness"
CACHE = f"{OUT}/sample_global_h.npz"
CHROMS = [f"Chr{i}" for i in range(1, 6)]
SURV = "/global/home/users/tbellg/scratch/grene/data/survival.csv"
GENS_CENSUS = (1, 2, 3, 4, 5)


def load_census():
    """Per-plot lifetime flowering census (cum_flowers) — the demographic weight.

    Mirrors build_fitness_table.load_fitness: sum flowerstotal over gens with survival in
    {0,1} (death=real 0; drop nan/-1). Returns dict (site,plot) -> cum_flowers.
    """
    sv = pd.read_csv(SURV)
    sv = sv[pd.to_numeric(sv["site"], errors="coerce").notna()
            & pd.to_numeric(sv["plot"], errors="coerce").notna()].copy()
    sv["site"] = sv["site"].astype(float).astype(int)
    sv["plot"] = sv["plot"].astype(float).astype(int)
    cum = np.zeros(len(sv)); ngen = np.zeros(len(sv))
    for g in GENS_CENSUS:
        s = pd.to_numeric(sv.get(f"{g}_survival"), errors="coerce").to_numpy()
        f = pd.to_numeric(sv.get(f"{g}_flowerstotal"), errors="coerce").to_numpy()
        use = np.isin(s, [0, 1]) & np.isfinite(f)
        cum[use] += f[use]; ngen[use] += 1
    return {(int(s), int(p)): c for s, p, c, n in
            zip(sv["site"], sv["plot"], cum, ngen) if n > 0}


def load_h_matrix(sample_ids, base):
    """founder x sample h, averaged over the 5 chromosomes; + cross-chrom SD per cell.

    Each *_Chr{N}.h_per_chrom.npz holds that chrom's 231-vector h (sums to 1). We average
    across chroms for a genome-level founder frequency, and keep the across-chrom SD as the
    per-(founder,sample) identifiability wobble.
    """
    founders = None
    H = {}      # sample -> (231,) mean-over-chrom h
    SD = {}     # sample -> (231,) sd-over-chrom h
    dropped = []
    for s in sample_ids:
        chr_hs = []
        ok = True
        for c in CHROMS:
            p = f"{base}/{s}_{c}.h_per_chrom.npz"
            if not os.path.exists(p):
                ok = False
                break
            z = np.load(p, allow_pickle=True)
            if founders is None:
                founders = z["founders"].astype(str)
            chr_hs.append(z[c].astype(float))
        if not ok or not chr_hs:
            continue
        M = np.vstack(chr_hs)               # 5 x 231
        # Drop degenerate samples: near-zero-coverage / failed libraries whose EM has no
        # data (all counts removed by the repeat-guard) produce all-NaN h. Including them
        # NaN-poisons every pool/site they touch. NaN is the honest "no data" output; the
        # sample is excluded from all downstream (via the cache -> pt.isin(smap) filter).
        if not np.isfinite(M).all():
            dropped.append(s)
            continue
        H[s] = M.mean(0)
        SD[s] = M.std(0)
    if dropped:
        print(f"  load_h_matrix: dropped {len(dropped)} degenerate (non-finite h) samples: "
              f"{', '.join(dropped)}", flush=True)
    return founders, H, SD


def build_pool_h(pt):
    """founder x sample global-h matrix (cached) + cross-chrom SD; aligned to pt rows.

    Cached to CACHE so the ~10k-npz Lustre pass runs once. Returns (founders, Hmat, SDmat)
    with rows aligned to the (filtered) pt.sampleid order; pt is modified in place to the
    samples that have h on disk.
    """
    if os.path.exists(CACHE):
        z = np.load(CACHE, allow_pickle=True)
        founders = z["founders"].astype(str)
        smap = {s: i for i, s in enumerate(z["samples"].astype(str))}
        keep = [s for s in pt.sampleid if s in smap]
        idx = [smap[s] for s in keep]
        return founders, z["H"][idx], z["SD"][idx], keep
    f_evo, He, SDe = load_h_matrix(list(pt.sampleid), lib.OUT)
    keep = [s for s in pt.sampleid if s in He]
    Hmat = np.vstack([He[s] for s in keep])
    SDmat = np.vstack([SDe[s] for s in keep])
    np.savez(CACHE, samples=np.array(keep), H=Hmat, SD=SDmat, founders=f_evo)
    return f_evo, Hmat, SDmat, keep


def _site_comp(pool_h, pm, clim, weight):
    """Per-site last-gen founder composition, plots weighted by `weight` (dict pool->w).

    Returns (sites, H_last[n_site x 231]). weight=None => equal plot weight (relative).
    """
    sites, H_last = [], []
    for site, g in pm.groupby("site"):
        if site not in clim.index:
            continue
        lg = g.generation.max()
        pools = list(g[g.generation == lg].pool)
        W = np.array([1.0 if weight is None else max(weight.get(p, 0.0), 0.0) for p in pools])
        if W.sum() <= 0:
            W = np.ones(len(pools))
        W = W / W.sum()
        H_last.append(sum(w * pool_h[p] for w, p in zip(W, pools)))
        sites.append(site)
    return np.array(sites), np.vstack(H_last)


def _axes(DH, bio1, zone, tag):
    """Per-founder fitness axes from a Δh[site x founder] matrix. Returns a dict of columns."""
    x = (bio1 - bio1.mean()) / bio1.std()
    out = {f"w_global_{tag}": DH.mean(0),
           f"c_climate_{tag}": (DH * x[:, None]).sum(0) / (x @ x)}
    for z in ("cold", "mid", "hot"):
        m = zone == z
        out[f"w_{z}_{tag}"] = DH[m].mean(0) if m.any() else np.full(DH.shape[1], np.nan)
    return out


def _cn_class(DH, zone, seed=0, nboot=2000):
    """Conditional-neutrality class per founder from a site x founder Δh matrix.

    Bootstrap sites within each zone; a zone effect is 'sig' if its 90% CI excludes 0.
    Classes: global (cold&hot same sign, both sig) | antagonistic (opposite sign, both sig)
    | cond_hot / cond_cold (one sig, other not) | neutral.
    """
    rng = np.random.default_rng(seed)
    zi = {z: np.where(zone == z)[0] for z in ("cold", "hot")}
    est, lo, hi = {}, {}, {}
    for z, ix in zi.items():
        if len(ix) == 0:
            est[z] = lo[z] = hi[z] = np.zeros(DH.shape[1]); continue
        boot = np.stack([DH[rng.choice(ix, len(ix))].mean(0) for _ in range(nboot)])
        est[z] = DH[ix].mean(0)
        lo[z], hi[z] = np.percentile(boot, [5, 95], axis=0)
    sig = {z: (lo[z] > 0) | (hi[z] < 0) for z in ("cold", "hot")}
    cls = np.full(DH.shape[1], "neutral", dtype=object)
    for i in range(DH.shape[1]):
        sc, sh = sig["cold"][i], sig["hot"][i]
        if sc and sh:
            cls[i] = "global" if np.sign(est["cold"][i]) == np.sign(est["hot"][i]) else "antagonistic"
        elif sh:
            cls[i] = "cond_hot"
        elif sc:
            cls[i] = "cond_cold"
    return cls


def main():
    os.makedirs(OUT, exist_ok=True)
    pt = lib.pool_table()                         # sample -> pool/site/plot/gen/flowers
    clim = lib.load_climate()["bio1"]
    census = load_census()                        # (site,plot) -> cum_flowers

    # --- founding h0: mean over the 8 SEEDMIX reps -----------------------------------
    seed_ids = sorted({os.path.basename(p).split("_Chr")[0]
                       for p in glob.glob(f"{lib.SEEDMIX}/*_Chr1.h_per_chrom.npz")})
    founders, Hs, _ = load_h_matrix(seed_ids, lib.SEEDMIX)
    h0 = np.mean([Hs[s] for s in Hs], axis=0)     # (231,)

    # --- evolved pool h (cached) -----------------------------------------------------
    f_evo, Hmat, SDmat, keep = build_pool_h(pt)
    assert (founders == f_evo).all(), "founder order mismatch seedmix vs evolved"
    pt = pt.set_index("sampleid").loc[keep].reset_index()

    # collapse timepoint samples -> pool (site_gen_plot), flower-weighted mean h
    w = pt.flowerscollected.to_numpy(float)
    w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
    pool_h, pool_meta, pool_plot = {}, {}, {}
    for pool, idx in pd.Series(range(len(pt))).groupby(pt.pool.to_numpy()).groups.items():
        idx = np.asarray(idx)
        ww = w[idx] / w[idx].sum()
        pool_h[pool] = (Hmat[idx] * ww[:, None]).sum(0)
        r = pt.iloc[idx[0]]
        pool_meta[pool] = (int(r["site"]), int(r["generation"]))
        pool_plot[pool] = (int(r["site"]), int(r["plot"]))
    pm = pd.DataFrame([(k, *v) for k, v in pool_meta.items()],
                      columns=["pool", "site", "generation"])
    pool_census = {p: census.get(pool_plot[p], 0.0) for p in pool_h}

    # --- per-site composition: relative (equal plots) and census (cum_flowers wt) -----
    sites, H_rel = _site_comp(pool_h, pm, clim, weight=None)
    _, H_cen = _site_comp(pool_h, pm, clim, weight=pool_census)
    bio1 = clim.loc[sites].to_numpy(float)
    zone = pd.qcut(bio1, 3, labels=["cold", "mid", "hot"]).astype(str)  # bio1 terciles
    DH_rel = H_rel - h0[None, :]
    DH_cen = H_cen - h0[None, :]

    # --- assemble per-founder table --------------------------------------------------
    cols = {"founder": founders, "h0": h0, "hbar_last": H_rel.mean(0),
            "xchrom_sd": SDmat.mean(0), "n_sites": len(sites)}
    cols.update(_axes(DH_rel, bio1, zone, "rel"))
    cols.update(_axes(DH_cen, bio1, zone, "cen"))
    cols["cn_class_rel"] = _cn_class(DH_rel, zone)
    cols["cn_class_cen"] = _cn_class(DH_cen, zone)
    df = pd.DataFrame(cols).sort_values("w_global_rel", ascending=False)
    df.to_csv(f"{OUT}/ecotype_fitness.csv", index=False)
    np.savez(f"{OUT}/founder_site_dh.npz", founders=founders, sites=sites, bio1=bio1,
             zone=np.array(zone), DH_rel=DH_rel, DH_cen=DH_cen, H0=h0)

    # --- summary ---------------------------------------------------------------------
    print(f"founders={len(founders)}  sites={len(sites)}  bio1 {bio1.min():.1f}..{bio1.max():.1f}")
    zc = pd.Series(zone).value_counts()
    print(f"zones (sites): cold={zc.get('cold',0)} mid={zc.get('mid',0)} hot={zc.get('hot',0)}")
    eff0 = 1/np.sum(h0**2); effL = 1/np.sum(H_rel.mean(0)**2)
    print(f"eff_n founders  founding {eff0:.1f} -> evolved(last,mean) {effL:.1f}")
    print("\nConditional-neutrality class (relative flavour):")
    print(pd.Series(df.cn_class_rel).value_counts().to_string())
    print("\nrelative vs census GLOBAL fitness corr: %.3f" %
          np.corrcoef(df.w_global_rel, df.w_global_cen)[0, 1])
    print("cold-zone vs hot-zone relative fitness corr: %.3f (near 0/neg => antagonistic; "
          ">0 => shared winners)" % np.corrcoef(df.w_cold_rel, df.w_hot_rel)[0, 1])
    show = ["founder", "h0", "w_cold_rel", "w_hot_rel", "cn_class_rel", "xchrom_sd"]
    print("\nTop 8 by cold-zone fitness:")
    print(df.sort_values("w_cold_rel").tail(8)[show].to_string(index=False))
    print("\nTop 8 by hot-zone fitness:")
    print(df.sort_values("w_hot_rel").tail(8)[show].to_string(index=False))
    print(f"\nwrote {OUT}/ecotype_fitness.csv  + founder_site_dh.npz  + {os.path.basename(CACHE)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Two-stage SV climate-GEA with an EMPIRICAL drift null (the principled fix).

The naive Kendall (GIF~3) and LFMM (GIF 2-3) inflate because climate is a
SITE-level (Level-2) predictor: the honest N is #sites (20), not #pools (193).
The MixedLM `~bio1+(1|site)` swung to null (GIF~0) because the site random
intercept IS the site mean = what climate explains. The fix is to make the SITE
the unit and MEASURE the drift null from the replicate plots:

  Stage 1 (per SV, per unit): a per-plot/-lineage statistic y that differences
    out the shared SEEDMIX founding p0:
      --stat dp     y = p3_plot - p0                    (gen-3 endpoint Δp)
      --stat scoef  y = logit-slope of [p0,p1,p2,p3]    (selection coefficient,
                        per lineage present in all gens; t=0,1,2,3)
  Stage 1.5: collapse plots within a site -> site mean ȳ_s; the POOLED
    within-site variance V_w (same climate + same p0, so plot spread = drift +
    pool-sampling noise) is the empirical noise floor. df = N_units - N_sites.
  Stage 2: weighted meta-regression of ȳ_s on standardized climate across the 20
    sites, weights w_s = n_plots_s (inverse-variance under homoscedastic drift).
    β1 = climate effect. Calibrated by a SITE-LEVEL PERMUTATION null (permute the
    climate label across sites) — the only valid null for a Level-2 predictor.

Outputs (--out dir), per stat:
  twostage_{stat}_{climate}.npz  beta, se, p_model, z_emp, p_perm, V_within,
      n_sites_used, gif_model, gif_emp + meta (chrom,pos,ref_len,alt_len,sv_size,
      p0, site_dp_min/max) for the SV set (n_called>=NC_MIN & |Δlen|>SV_MIN_BP)
  twostage_{stat}_{climate}.top.csv  top hits by |z_emp|, signed (up-in-warm +)

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY analysis/grenenet_gea/build_two_stage_gea.py --stat dp scoef --climate bio1
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

PM = f"{lib.GEA}/pool_matrices"
STORE = lib.AF_STORE
CACHE = f"{STORE}/sv_support_cache.npz"
NC_MIN = 150            # SV support filter (see notebooks/06_sv_support_filter)
SV_MIN_BP = 50
EPS = 1e-3              # logit clip for the selection coefficient
N_PERM = 2000


def sv_mask_and_meta():
    """Boolean SV mask over the 2.25M non-SNP records + per-SV meta (n_called≥150)."""
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    rl = idx["ref_len"].astype(np.int64); al = idx["alt_len"].astype(np.int64)
    size = np.abs(al - rl)
    nc = np.asarray(np.load(sorted(__import__("glob").glob(f"{STORE}/nc_nonsnp/*.npy"))[0]))
    mask = (size > SV_MIN_BP) & (nc >= NC_MIN)
    meta = pd.DataFrame(dict(chrom=idx["chrom"].astype("U5")[mask], pos=idx["pos"][mask],
                             ref_len=rl[mask], alt_len=al[mask], sv_size=size[mask],
                             n_called=nc[mask]))
    return mask, meta


def load_gen(g, sv_idx):
    """SV-subset gen-g pool AF matrix [n_pools x nSV] + its meta (site/plot/bio).

    Read the whole matrix into RAM (sequential, seconds) then subset columns
    in-memory: fancy-indexing 100k+ scattered columns straight off the memmap is
    I/O-pathological (strided gathers across a 1.7-2.9 GB file).
    """
    P = np.load(f"{PM}/pool_gen{g}_nonsnp_af.npy")[:, sv_idx].astype(np.float32)
    m = pd.read_csv(f"{PM}/pool_gen{g}_nonsnp.meta.csv")
    return P, m


def stage1_dp(sv_idx, p0):
    """Endpoint Δp per gen-3 plot. Returns y [n_units x nSV], site[n_units]."""
    P3, m3 = load_gen(3, sv_idx)
    y = P3 - p0[None, :]
    return y, m3["site"].to_numpy(), m3, "gen3-plot"


def stage1_scoef(sv_idx, p0):
    """Selection coefficient = logit-slope over [p0,p1,p2,p3] per lineage present
    in all gens. Returns s [n_lineages x nSV], site[n_lineages]."""
    gens = {g: load_gen(g, sv_idx) for g in (1, 2, 3)}
    # lineage = site_plot present in g1 & g2 & g3
    keys = {g: (gens[g][1].site.astype(str) + "_" + gens[g][1]["plot"].astype(str))
            for g in (1, 2, 3)}
    common = sorted(set(keys[1]) & set(keys[2]) & set(keys[3]))
    rows = {g: {k: i for i, k in enumerate(keys[g])} for g in (1, 2, 3)}
    nL = len(common); nSV = len(sv_idx)
    P = {g: gens[g][0] for g in (1, 2, 3)}
    # stack freqs [4 x nL x nSV]: t=0 (p0, shared), 1,2,3
    F = np.empty((4, nL, nSV), dtype=np.float32)
    F[0] = p0[None, :]
    site = np.empty(nL, dtype=int)
    m3 = gens[3][1]; site3 = {k: s for k, s in zip(keys[3], m3.site)}
    for i, k in enumerate(common):
        F[1, i] = P[1][rows[1][k]]; F[2, i] = P[2][rows[2][k]]; F[3, i] = P[3][rows[3][k]]
        site[i] = site3[k]
    L = np.log(np.clip(F, EPS, 1 - EPS) / (1 - np.clip(F, EPS, 1 - EPS)))   # logit
    t = np.array([0, 1, 2, 3.0]); tc = t - t.mean()                         # Σtc²=5
    s = np.tensordot(tc, L, axes=(0, 0)) / np.sum(tc ** 2)                   # [nL x nSV]
    m_lin = pd.DataFrame(dict(site=site))
    return s, site, m_lin, f"{nL}-lineage-trajectory"


def stage2(y, site, x_site_df, climate, n_perm, rng_seed=0):
    """Empirical-drift two-stage meta-regression of site-mean(y) on climate.

    y [n_units x nSV], site [n_units]; x_site_df indexed by site -> climate value.
    Weights w_s = n_plots_s. Calibrated by a site-permutation null.
    """
    sites = np.array(sorted(pd.unique(site)))
    S = len(sites); nSV = y.shape[1]
    x = x_site_df.loc[sites, climate].to_numpy(float)
    x = (x - x.mean()) / x.std()                       # standardized climate

    # site means + per-site counts (NaN-aware over plots/lineages)
    site_mean = np.full((S, nSV), np.nan, np.float64)
    n_s = np.zeros(S)
    ss_within = np.zeros(nSV); df_within = 0
    for k, s in enumerate(sites):
        rows = y[site == s]
        cnt = np.sum(np.isfinite(rows), axis=0)
        mean = np.nanmean(rows, axis=0)
        site_mean[k] = mean
        n_s[k] = rows.shape[0]
        resid = rows - mean[None, :]
        ss_within += np.nansum(resid ** 2, axis=0)
        df_within += np.sum(np.isfinite(rows[:, 0])) - 1   # df per SV ≈ units-sites
    V_w = ss_within / max(df_within, 1)                    # pooled within-site var

    # Stage-2 weighted least squares (weights w_s = n_s, shared across SVs)
    w = n_s
    xw = np.sum(w * x) / np.sum(w)
    dx = x - xw                                            # centered climate
    D = np.sum(w * dx ** 2)                                # Σ w (x-x̄)²
    u = (w * dx) / D                                       # β1 = u · site_mean
    beta = u @ site_mean                                   # [nSV]
    b0 = (np.sum(w[:, None] * site_mean, 0) / np.sum(w)) - beta * xw
    fitted = b0[None, :] + beta[None, :] * x[:, None]
    resid = site_mean - fitted
    mse = np.sum(w[:, None] * resid ** 2, 0) / max(S - 2, 1)
    se = np.sqrt(mse / D)
    with np.errstate(invalid="ignore", divide="ignore"):
        from scipy.stats import t as tdist
        tstat = beta / se
        p_model = 2 * tdist.sf(np.abs(tstat), df=max(S - 2, 1))

    # site-permutation null: permute x across sites, recompute β1 for all SVs
    rng = np.random.default_rng(rng_seed)
    U = np.empty((n_perm, S))
    for b in range(n_perm):
        xp = x[rng.permutation(S)]
        xwp = np.sum(w * xp) / np.sum(w); dxp = xp - xwp
        U[b] = (w * dxp) / np.sum(w * dxp ** 2)
    beta_perm = U @ site_mean                              # [n_perm x nSV]
    sd_perm = beta_perm.std(0); sd_perm[sd_perm == 0] = np.nan
    z_emp = beta / sd_perm
    # genome-wide pooled null of the standardized statistic
    z_null = (beta_perm / sd_perm[None, :]).ravel()
    z_null = z_null[np.isfinite(z_null)]
    nN = len(z_null)
    sorted_null = np.sort(z_null)
    az = np.abs(z_emp)
    order = np.searchsorted(np.sort(np.abs(z_null)), az)
    p_perm = np.clip(1.0 - order / nN, 1.0 / nN, 1.0)       # two-sided
    p_perm[~np.isfinite(z_emp)] = np.nan
    # one-tailed UP-in-warm p: fraction of null z >= observed signed z_emp
    p_perm_up = np.clip((nN - np.searchsorted(sorted_null, z_emp)) / nN, 1.0 / nN, 1.0)
    p_perm_up[~np.isfinite(z_emp)] = np.nan

    from scipy.stats import chi2, norm
    def gif(p):
        ok = np.isfinite(p) & (p > 0)
        return float(np.median(norm.isf(p[ok] / 2) ** 2) / chi2.ppf(0.5, 1))
    return dict(beta=beta, se=se, p_model=p_model, z_emp=z_emp, p_perm=p_perm,
                p_perm_up=p_perm_up, V_within=V_w, site_mean=site_mean, sites=sites,
                x=x, n_s=n_s, gif_model=gif(p_model), gif_emp=gif(p_perm), n_sites=S)


def run(stat, climate, mask, sv_meta, sv_idx, p0, clim_by_site, n_perm, outdir):
    f = {"dp": stage1_dp, "scoef": stage1_scoef}[stat]
    y, site, _, unit = f(sv_idx, p0)
    print(f"[{stat}] Stage-1 unit = {unit}: {y.shape[0]} units x {y.shape[1]:,} SVs", flush=True)
    R = stage2(y, site, clim_by_site, climate, n_perm)
    sdp = R["site_mean"]
    out = dict(beta=R["beta"].astype(np.float32), se=R["se"].astype(np.float32),
               p_model=R["p_model"].astype(np.float32), z_emp=R["z_emp"].astype(np.float32),
               p_perm=R["p_perm"].astype(np.float32),
               p_perm_up=R["p_perm_up"].astype(np.float32),
               V_within=R["V_within"].astype(np.float32),
               site_dp_min=np.nanmin(sdp, 0).astype(np.float32),
               site_dp_max=np.nanmax(sdp, 0).astype(np.float32),
               chrom=sv_meta.chrom.to_numpy().astype("U5"), pos=sv_meta.pos.to_numpy(),
               ref_len=sv_meta.ref_len.to_numpy(), alt_len=sv_meta.alt_len.to_numpy(),
               sv_size=sv_meta.sv_size.to_numpy(), p0=p0.astype(np.float32),
               n_sites=R["n_sites"], gif_model=R["gif_model"], gif_emp=R["gif_emp"])
    path = f"{outdir}/twostage_{stat}_{climate}.npz"
    np.savez(path, **out)
    # top table (true SVs by |z_emp|, signed)
    d = pd.DataFrame(dict(chrom=out["chrom"], pos=out["pos"], ref_len=out["ref_len"],
                          alt_len=out["alt_len"], sv_size=out["sv_size"], p0=out["p0"],
                          beta=out["beta"], z_emp=out["z_emp"], p_perm=out["p_perm"]))
    d = d[np.isfinite(d.z_emp)].reindex(d.z_emp.abs().sort_values(ascending=False).index)
    d.head(200).to_csv(f"{outdir}/twostage_{stat}_{climate}.top.csv", index=False)
    bonf = 0.05 / np.isfinite(out["z_emp"]).sum()
    print(f"[{stat}] GIF model={R['gif_model']:.2f}  emp(perm)={R['gif_emp']:.2f}  "
          f"(naive Kendall ~3, LFMM ~2-3, mixedLM ~0)")
    print(f"[{stat}] p_perm<0.05: {(out['p_perm']<0.05).sum():,} | "
          f"<{bonf:.1e}(Bonf): {(out['p_perm']<bonf).sum()} | min p_perm {np.nanmin(out['p_perm']):.1e}")
    print(f"[{stat}] -> {path}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stat", nargs="+", default=["dp", "scoef"], choices=["dp", "scoef"])
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    ap.add_argument("--out", default=f"{lib.GEA}/gea")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    mask, sv_meta = sv_mask_and_meta()
    sv_idx = np.where(mask)[0]
    p0 = np.load(f"{STORE}/p0_nonsnp.npy")[mask].astype(np.float64)
    # climate per site (constant within site) from the gen-3 meta
    m3 = pd.read_csv(f"{PM}/pool_gen3_nonsnp.meta.csv")
    clim_by_site = m3.groupby("site")[lib.BIO_COLS].first()
    print(f"SVs (|Δlen|>{SV_MIN_BP} & n_called>={NC_MIN}): {mask.sum():,} | "
          f"climate={args.climate}, sites={clim_by_site.shape[0]}, perms={args.n_perm}\n")
    for stat in args.stat:
        run(stat, args.climate, mask, sv_meta, sv_idx, p0, clim_by_site, args.n_perm, args.out)


if __name__ == "__main__":
    main()

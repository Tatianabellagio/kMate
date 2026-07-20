#!/usr/bin/env python
"""TEMPORAL two-stage GEA — Stage 2 + Stage-1 selection summary.

Reads temporal_stage1.py output (per-site per-variant temporal selection coef s
and among-plot drift SE) and asks the two questions in the plan:

STAGE 2 (climate GEA):  the STRUCTURE-CORRECTED, baseline-comparable test — the
  exact single-timepoint permnull_site_gea.py machinery, but on the TEMPORAL
  Stage-1 selection-coefficient matrix s [31 x nv] instead of the single-timepoint
  Delta-p. Per variant we have 31 per-site temporal selection coefficients; the
  question is whether that trajectory-derived vector tracks climate better than the
  single Delta-p snapshot did (which was NULL under the same null).
    1. Yc = s centered per variant over sites; SVD over sites -> remove top-K
       site-PCs (latent-structure / shared-drift axes, LFMM analog); unit-norm each
       variant column. Per-locus statistic = cosine(Yr_locus, x_std) = partial
       correlation controlling for the K site-PCs, in [-1,1].
    2. SITE-LABEL PERMUTATION null (exchangeable unit = the 31 SITES; x permuted
       >= NPERM):  (a) genome-wide MAX |cos| over the tested family -> 5% = FWER thr;
       (b) per-locus empirical p = (1 + #{|perm| >= |obs|}) / (1 + NPERM);
       BH-FDR on emp p.  Run for K in {1,3}, axes {bio5,pc1}.
  Also reports (secondary, structure-UNCORRECTED, for direction only) the raw IV
  slope beta_iv = weighted regression of s on x, w_i = 1/se_i^2 — NOT a significance
  stat (single-site domination + shared among-site structure make its raw FDR
  invalid). site-MAF >= 0.05 (Stage-1 flower-weighted site-mean AF), >= MIN_SITES
  contributing sites. clq0.9 blocks; annotate hit genes.

STAGE-1-ONLY (selection, no climate):  per site, per variant t = s/se on
  df = n_plots-1 -> two-sided p; BH-FDR within (site). Counts variants/blocks
  under convincing temporal selection at ANY site and RECURRENTLY across sites,
  SNP vs non-SNP, with gene clustering. This has power even where the n=31
  climate step does not.

Outputs under analysis/grenenet_gea/gea_newpanel/results/temporal_twostage/.
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY temporal_stage2.py --class nonsnp
"""
from __future__ import annotations
import argparse, os, sys, json, time
import numpy as np
import pandas as pd
from scipy.stats import t as tdist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))            # analysis/grenenet_gea
sys.path.insert(0, HERE)                             # gea_newpanel
import lib
from blocks_clq09 import assign_clq09_blocks

GEA = lib.GEA
TW = f"{GEA}/gea_newpanel/results/temporal_twostage"
ENVD = f"{GEA}/gea_newpanel/results/env_site"
NSITES = 31
MAF_MIN = 0.05
MIN_SITES = 20                # min # of sites contributing a finite (s,se) to a locus
SE_FLOOR = 1e-6              # avoid inf weights when se==0
KS = [1, 3]                  # site-PCs removed (matches single-timepoint permnull baseline)


def bh_fdr(p):
    p = np.asarray(p, float); n = p.size
    o = np.argsort(p); r = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(q, 0, 1)
    return out


def std_env(axis):
    x = pd.read_csv(f"{ENVD}/env_site_{axis}.csv").iloc[:, 0].to_numpy(float)
    assert x.size == NSITES, (axis, x.size)
    x = x - x.mean(); x = x / np.sqrt((x * x).sum())   # unit-norm (permnull convention)
    return x


def residualize(Y, K):
    """permnull_site_gea.residualize: center per variant over sites, remove top-K
    site-PCs (SVD), unit-norm each variant column. Y = [nsites x nv]."""
    Yc = Y - Y.mean(axis=0, keepdims=True)
    U, _, _ = np.linalg.svd(Yc, full_matrices=False)   # U: [nsites x nsites]
    Uk = U[:, :K]
    Yr = Yc - Uk @ (Uk.T @ Yc)
    nrm = np.sqrt((Yr * Yr).sum(axis=0)); nrm[nrm == 0] = np.inf
    return Yr / nrm[None, :]


def stage1_selection(cls, z, meta_df, genes, axes_out):
    """Stage-1-only per-site drift-null selection scan. z carries s[31xnv],
    se[31xnv], n_plots[31]. Returns per-variant recurrence + writes a summary."""
    s = z["s"]; se = z["se"]; nplots = z["n_plots"]
    sites = z["sites"]; nS, nv = s.shape
    maf_ok = meta_df["site_maf"].to_numpy() >= MAF_MIN
    # per-site t-stat / two-sided p (df = n_plots-1); only sites with >=2 plots
    p05 = np.zeros(nv, np.int32)          # # sites with raw p<0.05
    fdr_sig = np.zeros(nv, np.int32)      # # sites BH-q<0.05 (within site, MAF loci)
    calib = []
    for i in range(nS):
        df = int(nplots[i]) - 1
        if df < 1:
            continue
        with np.errstate(invalid="ignore", divide="ignore"):
            t = s[i] / se[i]
            pv = 2.0 * tdist.sf(np.abs(t), df)
        fin = np.isfinite(pv) & maf_ok
        p05 += (fin & (pv < 0.05)).astype(np.int32)
        # BH within this site's MAF-passing loci
        q = np.ones(nv)
        if fin.sum() > 0:
            q[fin] = bh_fdr(pv[fin])
        fdr_sig += (fin & (q < 0.05)).astype(np.int32)
        calib.append((int(sites[i]), int(fin.sum()), float((pv[fin] < 0.05).mean())
                      if fin.sum() else np.nan, int((fin & (q < 0.05)).sum())))
    cal = pd.DataFrame(calib, columns=["site", "n_maf_loci", "frac_p<.05", "n_fdr_sig"])
    cal.to_csv(f"{TW}/stage1_selection_calib_{cls}.csv", index=False)
    return p05, fdr_sig, cal


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True, choices=["snp", "nonsnp"])
    ap.add_argument("--axes", default="bio5,pc1")
    ap.add_argument("--n-perm", type=int, default=2000)
    args = ap.parse_args()
    axes = args.axes.split(",")
    cls = args.cls
    os.makedirs(TW, exist_ok=True)
    genes = lib.load_genes()

    t0 = time.time()
    z = np.load(f"{TW}/stage1_{cls}.npz", allow_pickle=False)
    chrom = z["chrom"].astype("U5"); pos = z["pos"]
    rl = z["ref_len"].astype(np.int64); al = z["alt_len"].astype(np.int64)
    size = np.abs(al - rl)
    s = z["s"].astype(np.float64); se = z["se"].astype(np.float64)   # [31 x nv]
    sf = z["site_freq"].astype(np.float64)                          # [31 x nv]
    sites = z["sites"]; nS, nv = s.shape
    print(f"[{cls}] loaded stage1 {s.shape} ({time.time()-t0:.0f}s)", flush=True)

    # site-MAF = folded mean over sites of flower-weighted site AF
    with np.errstate(invalid="ignore"):
        mean_af = np.nanmean(sf, axis=0)
    site_maf = np.minimum(mean_af, 1.0 - mean_af)

    # blocks (once)
    blocks = assign_clq09_blocks(chrom, pos)
    meta_df = pd.DataFrame(dict(chrom=chrom, pos=pos, ref_len=rl, alt_len=al,
                                sv_size=size, p0=z["p0"], site_maf=site_maf,
                                block_clq09=blocks))

    # ---- tested family: MAF>=.05 & >=MIN_SITES contributing sites ----
    good = np.isfinite(s) & np.isfinite(se) & (se > 0)
    n_eff_sites = good.sum(0)
    tested = (site_maf >= MAF_MIN) & (n_eff_sites >= MIN_SITES)
    idx = np.where(tested)[0]
    print(f"[{cls}] tested loci (MAF>=.05 & >={MIN_SITES} sites): {idx.size:,}/{nv:,}",
          flush=True)

    # complete s-matrix for SVD: center per variant over finite sites, fill missing 0
    with np.errstate(invalid="ignore"):
        smean = np.nansum(np.where(good, s, 0.0), 0) / np.maximum(good.sum(0), 1)
    Sc = np.where(good, s - smean[None, :], 0.0)                        # [31 x nv], centered
    ScT = np.ascontiguousarray(Sc[:, idx])                             # tested family

    # secondary (direction only): raw IV slope beta_iv = Σw(x-x̄)(s-s̄)/Σw(x-x̄)²
    w = np.where(good, 1.0 / (np.maximum(se, SE_FLOOR) ** 2), 0.0)
    wT = np.ascontiguousarray(w[:, idx]); s0T = np.ascontiguousarray(np.where(good, s, 0.0)[:, idx])
    wsT = wT * s0T
    SwT = wT.sum(0); SwsT = wsT.sum(0)

    def iv_slope(xp):
        Swx = xp @ wT; Swxx = (xp * xp) @ wT; Swxs = xp @ wsT
        num = Swxs - Swx * SwsT / SwT
        varx = Swxx - Swx * Swx / SwT
        return np.divide(num, varx, out=np.zeros_like(num), where=varx > 0)

    summary = []
    for axis in axes:
        x = std_env(axis)
        for K in KS:
            # residualize the TESTED family only (SVD structure from tested loci)
            Yr = residualize(ScT, K)                    # [31 x nT], unit-norm cols
            cos_obs = Yr.T @ x                          # [nT] partial corr in [-1,1]
            a_obs = np.abs(cos_obs)
            ge = np.zeros(idx.size, np.int32)
            maxnull = np.empty(args.n_perm)
            rng = np.random.default_rng(20260702 + (cls == "snp") + K * 7 + hash(axis) % 101)
            for b in range(args.n_perm):
                xp = x[rng.permutation(NSITES)]
                ap = np.abs(Yr.T @ xp)
                maxnull[b] = ap.max()
                ge += (ap >= a_obs)
            emp_p = (1 + ge) / (1 + args.n_perm)
            q = bh_fdr(emp_p)
            thr_fwer = np.quantile(maxnull, 0.95)
            pass_fwer = a_obs >= thr_fwer
            pass_fdr = q < 0.05
            beta = iv_slope(x)                          # direction/effect size

            sub = meta_df.iloc[idx].reset_index(drop=True).copy()
            sub["cos"] = cos_obs; sub["beta_iv"] = beta
            sub["emp_p"] = emp_p; sub["fdr_q"] = q
            sub["pass_fwer"] = pass_fwer; sub["pass_fdr"] = pass_fdr
            sub["n_eff_sites"] = n_eff_sites[idx]
            hits = sub[pass_fwer | pass_fdr].copy()
            if len(hits):
                hits = lib.annotate_svs(hits, flank=2000, genes=genes)
            hits.sort_values("emp_p").to_csv(
                f"{TW}/stage2_hits_{cls}_{axis}_K{K}.csv", index=False)
            frac_p05 = float((emp_p < 0.05).mean())
            row = dict(cls=cls, axis=axis, K=K, n_loci=int(idx.size),
                       thr_fwer=float(thr_fwer), max_obs=float(a_obs.max()),
                       n_fwer=int(pass_fwer.sum()), n_fdr=int(pass_fdr.sum()),
                       n_blocks_fwer=int(sub.loc[pass_fwer, "block_clq09"].nunique()),
                       n_blocks_fdr=int(sub.loc[pass_fdr, "block_clq09"].nunique()),
                       min_emp_p=float(emp_p.min()), min_fdr_q=float(q.min()),
                       frac_emp_p_lt05=frac_p05)
            summary.append(row)
            print(f"  {cls} {axis} K{K}: n={idx.size:,} thrFWER={thr_fwer:.4f} "
                  f"maxObs={a_obs.max():.4f} nFWER={pass_fwer.sum()} nFDR={pass_fdr.sum()} "
                  f"minP={emp_p.min():.2e} minQ={q.min():.3f} frac_p<.05={frac_p05:.3f}",
                  flush=True)

    # ---- Stage-1-only selection scan ----
    p05, fdr_sig, cal = stage1_selection(cls, z, meta_df, genes, axes)
    sel = meta_df.copy()
    sel["n_sites_p05"] = p05; sel["n_sites_fdr"] = fdr_sig
    sel_tested = sel[(sel.site_maf >= MAF_MIN)].copy()
    # recurrence: expected #sites p<0.05 under null ~ Binomial(nsite_eff, 0.05)
    top_rec = sel_tested.sort_values("n_sites_fdr", ascending=False).head(200).copy()
    top_rec = lib.annotate_svs(top_rec, flank=2000, genes=genes)
    top_rec.to_csv(f"{TW}/stage1_recurrent_{cls}.csv", index=False)
    n_any_fdr = int((sel_tested.n_sites_fdr >= 1).sum())
    n_blocks_any = int(sel_tested.loc[sel_tested.n_sites_fdr >= 1, "block_clq09"].nunique())
    n_rec2 = int((sel_tested.n_sites_fdr >= 2).sum())
    print(f"\n[{cls}] STAGE-1 SELECTION: variants FDR-sig in >=1 site: {n_any_fdr:,} "
          f"({n_blocks_any:,} blocks) | in >=2 sites: {n_rec2:,}", flush=True)
    print(cal.to_string(index=False), flush=True)

    sdf = pd.DataFrame(summary)
    sdf.to_csv(f"{TW}/stage2_summary_{cls}.csv", index=False)
    stage1_stats = dict(cls=cls, n_any_fdr=n_any_fdr, n_blocks_any=n_blocks_any,
                        n_recurrent2=n_rec2,
                        median_frac_p05=float(np.nanmedian(cal["frac_p<.05"])))
    json.dump(stage1_stats, open(f"{TW}/stage1_stats_{cls}.json", "w"), indent=2)
    print(f"\n[{cls}] done ({time.time()-t0:.0f}s) -> {TW}", flush=True)


if __name__ == "__main__":
    main()

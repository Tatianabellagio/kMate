#!/usr/bin/env python
"""TEMPORAL two-stage GEA — Stage 1: per-site, per-variant temporal selection
coefficient s_site and its AMONG-PLOT drift SE, for all 31 sites.

Adapts analysis/grenenet_gea/site_variant_temporal_scoef.py to the model in its
own docstring (which the shipped code only *pooled*): here we actually fit a
per-PLOT logit-slope and use the spread of replicate plots as the drift null.

Model (per variant v, per plot j at a site, generation t):
    logit(p_{t,j,v}) ~ a_j + s_{j,v}*t ,  t in {0}∪{gens the plot was sampled}
    s_{j,v} = Σ_t (t-t̄_j) logit(p_{t,j,v}) / Σ_t (t-t̄_j)^2      (OLS slope)
  gen0 (t=0) = the SHARED SEEDMIX founding p0 (same for every plot & site).
Site statistic (replicate PLOTS = drift null: same p0, same climate):
    s_v  = mean_j s_{j,v}
    se_v = sd_j(s_{j,v}) / sqrt(n_plots)       (empirical among-plot drift SE)

A variant that moves the SAME way across independent replicate plots (small
among-plot spread vs the mean) is under selection; pure drift scatters the
plot slopes around 0.

INPUT (GLOBAL-mode kMate AF): analysis/grenenet_gea/pool_matrices/
    pool_gen{1,2,3}_{snp,nonsnp}_af.npy   [n_pools x n_variants]
    pool_gen{g}_{snp,nonsnp}.meta.csv     (site, plot, total_flowers, ...)
    af_store/p0_{snp,nonsnp}.npy          founding gen-0 AF
    af_store/index_{snp,nonsnp}.npz       chrom/pos/ref_len/alt_len

OUTPUT: <out>/stage1_{class}.npz
    chrom,pos,ref_len,alt_len,size,p0, sites[31],
    s[31 x nvar], se[31 x nvar], site_freq[31 x nvar], n_plots[31]

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY temporal_stage1.py --class nonsnp --out .../temporal_twostage
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PM = f"{lib.GEA}/pool_matrices"
STORE = lib.AF_STORE
EPS = 1e-3
GENS = (1, 2, 3)


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True, choices=["snp", "nonsnp"])
    ap.add_argument("--out", default=f"{lib.GEA}/gea_newpanel/results/temporal_twostage")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    cls = args.cls

    idx = np.load(f"{STORE}/index_{cls}.npz")
    chrom = idx["chrom"].astype("U5"); pos = idx["pos"]
    rl = idx["ref_len"].astype(np.int64); al = idx["alt_len"].astype(np.int64)
    size = np.abs(al - rl)
    p0 = np.load(f"{STORE}/p0_{cls}.npy").astype(np.float64)
    nvar = p0.shape[0]
    L0 = logit(p0).astype(np.float32)               # shared t=0 log-odds
    print(f"[{cls}] nvar={nvar:,}", flush=True)

    # load metas + open memmaps
    meta = {g: pd.read_csv(f"{PM}/pool_gen{g}_{cls}.meta.csv") for g in GENS}
    mm = {g: np.load(f"{PM}/pool_gen{g}_{cls}_af.npy", mmap_mode="r") for g in GENS}
    for g in GENS:
        meta[g]["key"] = meta[g].site.astype(str) + "_" + meta[g]["plot"].astype(str)

    sites = np.sort(np.unique(np.concatenate([meta[g].site.to_numpy() for g in GENS])))
    nS = len(sites)
    print(f"[{cls}] sites={nS}: {list(sites)}", flush=True)

    s_out = np.full((nS, nvar), np.nan, np.float32)
    se_out = np.full((nS, nvar), np.nan, np.float32)
    sf_out = np.full((nS, nvar), np.nan, np.float32)   # flower-weighted site-mean AF
    nplots = np.zeros(nS, np.int32)

    for si, s in enumerate(sites):
        # gather this site's plot rows across gens: plot -> {gen: logit_row}
        plot_L = {}                # plot -> dict(t -> logit vector)
        wsum = np.zeros(nvar, np.float64); wtot = 0.0   # for site-mean AF
        for g in GENS:
            m = meta[g]
            sel = np.where(m.site.to_numpy() == s)[0]
            if len(sel) == 0:
                continue
            rows = np.asarray(mm[g][sel]).astype(np.float32)     # [n_sel x nvar]
            w = m["total_flowers"].to_numpy(float)[sel]
            w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
            keys = m["key"].to_numpy()[sel]
            Lg = logit(rows.astype(np.float64)).astype(np.float32)
            for i, k in enumerate(keys):
                plot_L.setdefault(k, {})[float(g)] = Lg[i]
            wsum += (w[:, None] * rows.astype(np.float64)).sum(0)
            wtot += w.sum()
        sf_out[si] = (wsum / wtot).astype(np.float32) if wtot > 0 else np.nan

        # per-plot OLS logit-slope over [t=0 (shared L0)] + its sampled gens
        plot_slopes = []
        for k, tmap in plot_L.items():
            tset = np.array([0.0] + sorted(tmap.keys()))
            if tset.size < 2:
                continue
            tc = tset - tset.mean(); sst = float(np.sum(tc ** 2))
            acc = tc[0] * L0.astype(np.float64)                  # t=0 term
            for j, t in enumerate(tset[1:], start=1):
                acc = acc + tc[j] * tmap[t].astype(np.float64)
            plot_slopes.append((acc / sst).astype(np.float32))
        if len(plot_slopes) < 2:
            print(f"  site {s}: {len(plot_slopes)} plots (<2) — no SE", flush=True)
            nplots[si] = len(plot_slopes)
            if plot_slopes:
                s_out[si] = plot_slopes[0]
            continue
        S = np.vstack(plot_slopes)                               # [n_plots x nvar]
        n = S.shape[0]; nplots[si] = n
        s_out[si] = np.nanmean(S, axis=0)
        se_out[si] = np.nanstd(S, axis=0, ddof=1) / np.sqrt(n)
        print(f"  site {s}: {n} plots | median|s|={np.nanmedian(np.abs(s_out[si])):.4f}"
              f" | median se={np.nanmedian(se_out[si]):.4f}", flush=True)

    out = f"{args.out}/stage1_{cls}.npz"
    np.savez(out, chrom=chrom, pos=pos, ref_len=rl, alt_len=al, size=size,
             p0=p0.astype(np.float32), sites=sites.astype(np.int64),
             s=s_out, se=se_out, site_freq=sf_out, n_plots=nplots)
    print(f"[{cls}] -> {out}", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Per-VARIANT temporal selection coefficients at ONE site, by variant class.

The single-site analogue of build_two_stage_gea.stage1_scoef, but (i) restricted to
one common-garden SITE, (ii) kept at INDIVIDUAL-VARIANT resolution (not blocks / not
haplotypes), and (iii) SNP / small-indel / SV classes retained as a column — so we can
ask whether SVs are enriched among the selected variants/regions at that site.

Model (per variant v, per persistent plot j, generation t in {0,1,2,3}):
    logit(p_{t,j,v}) ~ a_j + s_{j,v}*t          (t=0 is the shared SEEDMIX p0)
    s_{j,v} = Sum_t (t-1.5) logit(p) / Sum_t (t-1.5)^2   (WLS-free OLS slope, weights equal)
Site statistic (the replicate PLOTS are the drift null — same p0, same climate):
    s_v  = mean_j s_{j,v}          (site selection coefficient)
    se_v = sd_j(s_{j,v}) / sqrt(n) (among-plot SE = empirical drift floor)
    z_v  = s_v / se_v ,  p_v = 2*t.sf(|z|, df=n-1)
A variant whose frequency moves the SAME direction across independent replicate plots
(small among-plot spread relative to the mean) is under selection; pure drift gives
plot slopes scattered around 0.

INPUT (GLOBAL-mode kMate AF, as requested): analysis/grenenet_selection/pool_matrices/
    pool_gen{1,2,3}_{snp,nonsnp}_af.npy  [n_pools x n_variants] float32 AF in [0,1]
    pool_gen{g}_{snp,nonsnp}.meta.csv     row-aligned pool metadata (site/plot/flowers)
    af_store/p0_{snp,nonsnp}.npy          founding (gen-0) AF, column-aligned
    af_store/index_nonsnp.npz             ref_len/alt_len for the non-SNP class split

CLASSES:  snp  (ref_len==1 & alt_len==1)
          smallindel  (non-SNP, |alt_len-ref_len| <= 50, incl equal-length MNPs)
          sv          (non-SNP, |alt_len-ref_len| >  50)

OUTPUT (--out, default analysis/grenenet_selection/site_temporal):
    site{S}_scoef_{class}.npz   chrom,pos,ref_len,alt_len,size,p0,s,se,z,pval,
                                n_plots,block,keep  (keep = passed the QC/reach filter)

CAVEAT (state in any writeup): GLOBAL mode projects AF through the INTACT founder panel,
so variants in perfect founder-LD share one trajectory and get IDENTICAL s. Per-variant
SV enrichment therefore reflects whether SELECTED FOUNDER HAPLOTYPES carry SVs, not
independent SV-level selection -> do the enrichment test at the LD-block level.

  PY=<plotting env>; $PY site_variant_temporal_scoef.py --site 4
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
from scipy.stats import t as tdist
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PM = f"{lib.GEA}/pool_matrices"
STORE = lib.AF_STORE
EPS = 1e-3                 # logit clip
SV_MIN_BP = 50
T = np.array([0.0, 1.0, 2.0, 3.0]); TC = T - T.mean(); SST = float(np.sum(TC ** 2))  # =5


def site_freq_per_gen(site: int, kind: str):
    """Flower-weighted site frequency per available generation (Pipeline-B Step 1).

    Returns (t_list, F) where t_list = [0] + [available gens in {1,2,3}] and F is
    [n_t x nvar]: row 0 = founding p0, then the flower-weighted mean AF over that
    generation's plots (NaN-aware). No persistent-plot requirement -> plot turnover
    does not drop a site; variable-length -> sites missing a generation still fit."""
    p0 = np.load(f"{STORE}/p0_{kind}.npy").astype(np.float32)
    nvar = p0.shape[0]
    t_list = [0.0]; F_rows = [p0.astype(np.float64)]
    for g in (1, 2, 3):
        m = pd.read_csv(f"{PM}/pool_gen{g}_{kind}.meta.csv")
        m4 = m[m.site == site]
        if len(m4) == 0:
            continue
        rows = m4.index.to_numpy()
        w = m4["total_flowers"].to_numpy(float)
        w = np.where(np.isfinite(w) & (w > 0), w, 1.0)[:, None]
        mat = np.load(f"{PM}/pool_gen{g}_{kind}_af.npy", mmap_mode="r")
        sub = np.asarray(mat[rows]).astype(np.float64)                 # [nplots x nvar]
        ok = np.isfinite(sub)
        num = np.nansum(np.where(ok, sub, 0.0) * w, axis=0)
        den = np.nansum(np.where(ok, w, 0.0), axis=0)
        F_rows.append(num / np.where(den > 0, den, np.nan))
        t_list.append(float(g))
    return t_list, np.vstack(F_rows), p0


def site_scoef(site: int, kind: str):
    """Per-variant temporal selection coefficient = OLS logit-slope over the site's
    available trajectory [p0, gen...] (plot-pooled, variable length)."""
    t_list, F, p0 = site_freq_per_gen(site, kind)
    t = np.array(t_list)
    if t.size < 2:
        raise ValueError(f"site {site} {kind}: <2 timepoints, cannot fit a slope")
    tc = t - t.mean(); sst = float(np.sum(tc ** 2))
    Fc = np.clip(F, EPS, 1 - EPS)
    L = np.log(Fc / (1 - Fc))                                          # [n_t x nvar]
    s = np.tensordot(tc, L, axes=(0, 0)) / sst                         # [nvar]
    return dict(p0=p0, s=s.astype(np.float32), n_gens=int(t.size - 1),
                gens=t_list[1:])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=int, default=4)
    ap.add_argument("--out", default=f"{lib.GEA}/site_temporal")
    ap.add_argument("--min-p0", type=float, default=0.02,
                    help="reachability: keep variants with min-p0 <= p0 <= 1-min-p0")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    idx_snp = np.load(f"{STORE}/index_snp.npz")
    idx_non = np.load(f"{STORE}/index_nonsnp.npz")

    # compute both class-matrices once
    res = {"snp": site_scoef(args.site, "snp"),
           "nonsnp": site_scoef(args.site, "nonsnp")}

    def emit(name, idx, R, sel):
        chrom = idx["chrom"].astype("U5")[sel]; pos = idx["pos"][sel]
        rl = idx["ref_len"].astype(np.int64)[sel]; al = idx["alt_len"].astype(np.int64)[sel]
        size = np.abs(al - rl)
        p0 = R["p0"][sel]
        reach = (p0 >= args.min_p0) & (p0 <= 1 - args.min_p0)
        keep = reach & np.isfinite(R["s"][sel])
        block = lib.assign_ld_blocks(chrom, pos)
        out = f"{args.out}/site{args.site}_scoef_{name}.npz"
        np.savez(out, chrom=chrom, pos=pos, ref_len=rl, alt_len=al, size=size,
                 p0=p0.astype(np.float32), s=R["s"][sel], block=block.astype("U16"),
                 keep=keep, n_gens=R["n_gens"])
        print(f"[{name}] {sel.sum() if sel.dtype==bool else len(pos):,} variants | "
              f"keep(reach&finite)={int(keep.sum()):,} | "
              f"median|s|(keep)={np.nanmedian(np.abs(R['s'][sel][keep])):.4f} | -> {out}")
        return keep.sum()

    all_true = np.ones(res["snp"]["s"].shape[0], bool)
    emit("snp", idx_snp, res["snp"], all_true)

    dlen = np.abs(idx_non["alt_len"].astype(np.int64) - idx_non["ref_len"].astype(np.int64))
    emit("smallindel", idx_non, res["nonsnp"], dlen <= SV_MIN_BP)
    emit("sv", idx_non, res["nonsnp"], dlen > SV_MIN_BP)
    print("done.")


if __name__ == "__main__":
    main()

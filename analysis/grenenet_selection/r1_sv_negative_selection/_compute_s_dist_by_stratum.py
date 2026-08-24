#!/usr/bin/env python
"""Precompute the DISTRIBUTIONAL s comparison (SNP/indel/SV) within initial-frequency strata, PER
SITE, for the notebook (user 2026-07-02: no threshold, no null — just the s distributions and their
shift). s = plot-replicate mean logit-slope per variant per site (plots as replicates).

For each site, each initial-frequency (p0) decile, each class: median s + q25/q75 (distribution
spread) + n, and a frequency-matched SV-vs-SNP / indel-vs-SNP median shift. Saved compact for a
30-panel (one-per-site) figure. Also a small per-variant subsample per site (all SVs + matched SNP/
indel sample) for optional violins.

Env: kmate. Writes analysis/grenenet_selection/sv_adaptive/s_dist_by_stratum.{npz,csv}.
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import importlib.util
_sp = importlib.util.spec_from_file_location(
    "plotsmod", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "_temporal_s_plots_snp_vs_nonsnp.py"))
plotsmod = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(plotsmod)

os.chdir("/global/scratch/users/tbellg/kmate")
STORE = lib.AF_STORE; PM = f"{lib.GEA}/pool_matrices"
MIN_P0, SV_BP, MIN_MAC, NBIN = 0.02, 50, 12, 10


def main():
    idx_snp = np.load(f"{STORE}/index_snp.npz"); idx_non = np.load(f"{STORE}/index_nonsnp.npz")
    ch_non = idx_non["chrom"].astype("U5"); pos_non = idx_non["pos"].astype(np.int64)
    dlen = np.abs(idx_non["alt_len"].astype(np.int64) - idx_non["ref_len"].astype(np.int64))
    is_del = idx_non["ref_len"].astype(np.int64) > idx_non["alt_len"].astype(np.int64)
    ch_snp = idx_snp["chrom"].astype("U5"); pos_snp = idx_snp["pos"].astype(np.int64)
    common_non = lib.founder_panel_keep(ch_non, pos_non, min_mac=MIN_MAC)
    common_snp = lib.founder_panel_keep(ch_snp, pos_snp, min_mac=MIN_MAC)
    p0_non = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    p0_snp = np.load(f"{STORE}/p0_snp.npy").astype(np.float64)
    snp_m = common_snp & (p0_snp >= MIN_P0) & (p0_snp <= 1 - MIN_P0)
    non_m = (dlen >= 1) & common_non & (p0_non >= MIN_P0) & (p0_non <= 1 - MIN_P0)
    cols_snp = np.where(snp_m)[0]; cols_non = np.where(non_m)[0]
    p0s = p0_snp[cols_snp]; p0n = p0_non[cols_non]
    dl = dlen[cols_non]; sv_k = dl > SV_BP; ind_k = (dl >= 1) & (dl <= SV_BP)
    p0q = np.quantile(p0s, np.linspace(0, 1, NBIN + 1)); p0q[-1] += 1e-6
    binf = lambda p: np.clip(np.digitize(p, p0q[1:-1]), 0, NBIN - 1)
    bsnp = binf(p0s); bnon = binf(p0n)

    clim = lib.load_climate()
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                    for s in pd.read_csv(f).site.unique()})
    rng = np.random.default_rng(0)

    long_rows = []            # (site,bio1,stratum,class,median,q25,q75,n)
    site_meta = []            # (site,bio1,n_plots, shift_sv, shift_ind, wil_p_sv)
    sub = {}                  # site -> dict of subsample arrays for violins
    for site in sites:
        s_sn, _, _, n = plotsmod.plot_replicate_sz(site, "snp", cols_snp, p0s)
        s_no, _, _, _ = plotsmod.plot_replicate_sz(site, "nonsnp", cols_non, p0n)
        if n < 2:
            print(f"  site {site}: {n} plots -> skip"); continue
        b1 = float(clim["bio1"].reindex([site]).iloc[0]) if site in clim.index else np.nan
        shifts_sv, shifts_ind = [], []
        for b in range(NBIN):
            sn = s_sn[bsnp == b]; no = s_no[bnon == b]
            svv = no[sv_k[bnon == b]]; inv = no[ind_k[bnon == b]]
            allb = np.concatenate([sn, no])                     # kept only for the "ALL" reference row
            base = float(np.median(sn)) if sn.size else np.nan  # per-bin SNP-only de-trending baseline
            for nm, arr in (("ALL", allb), ("SNP", sn), ("indel", inv), ("SV", svv)):
                if arr.size == 0:
                    continue
                long_rows.append(dict(site=site, bio1=round(b1, 1), stratum=b,
                                      p0_lo=round(p0q[b], 4), p0_hi=round(p0q[b + 1], 4),
                                      cls=nm, median=float(np.median(arr)),
                                      q25=float(np.quantile(arr, .25)),
                                      q75=float(np.quantile(arr, .75)), n=int(arr.size),
                                      base=base))              # per-bin SNP median for de-trending
            if svv.size:
                shifts_sv.append(np.median(svv) - base)         # SV excess vs same-freq baseline
            if inv.size:
                shifts_ind.append(np.median(inv) - base)
        # per-site frequency-matched shift = mean over strata of (median_class - median_SNP)
        wil = stats.wilcoxon(shifts_sv).pvalue if len(shifts_sv) >= 6 else np.nan
        site_meta.append(dict(site=site, bio1=round(b1, 1), n_plots=n,
                              shift_sv=round(float(np.mean(shifts_sv)), 4),
                              shift_ind=round(float(np.mean(shifts_ind)), 4),
                              wil_p_sv=round(float(wil), 4)))
        # subsample s AND p0 per class (for MAF-filterable per-site histograms). all SV; 25k SNP/indel.
        def sidx(n, k): return np.arange(n) if n <= k else rng.choice(n, k, replace=False)
        i_sn = sidx(s_sn.size, 25000)
        ii = np.where(ind_k)[0]; i_in = ii if ii.size <= 25000 else rng.choice(ii, 25000, replace=False)
        i_sv = np.where(sv_k)[0]
        sub[f"{site}_s_SNP"] = s_sn[i_sn].astype(np.float32); sub[f"{site}_p0_SNP"] = p0s[i_sn].astype(np.float32)
        sub[f"{site}_s_indel"] = s_no[i_in].astype(np.float32); sub[f"{site}_p0_indel"] = p0n[i_in].astype(np.float32)
        sub[f"{site}_s_SV"] = s_no[i_sv].astype(np.float32); sub[f"{site}_p0_SV"] = p0n[i_sv].astype(np.float32)
        print(f"  site {site:>2} b1={b1:>5.1f} plots={n}: shift SV={site_meta[-1]['shift_sv']:+.3f} "
              f"indel={site_meta[-1]['shift_ind']:+.3f} (wilcoxon p={wil:.3f})", flush=True)

    long = pd.DataFrame(long_rows); meta = pd.DataFrame(site_meta)
    long.to_csv(f"{lib.GEA}/sv_adaptive/s_dist_by_stratum.csv", index=False)
    meta.to_csv(f"{lib.GEA}/sv_adaptive/s_dist_by_stratum_sitemeta.csv", index=False)
    np.savez_compressed(f"{lib.GEA}/sv_adaptive/s_dist_by_stratum.npz",
                        p0q=p0q, **sub)
    print(f"\n[wrote] s_dist_by_stratum.{{csv,npz}} + _sitemeta.csv  ({len(meta)} sites)")
    print("\n=== per-site SV s-shift vs SNP (freq-matched, mean over strata) ===")
    print(f"median across sites: SV shift = {meta.shift_sv.median():+.4f}  "
          f"indel shift = {meta.shift_ind.median():+.4f}")
    print(f"sites with SV shifted DOWN (more negative than SNP): "
          f"{int((meta.shift_sv<0).sum())}/{len(meta)}  "
          f"(sign p={stats.binomtest((meta.shift_sv<0).sum(), len(meta)).pvalue:.4f})")


if __name__ == "__main__":
    main()

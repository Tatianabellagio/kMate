#!/usr/bin/env python
"""FOLDED site-frequency-spectrum shift per site, per class (SNP null / indel / SV), per p0-decile
bin (user 2026-07-14): does the SFS shape change (extinction / sweep-past-half / spread) more for
SVs than frequency-matched SNPs, between founding (p0) and each site's last-observed generation?

Complementary to the plot-replicate logit-slope `s` approach (_compute_s_dist_by_stratum.py): here we
track the RAW (not logit) frequency of the founder-relative MINOR allele, so it (a) needs no logit
transform (no Jensen/boundary artifact) and (b) needs no ref/alt polarization -- it folds on the p0
label, sidestepping the open insertion/deletion ancestral-polarity caveat in
SV_TEMPORAL_PURGING_SUMMARY.md entirely.

Per site: use the LAST generation (of {1,2,3}) with pool data there (same generation set
plot_replicate_sz uses); final AF per variant = mean over the site's plots at that generation
(plots as replicates, same as the existing analysis). Fold both p0 and final AF on the p0<0.5 label
("minor-at-founding" allele): folded_p0 in [0,0.5] by construction; folded_final tracks where that
SAME allele ended up (~0 = lost, >0.5 = swept past the other allele, ~1 = near-total sweep).

Per site x p0-decile x class (SNP null / indel / SV): n, median/shift of folded_final, extinction
rate (folded_final<EXT_THR), sweep-past-half rate (>0.5), near-fixation rate (>FIX_THR), variance,
and a frequency-matched-SNP-null KS distance on folded_final. Per-site SV/indel-vs-SNP mean-shift
delta + sign test across sites, mirroring _compute_s_dist_by_stratum.py's site_meta.

Missingness, two separate levels (user 2026-07-14):
  - Founder-panel completeness PER VARIANT: `lib.founder_panel_keep(called_min=0.9)` already
    requires >=90% of the 231 founders called at that (chrom,pos) -- reused unchanged here.
  - Pool coverage PER PLOT: the pool AF matrices carry no NaNs (kMate's h-projection always
    returns a value regardless of local depth), so a low-coverage plot doesn't drop out on its
    own -- it just contributes a noisier AF into the site mean. That noise lands disproportionately
    on this script's TAIL statistics (ext_rate/sweep_rate/fix_rate are far more sensitive to a
    single noisy pool than a median slope is), so plots below MIN_COV are excluded from the
    final-generation mean (checked: a 4x floor drops ~20% of plots and only 1-2 of 31 sites, no
    site left with <2 surviving plots at its final generation).

Env: kmate. Writes analysis/grenenet_gea/sv_adaptive/sfs_shift_by_site.{npz,csv} + _sitemeta.csv.
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

os.chdir("/global/scratch/users/tbellg/kmate")
STORE = lib.AF_STORE
PM = f"{lib.GEA}/pool_matrices"
MIN_P0, SV_BP, MIN_MAC, NBIN = 0.02, 50, 12, 10
EXT_THR, FIX_THR = 0.01, 0.95
MIN_COV = 4.0                    # pool-seq coverage floor (per plot); see missingness note above
GENS = (1, 2, 3)


def site_final_af(site, kind, cols):
    """Mean AF over the site's covered (>=MIN_COV) plots at its LAST available generation."""
    g_use, rows, n_dropped = None, None, 0
    for g in GENS:
        meta = pd.read_csv(f"{PM}/pool_gen{g}_{kind}.meta.csv")
        m4 = meta[meta.site == site]
        if len(m4):
            m4c = m4[m4.mean_coverage >= MIN_COV]
            if len(m4c):
                g_use, rows, n_dropped = g, m4c.index.to_numpy(), len(m4) - len(m4c)
    if g_use is None:
        return None, None, 0
    mat = np.load(f"{PM}/pool_gen{g_use}_{kind}_af.npy", mmap_mode="r")
    af = np.asarray(mat[rows][:, cols]).astype(np.float64)
    return af.mean(0), g_use, n_dropped


def main():
    idx_snp = np.load(f"{STORE}/index_snp.npz"); idx_non = np.load(f"{STORE}/index_nonsnp.npz")
    ch_non = idx_non["chrom"].astype("U5"); pos_non = idx_non["pos"].astype(np.int64)
    dlen = np.abs(idx_non["alt_len"].astype(np.int64) - idx_non["ref_len"].astype(np.int64))
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
    fold_p0_snp = np.minimum(p0s, 1 - p0s); fold_p0_non = np.minimum(p0n, 1 - p0n)
    minor_is_alt_snp = p0s < 0.5; minor_is_alt_non = p0n < 0.5

    clim = lib.load_climate()
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                    for s in pd.read_csv(f).site.unique()})

    long_rows, site_meta = [], []
    for site in sites:
        af_sn, g_sn, drop_sn = site_final_af(site, "snp", cols_snp)
        af_no, g_no, drop_no = site_final_af(site, "nonsnp", cols_non)
        if af_sn is None or af_no is None:
            print(f"  site {site}: no plot passes coverage floor -> skip"); continue
        ffin_snp = np.where(minor_is_alt_snp, af_sn, 1 - af_sn)
        ffin_non = np.where(minor_is_alt_non, af_no, 1 - af_no)
        b1 = float(clim["bio1"].reindex([site]).iloc[0]) if site in clim.index else np.nan
        shifts_sv, shifts_ind = [], []
        for b in range(NBIN):
            sn, sn0 = ffin_snp[bsnp == b], fold_p0_snp[bsnp == b]
            no, no0 = ffin_non[bnon == b], fold_p0_non[bnon == b]
            m_sv, m_in = sv_k[bnon == b], ind_k[bnon == b]
            svv, svv0 = no[m_sv], no0[m_sv]
            inv, inv0 = no[m_in], no0[m_in]
            for nm, arr, arr0 in (("SNP", sn, sn0), ("indel", inv, inv0), ("SV", svv, svv0)):
                if arr.size == 0:
                    continue
                shift = arr - arr0
                ks = np.nan
                if nm != "SNP" and sn.size >= 20 and arr.size >= 10:
                    ks = float(stats.ks_2samp(arr, sn).statistic)
                long_rows.append(dict(
                    site=site, bio1=round(b1, 1), stratum=b,
                    p0_lo=round(p0q[b], 4), p0_hi=round(p0q[b + 1], 4), cls=nm, n=int(arr.size),
                    median_final=float(np.median(arr)), mean_shift=float(np.mean(shift)),
                    var_final=float(np.var(arr)), ext_rate=float(np.mean(arr < EXT_THR)),
                    sweep_rate=float(np.mean(arr > 0.5)), fix_rate=float(np.mean(arr > FIX_THR)),
                    ks_vs_snp=ks))
            if svv.size and sn.size:
                shifts_sv.append(np.mean(svv - svv0) - np.mean(sn - sn0))
            if inv.size and sn.size:
                shifts_ind.append(np.mean(inv - inv0) - np.mean(sn - sn0))
        wil = stats.wilcoxon(shifts_sv).pvalue if len(shifts_sv) >= 6 else np.nan
        site_meta.append(dict(
            site=site, bio1=round(b1, 1), gen_snp=g_sn, gen_non=g_no,
            n_plots_dropped_lowcov=drop_sn + drop_no,
            shift_sv=round(float(np.mean(shifts_sv)), 4) if shifts_sv else np.nan,
            shift_ind=round(float(np.mean(shifts_ind)), 4) if shifts_ind else np.nan,
            wil_p_sv=round(float(wil), 4) if wil == wil else np.nan))
        print(f"  site {site:>2} b1={b1:>5.1f} gen(snp/non)={g_sn}/{g_no} "
              f"(dropped {drop_sn+drop_no} low-cov plots): "
              f"ΔSV={site_meta[-1]['shift_sv']:+.4f} Δindel={site_meta[-1]['shift_ind']:+.4f}",
              flush=True)

    long = pd.DataFrame(long_rows); meta = pd.DataFrame(site_meta)
    os.makedirs(f"{lib.GEA}/sv_adaptive", exist_ok=True)
    long.to_csv(f"{lib.GEA}/sv_adaptive/sfs_shift_by_site.csv", index=False)
    meta.to_csv(f"{lib.GEA}/sv_adaptive/sfs_shift_by_site_sitemeta.csv", index=False)
    np.savez_compressed(f"{lib.GEA}/sv_adaptive/sfs_shift_by_site.npz", p0q=p0q)
    print(f"\n[wrote] sfs_shift_by_site.{{csv,npz}} + _sitemeta.csv  ({len(meta)} sites)")
    m = meta.dropna(subset=["shift_sv"])
    print(f"median across sites: SV shift = {m.shift_sv.median():+.4f}  "
          f"indel shift = {m.shift_ind.median():+.4f}")
    nneg = int((m.shift_sv < 0).sum())
    print(f"sites with SV shifted DOWN (more negative than SNP): {nneg}/{len(m)}  "
          f"(sign p={stats.binomtest(nneg, len(m)).pvalue:.4f})")


if __name__ == "__main__":
    main()

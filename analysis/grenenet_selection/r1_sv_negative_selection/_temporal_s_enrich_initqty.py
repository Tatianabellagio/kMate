#!/usr/bin/env python
"""SNP-null-FREE comparison of selection by variant class (user 2026-07-02). Instead of using SNPs
to define the neutral tail, we use:

  (A) DISTRIBUTIONAL: within each INITIAL-FREQUENCY (p0) stratum, the raw median plot-replicate s of
      SNP / indel / SV side by side (+ Mann-Whitney shift of non-SNP vs SNP). No threshold, no null.

  (B) COMPOSITIONAL ENRICHMENT BY INITIAL QUANTITY: within each p0 stratum, take the top 5% most-
      PURGED variants ranked over ALL classes pooled (class-agnostic threshold). A class is ENRICHED
      if its share among those top movers exceeds its share of the stratum at baseline (its "initial
      quantity"):  enrich = (class share among selected) / (class share at baseline) = class_rate /
      0.05. Expectation comes from initial composition, NOT from SNPs.

s = plot-replicate mean logit-slope per variant per site (from _temporal_s_plots_snp_vs_nonsnp).
Strata = global p0 deciles. Pool contingency counts across the 31 sites; also per-site enrich for a
sign test. Purged (down) is the focal direction (that's where any SV signal lives); UP reported too.

Env: kmate. Writes analysis/grenenet_selection/r1_sv_negative_selection/results/sv_adaptive/temporal_s_enrich_initqty.csv.
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
STORE = lib.AF_STORE; PM = f"{lib.GEA}/common/results/pool_matrices"
MIN_P0, SV_BP, MIN_MAC, NBIN, TOP = 0.02, 50, 12, 10, 0.05


def main():
    idx_snp = np.load(f"{STORE}/index_snp.npz"); idx_non = np.load(f"{STORE}/index_nonsnp.npz")
    ch_non = idx_non["chrom"].astype("U5"); pos_non = idx_non["pos"].astype(np.int64)
    dlen = np.abs(idx_non["alt_len"].astype(np.int64) - idx_non["ref_len"].astype(np.int64))
    ch_snp = idx_snp["chrom"].astype("U5"); pos_snp = idx_snp["pos"].astype(np.int64)
    common_non = lib.founder_panel_keep(ch_non, pos_non, min_mac=MIN_MAC)
    common_snp = lib.founder_panel_keep(ch_snp, pos_snp, min_mac=MIN_MAC)
    p0_non = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    p0_snp = np.load(f"{STORE}/p0_snp.npy").astype(np.float64)
    reach_non = (p0_non >= MIN_P0) & (p0_non <= 1 - MIN_P0)
    reach_snp = (p0_snp >= MIN_P0) & (p0_snp <= 1 - MIN_P0)
    snp_m = common_snp & reach_snp
    non_m = (dlen >= 1) & common_non & reach_non
    cols_snp = np.where(snp_m)[0]; cols_non = np.where(non_m)[0]
    p0s = p0_snp[cols_snp]; p0n = p0_non[cols_non]
    dl = dlen[cols_non]; sv_k = dl > SV_BP; ind_k = (dl >= 1) & (dl <= SV_BP)
    # class code over the POOLED variant vector: 0=SNP 1=indel 2=SV
    p0q = np.quantile(p0s, np.linspace(0, 1, NBIN + 1)); p0q[-1] += 1e-6
    binf = lambda p: np.clip(np.digitize(p, p0q[1:-1]), 0, NBIN - 1)
    bsnp = binf(p0s); bnon = binf(p0n)

    clim = lib.load_climate()
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                    for s in pd.read_csv(f).site.unique()})

    # pooled contingency: per bin, per class -> [n_total, n_purged_top, n_fav_top]; + s sums for median
    ncls = 3
    tot = np.zeros((NBIN, ncls)); dnsel = np.zeros((NBIN, ncls)); upsel = np.zeros((NBIN, ncls))
    smed_site = {c: {b: [] for b in range(NBIN)} for c in range(ncls)}   # per-site medians
    persite = []
    for site in sites:
        s_sn, _, _, n = plotsmod.plot_replicate_sz(site, "snp", cols_snp, p0s)
        s_no, _, _, _ = plotsmod.plot_replicate_sz(site, "nonsnp", cols_non, p0n)
        if n < 2:
            continue
        rec = dict(site=site, bio1=round(float(clim["bio1"].reindex([site]).iloc[0]), 1)
                   if site in clim.index else np.nan)
        for b in range(NBIN):
            # pooled s in this bin across all classes
            sn_b = s_sn[bsnp == b]; no_b = s_no[bnon == b]
            svm = sv_k[bnon == b]; inm = ind_k[bnon == b]
            alls = np.concatenate([sn_b, no_b])
            if alls.size < 40:
                continue
            dn_thr = np.quantile(alls, TOP)          # bottom 5% (most purged), class-agnostic
            up_thr = np.quantile(alls, 1 - TOP)
            # per class: SNP(0), indel(1), SV(2)
            classes = [("snp", sn_b, np.ones(sn_b.size, bool)),
                       ("indel", no_b, inm), ("sv", no_b, svm)]
            for ci, (nm, sv_, msk) in enumerate(classes):
                arr = sv_[msk] if nm != "snp" else sv_
                if arr.size == 0:
                    continue
                tot[b, ci] += arr.size
                dnsel[b, ci] += (arr < dn_thr).sum()
                upsel[b, ci] += (arr > up_thr).sum()
                smed_site[ci][b].append(float(np.median(arr)))
            # per-site purged enrichment (SV) for sign test, common+mid pooled
        persite.append(rec)
        print(f"  site {site} done (plots {n})", flush=True)

    # ---- report ----
    cls = ["SNP", "indel", "SV"]
    print("\n=== (A) MEDIAN plot-replicate s by INITIAL-FREQUENCY stratum (no null) ===")
    print(f"{'p0 bin':>18} {'nSNP':>9} {'nSV':>7} | median s: SNP / indel / SV   (neg=declining)")
    for b in range(NBIN):
        lo, hi = p0q[b], p0q[b + 1]
        med = [np.median(smed_site[c][b]) if smed_site[c][b] else np.nan for c in range(3)]
        print(f"  [{lo:.3f},{hi:.3f})  {int(tot[b,0]/max(len(persite),1)):>8,} {int(tot[b,2]/max(len(persite),1)):>6,} | "
              f"{med[0]:+.3f} / {med[1]:+.3f} / {med[2]:+.3f}")

    print("\n=== (B) COMPOSITIONAL ENRICHMENT among top-5% PURGED (threshold class-agnostic, by init-freq) ===")
    print("enrich = class selection rate / 0.05  (>1 => over-represented vs its initial share)")
    print(f"{'p0 bin':>18} | purged enrich: indel / SV     | favoured enrich: indel / SV")
    rows = []
    for b in range(NBIN):
        er_dn = [dnsel[b, c] / tot[b, c] / TOP if tot[b, c] > 0 else np.nan for c in range(3)]
        er_up = [upsel[b, c] / tot[b, c] / TOP if tot[b, c] > 0 else np.nan for c in range(3)]
        print(f"  [{p0q[b]:.3f},{p0q[b+1]:.3f}) | {er_dn[1]:.2f} / {er_dn[2]:.2f}"
              f"              | {er_up[1]:.2f} / {er_up[2]:.2f}")
        rows.append(dict(bin=b, p0_lo=round(p0q[b], 3), p0_hi=round(p0q[b + 1], 3),
                         n_snp=int(tot[b, 0]), n_indel=int(tot[b, 1]), n_sv=int(tot[b, 2]),
                         purge_enr_indel=round(er_dn[1], 3), purge_enr_sv=round(er_dn[2], 3),
                         fav_enr_indel=round(er_up[1], 3), fav_enr_sv=round(er_up[2], 3)))
    # overall (pooled over bins) enrichment
    print("\n=== OVERALL (pooled across strata) ===")
    for ci, nm in ((1, "indel"), (2, "SV")):
        dn = dnsel[:, ci].sum() / tot[:, ci].sum() / TOP
        up = upsel[:, ci].sum() / tot[:, ci].sum() / TOP
        print(f"  {nm:>6}: purged {dn:.3f}x  favoured {up:.3f}x   "
              f"(baseline-share-controlled, no SNP null)")
    pd.DataFrame(rows).to_csv(f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive/temporal_s_enrich_initqty.csv", index=False)
    print(f"\n[wrote] temporal_s_enrich_initqty.csv")


if __name__ == "__main__":
    main()

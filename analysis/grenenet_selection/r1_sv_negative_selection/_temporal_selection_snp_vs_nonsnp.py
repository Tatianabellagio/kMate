#!/usr/bin/env python
"""PURE TEMPORAL, per-variant, per-site (user 2026-07-02): forget haploblocks/GWAS. For each
variant, selection = its projected allele-frequency CHANGE over the generations (Δp = AF_lastgen -
p0). Question: do NON-SNPs (indels+SVs) fall "under selection" MORE OFTEN than SNPs?

FREQUENCY-MATCHED by construction (rarer variants drift harder -> must match p0). The neutral drift
envelope is defined FROM THE SNPs within each starting-frequency (p0) decile: a variant is "under
selection" at a site if its |Δp| exceeds the TAILQ percentile of the SNP |Δp| distribution in its
own p0 bin. By definition (1-TAILQ) of SNPs are outliers per bin; ENRICHMENT = (non-SNP outlier
rate) / (SNP outlier rate) at matched p0. >1 => non-SNPs more often under selection than same-
frequency SNPs. Also split by DIRECTION (up = favoured, down = purged) and by sub-class (indel/SV).

Per site (run once): outlier ENRICHMENT for nonsnp/indel/SV vs the SNP-defined p0-matched null,
two-sided |Δp| and directional. Aggregate across the 31 sites: median fold + sign test + pooled.

CAVEAT (unavoidable here): pool AF is GLOBAL-mode founder projection, so a non-SNP and a SNP on the
same founders share ONE trajectory. This test therefore measures "are non-SNPs found on the more
strongly-MOVING haplotypes than matched SNPs" — it CANNOT separate non-SNP-SPECIFIC selection from
non-SNPs riding selected haplotypes (passenger co-occurrence). It is the honest answer to the literal
frequency question; an independent per-variant SV AF (local-mode/vg) is the only way past that.

Env: kmate. Writes analysis/grenenet_selection/sv_adaptive/temporal_selection_snp_vs_nonsnp.csv (per site)
and _summary.csv (across sites).
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from site_variant_temporal_scoef import site_freq_per_gen

os.chdir("/global/scratch/users/tbellg/kmate")
STORE = lib.AF_STORE
PM = f"{lib.GEA}/pool_matrices"
MIN_P0 = 0.02
SV_BP = 50
MIN_MAC = 12
NBIN = 10
TAILQ = 0.95   # SNP within-bin percentile of |Δp| defining the "under selection" tail (5% of SNPs)


def enrich(dv, p0v, snp_thr, binf):
    """outlier rate of a focal category vs the SNP-defined per-bin threshold, matched by p0.
    Returns (obs_rate, exp_rate, fold, n) two-sided (|Δp|>thr)."""
    b = binf(p0v)
    thr = snp_thr[b]
    out = np.abs(dv) > thr
    obs = out.mean()
    exp = (1 - TAILQ)              # SNP tail rate by construction
    return float(obs), float(exp), float(obs / exp if exp > 0 else np.nan), int(len(dv))


def dir_rate(dv, p0v, snp_up, snp_dn, binf):
    """fraction in the UP tail (Δp>snp_up_thr) and DOWN tail (Δp<snp_dn_thr), p0-matched."""
    b = binf(p0v)
    up = (dv > snp_up[b]).mean()
    dn = (dv < snp_dn[b]).mean()
    return float(up), float(dn)


def main():
    idx_snp = np.load(f"{STORE}/index_snp.npz")
    idx_non = np.load(f"{STORE}/index_nonsnp.npz")
    ch_non = idx_non["chrom"].astype("U5"); pos_non = idx_non["pos"].astype(np.int64)
    rl = idx_non["ref_len"].astype(np.int64); al = idx_non["alt_len"].astype(np.int64)
    dlen = np.abs(al - rl)
    ch_snp = idx_snp["chrom"].astype("U5"); pos_snp = idx_snp["pos"].astype(np.int64)

    print(f"founder-panel common masks (MAC>={MIN_MAC}) ...")
    common_non = lib.founder_panel_keep(ch_non, pos_non, min_mac=MIN_MAC)
    common_snp = lib.founder_panel_keep(ch_snp, pos_snp, min_mac=MIN_MAC)
    p0_non = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    p0_snp = np.load(f"{STORE}/p0_snp.npy").astype(np.float64)
    reach_non = (p0_non >= MIN_P0) & (p0_non <= 1 - MIN_P0)
    reach_snp = (p0_snp >= MIN_P0) & (p0_snp <= 1 - MIN_P0)
    nonsnp_m = (dlen >= 1) & common_non & reach_non
    sv_m = (dlen > SV_BP) & common_non & reach_non
    ind_m = (dlen >= 1) & (dlen <= SV_BP) & common_non & reach_non
    snp_m = common_snp & reach_snp
    print(f"common+reachable: SNP={int(snp_m.sum()):,}  non-SNP={int(nonsnp_m.sum()):,}  "
          f"indel={int(ind_m.sum()):,}  SV={int(sv_m.sum()):,}")

    clim = lib.load_climate()
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                    for s in pd.read_csv(f).site.unique()})
    p0q = np.quantile(p0_snp[snp_m], np.linspace(0, 1, NBIN + 1)); p0q[-1] += 1e-6
    binf = lambda p: np.clip(np.digitize(p, p0q[1:-1]), 0, NBIN - 1)

    rows = []
    for site in sites:
        try:
            _, Fn, _ = site_freq_per_gen(site, "nonsnp")
            _, Fs, _ = site_freq_per_gen(site, "snp")
        except Exception as e:
            print(f"  site {site}: skip ({e})"); continue
        dp_non = np.clip(Fn[-1], 0, 1) - p0_non
        dp_snp = np.clip(Fs[-1], 0, 1) - p0_snp
        b1 = float(clim["bio1"].reindex([site]).iloc[0]) if site in clim.index else np.nan
        # SNP-defined per-bin thresholds (the neutral drift envelope)
        ds = dp_snp[snp_m]; bs = binf(p0_snp[snp_m])
        snp_thr = np.array([np.quantile(np.abs(ds[bs == k]), TAILQ) if (bs == k).any() else np.inf
                            for k in range(NBIN)])
        snp_up = np.array([np.quantile(ds[bs == k], TAILQ) if (bs == k).any() else np.inf
                           for k in range(NBIN)])
        snp_dn = np.array([np.quantile(ds[bs == k], 1 - TAILQ) if (bs == k).any() else -np.inf
                           for k in range(NBIN)])
        rec = dict(site=site, bio1=round(b1, 1), last_gen=Fn.shape[0] - 1)
        # SNP self-check (should be ~1.0 fold, ~5% each tail)
        for name, m in (("snp", snp_m),):
            o, e, f, n = enrich(dp_snp[m], p0_snp[m], snp_thr, binf)
            rec["snp_tailrate"] = round(o, 4)
        for name, m, dp, p0 in (("nonsnp", nonsnp_m, dp_non, p0_non),
                                ("indel", ind_m, dp_non, p0_non),
                                ("sv", sv_m, dp_non, p0_non)):
            o, e, f, n = enrich(dp[m], p0[m], snp_thr, binf)
            up, dn = dir_rate(dp[m], p0[m], snp_up, snp_dn, binf)
            rec[f"{name}_fold"] = round(f, 3)
            rec[f"{name}_up_fold"] = round(up / (1 - TAILQ), 3)
            rec[f"{name}_dn_fold"] = round(dn / (1 - TAILQ), 3)
            rec[f"{name}_n"] = n
        rows.append(rec)
        print(f"  site {site:>2} bio1={b1:>5.1f} (SNP tail {100*rec['snp_tailrate']:.1f}%): "
              f"nonSNP {rec['nonsnp_fold']:.2f}x (up {rec['nonsnp_up_fold']:.2f} dn {rec['nonsnp_dn_fold']:.2f}) | "
              f"indel {rec['indel_fold']:.2f}x  SV {rec['sv_fold']:.2f}x")

    df = pd.DataFrame(rows).sort_values("bio1")
    out = f"{lib.GEA}/sv_adaptive/temporal_selection_snp_vs_nonsnp.csv"
    df.to_csv(out, index=False)
    print(f"\n[wrote] {out}")

    print("\n=== ACROSS 31 SITES: is the non-SNP category under selection MORE than matched SNPs? ===")
    print(f"(fold = category outlier rate / SNP outlier rate at matched p0; 1.0 = same as SNPs)\n")
    srows = []
    for name in ("nonsnp", "indel", "sv"):
        for tail in ("fold", "up_fold", "dn_fold"):
            col = df[f"{name}_{tail}"]
            med = col.median(); nsg = int((col > 1).sum())
            # sign test vs 0.5
            p = stats.binomtest(nsg, len(col), 0.5).pvalue
            srows.append(dict(category=name, tail=tail, median_fold=round(med, 3),
                              sites_gt1=f"{nsg}/{len(col)}", sign_p=round(p, 4)))
            lbl = {"fold": "|Δp| (either dir)", "up_fold": "UP (favoured)", "dn_fold": "DOWN (purged)"}[tail]
            print(f"  {name:>6} {lbl:<18}: median {med:.2f}x  |  >1 at {nsg}/{len(col)} sites  (sign p={p:.3f})")
    pd.DataFrame(srows).to_csv(
        f"{lib.GEA}/sv_adaptive/temporal_selection_snp_vs_nonsnp_summary.csv", index=False)
    print("\n>1 & most sites & sign-p<0.05 => category IS under selection more often than matched SNPs.")
    print("~1 => behaves like frequency-matched SNPs (no category-specific selection signal).")


if __name__ == "__main__":
    main()

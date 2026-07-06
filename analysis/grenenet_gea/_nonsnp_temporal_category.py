#!/usr/bin/env python
"""NON-SNP as a WHOLE CATEGORY vs SNP (user 2026-07-02): the SV-specific temporal test asked about
>50bp variants only; here we pool ALL non-SNP records (indels+SVs, |alt_len-ref_len|>=1) into one
category and ask the same question — does the non-SNP category move (change frequency / show
selection) differently than frequency-matched SNPs, per site and along the climate gradient?

Mirrors _sv_temporal_direct.py exactly (same p0-matching, same global-mode projected AF), but the
focal class is ALL-NON-SNP. Also reports the three sub-classes (indel 1-50 / SV >50 / all-nonSNP)
side by side so the category is decomposed.

  dp[v] = clip(AF_lastgen,0,1) - p0[v]           per-variant projected site AF change from founding
  focal - matched-SNP  = per-p0-bin median-signed dp difference, focal-count weighted.

CAVEAT (same as SV test): pool_matrices AF is GLOBAL-mode founder projection; this is a consistency
readout of the founder-h trajectories, NOT an independent instrument (that needs local-mode/vg AF).

Env: kmate.  Writes results/grenenet_gea/sv_adaptive/nonsnp_temporal_category.csv.
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from site_variant_temporal_scoef import site_freq_per_gen

os.chdir("/global/scratch/users/tbellg/kmate")
STORE = lib.AF_STORE
PM = f"{lib.GEA}/pool_matrices"
MIN_P0 = 0.02
SV_BP = 50
MIN_MAC = 12


def matched_diff(dfocal, p0focal, dsnp, p0snp, p0q):
    """per-p0-bin median-signed dp of focal minus matched SNP, focal-count weighted."""
    binv = lambda p: np.clip(np.digitize(p, p0q[1:-1]), 0, 9)
    bf = binv(p0focal); bs = binv(p0snp)
    med_s = np.array([np.median(dsnp[bs == k]) if (bs == k).any() else np.nan for k in range(10)])
    med_f = np.array([np.median(dfocal[bf == k]) if (bf == k).any() else np.nan for k in range(10)])
    w = np.array([(bf == k).sum() for k in range(10)], float)
    ok = np.isfinite(med_f) & np.isfinite(med_s) & (w > 0)
    if not ok.any():
        return np.nan, np.nan
    matched = float(np.sum(med_s[ok] * w[ok]) / w[ok].sum())
    diff = float(np.sum((med_f[ok] - med_s[ok]) * w[ok]) / w[ok].sum())
    return matched, diff


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
    # focal classes over the non-SNP index
    nonsnp_m = (dlen >= 1) & common_non & reach_non      # WHOLE non-SNP category
    sv_m = (dlen > SV_BP) & common_non & reach_non
    ind_m = (dlen >= 1) & (dlen <= SV_BP) & common_non & reach_non
    snp_m = common_snp & reach_snp
    print(f"common+reachable: non-SNP(all)={int(nonsnp_m.sum()):,}  indel={int(ind_m.sum()):,}  "
          f"SV={int(sv_m.sum()):,}  SNP={int(snp_m.sum()):,}")

    clim = lib.load_climate()
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                    for s in pd.read_csv(f).site.unique()})
    p0q = np.quantile(p0_snp[snp_m], np.linspace(0, 1, 11)); p0q[-1] += 1e-6

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
        rec = dict(site=site, bio1=round(b1, 1), last_gen=Fn.shape[0] - 1)
        for name, m in (("nonsnp", nonsnp_m), ("indel", ind_m), ("sv", sv_m)):
            matched, diff = matched_diff(dp_non[m], p0_non[m], dp_snp[snp_m], p0_snp[snp_m], p0q)
            rec[f"med_dp_{name}"] = round(float(np.median(dp_non[m])), 5)
            rec[f"{name}_minus_snp"] = round(diff, 5)
            rec[f"frac_{name}_down"] = round(float((dp_non[m] < 0).mean()), 3)
        rec["matched_snp_base"] = round(matched, 5)  # matched-SNP baseline (same for all foci)
        rows.append(rec)
        print(f"  site {site:>2} bio1={b1:>5.1f}: nonSNP-SNP={rec['nonsnp_minus_snp']:+.4f} "
              f"(indel {rec['indel_minus_snp']:+.4f}, SV {rec['sv_minus_snp']:+.4f}) | "
              f"{100*rec['frac_nonsnp_down']:.0f}% nonSNP down")

    df = pd.DataFrame(rows).sort_values("bio1")
    out = f"{lib.GEA}/sv_adaptive/nonsnp_temporal_category.csv"
    df.to_csv(out, index=False)
    print(f"\n[wrote] {out}")

    m = np.isfinite(df.bio1)
    print("\n=== SUMMARY: non-SNP category vs frequency-matched SNP ===")
    for name in ("nonsnp", "indel", "sv"):
        col = df[f"{name}_minus_snp"]
        r = stats.spearmanr(df.bio1[m], col[m])
        print(f"  {name:>6}: median {name}-SNP dp = {col.median():+.5f} | "
              f"more-negative-than-SNP at {int((col<0).sum())}/{len(df)} sites | "
              f"corr(bio1,{name}-SNP)= {r.statistic:+.2f} (p={r.pvalue:.3f})")
    print("\n(-ve median & majority of sites => category declines MORE than matched SNPs = category-"
          "specific selection; ~0 & ~half sites => category behaves like matched SNPs.)")


if __name__ == "__main__":
    main()

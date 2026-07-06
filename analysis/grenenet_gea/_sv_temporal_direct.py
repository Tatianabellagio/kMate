#!/usr/bin/env python
"""DIRECT temporal audit: do common SVs actually go DOWN in frequency over generations,
per site, and more so at hot sites? Grounds the founder-h beta result (_sv_winning_genetics.py)
in real per-variant projected AF trajectories.

For each site, per variant class (SV>50bp / small-indel / SNP):
    dp[v] = clip(AF_lastgen,0,1) - p0[v]                  (projected site AF change from founding)
Filter to founder-panel COMMON variants (MAC>=12, matching the beta test) that are pool-reachable
(0.02<=p0<=0.98). Compare SVs to p0-MATCHED SNPs (SVs start rarer -> a rare variant makes a bigger
frequency move; bin by p0, take per-bin median signed dp, SV-count-weight). Report per-site
median dp_SV, matched SNP baseline, their difference, and the correlation with site temperature.

AUDIT SPLITS:
  - insertion (alt>ref) vs deletion (alt<ref): if only one sign moves, suspect reference/mapping
    bias rather than 'SVs decline'.
  - example trajectories at the hottest site (concrete numbers).

CAVEAT: pool_matrices AF is GLOBAL-mode founder projection (dp[v] is algebraically tied to the
per-founder h change beta uses). So agreement here VALIDATES the pipeline + shows the real
magnitudes/trajectories; it is NOT an independent instrument (that needs local-mode / vg SV AF).

Env: kmate.  Writes results/grenenet_gea/sv_adaptive/sv_temporal_direct.csv.
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


def main():
    idx_snp = np.load(f"{STORE}/index_snp.npz")
    idx_non = np.load(f"{STORE}/index_nonsnp.npz")
    ch_non = idx_non["chrom"].astype("U5"); pos_non = idx_non["pos"].astype(np.int64)
    rl = idx_non["ref_len"].astype(np.int64); al = idx_non["alt_len"].astype(np.int64)
    dlen = np.abs(al - rl)
    ch_snp = idx_snp["chrom"].astype("U5"); pos_snp = idx_snp["pos"].astype(np.int64)

    print("founder-panel common masks (MAC>=12) ...")
    common_non = lib.founder_panel_keep(ch_non, pos_non, min_mac=12)
    common_snp = lib.founder_panel_keep(ch_snp, pos_snp, min_mac=12)
    p0_non = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    p0_snp = np.load(f"{STORE}/p0_snp.npy").astype(np.float64)
    reach_non = (p0_non >= MIN_P0) & (p0_non <= 1 - MIN_P0)
    reach_snp = (p0_snp >= MIN_P0) & (p0_snp <= 1 - MIN_P0)
    sv_m = (dlen > SV_BP) & common_non & reach_non
    ind_m = (dlen <= SV_BP) & (dlen >= 1) & common_non & reach_non
    snp_m = common_snp & reach_snp
    is_del = (rl > al); is_ins = (al > rl)
    print(f"common+reachable: SV={int(sv_m.sum()):,} indel={int(ind_m.sum()):,} SNP={int(snp_m.sum()):,}")

    clim = lib.load_climate()
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                    for s in pd.read_csv(f).site.unique()})
    # global p0 deciles for matching (from the reachable SNP pool)
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
        # p0-matched signed dp: per p0-bin median, SV-weighted
        binv = lambda p: np.clip(np.digitize(p, p0q[1:-1]), 0, 9)
        bsv = binv(p0_non[sv_m]); bsnp = binv(p0_snp[snp_m])
        dsv = dp_non[sv_m]; dsn = dp_snp[snp_m]
        med_snp_bin = np.array([np.median(dsn[bsnp == k]) if (bsnp == k).any() else np.nan for k in range(10)])
        med_sv_bin = np.array([np.median(dsv[bsv == k]) if (bsv == k).any() else np.nan for k in range(10)])
        wsv = np.array([(bsv == k).sum() for k in range(10)], float)
        ok = np.isfinite(med_sv_bin) & np.isfinite(med_snp_bin) & (wsv > 0)
        matched_snp = float(np.sum(med_snp_bin[ok] * wsv[ok]) / wsv[ok].sum())
        med_sv = float(np.median(dsv))
        diff = float(np.sum((med_sv_bin[ok] - med_snp_bin[ok]) * wsv[ok]) / wsv[ok].sum())  # SV - matched SNP
        # ins/del split (raw median dp)
        dp_del = float(np.median(dp_non[sv_m & is_del])) if (sv_m & is_del).any() else np.nan
        dp_ins = float(np.median(dp_non[sv_m & is_ins])) if (sv_m & is_ins).any() else np.nan
        frac_sv_down = float((dsv < 0).mean())
        rows.append(dict(site=site, bio1=round(b1, 1), last_gen=Fn.shape[0] - 1,
                         med_dp_sv=round(med_sv, 5), matched_snp=round(matched_snp, 5),
                         sv_minus_snp=round(diff, 5), frac_sv_down=round(frac_sv_down, 3),
                         dp_deletion=round(dp_del, 5), dp_insertion=round(dp_ins, 5)))
        print(f"  site {site:>2} bio1={b1:>5.1f}: med dp_SV={med_sv:+.4f} (matched SNP {matched_snp:+.4f}, "
              f"SV-SNP={diff:+.4f}) | {100*frac_sv_down:.0f}% SVs down | del{dp_del:+.4f} ins{dp_ins:+.4f}")

    df = pd.DataFrame(rows).sort_values("bio1")
    out = f"{lib.GEA}/sv_adaptive/sv_temporal_direct.csv"
    df.to_csv(out, index=False)
    print(f"\n[wrote] {out}")

    # cross-site summary + climate dependence
    m = np.isfinite(df.bio1)
    print("\n=== SUMMARY ===")
    print(f"sites where SV-SNP dp < 0 (SVs decline vs matched SNP): {int((df.sv_minus_snp<0).sum())}/{len(df)}")
    print(f"median across sites: dp_SV={df.med_dp_sv.median():+.4f}  SV-SNP={df.sv_minus_snp.median():+.4f}")
    r = stats.spearmanr(df.bio1[m], df.sv_minus_snp[m])
    print(f"climate dependence: corr(bio1, SV-SNP dp) = {r.statistic:+.2f} (p={r.pvalue:.3f}) "
          "-> -ve = SVs decline MORE (relative to SNPs) at HOT sites")
    # ins vs del consistency (reference-bias audit)
    print(f"deletion median dp across sites: {df.dp_deletion.median():+.4f} | "
          f"insertion: {df.dp_insertion.median():+.4f}  "
          f"(both negative => not a pure ref/mapping-direction artifact)")

    # concrete example trajectories at the hottest site
    hot = int(df.sort_values('bio1').iloc[-1].site)
    print(f"\n=== example SV trajectories at hottest site {hot} ===")
    _, Fn, _ = site_freq_per_gen(hot, "nonsnp")
    dp = np.clip(Fn[-1], 0, 1) - p0_non
    cand = np.where(sv_m)[0]
    worst = cand[np.argsort(dp[cand])[:12]]
    traj = pd.DataFrame(dict(chrom=ch_non[worst], pos=pos_non[worst],
                             ref_len=rl[worst], alt_len=al[worst], p0=np.round(p0_non[worst], 3),
                             **{f"gen{g}": np.round(Fn[i + 1][worst], 3) for i, g in enumerate(range(1, Fn.shape[0]))},
                             dp=np.round(dp[worst], 3)))
    print(traj.to_string(index=False))


if __name__ == "__main__":
    main()

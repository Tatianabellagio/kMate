#!/usr/bin/env python
"""SELECTION COEFFICIENT (not Δp) SNP vs non-SNP, per site (user 2026-07-02). With multi-generation
trajectories we estimate a true per-generation selection coefficient:

    s_v = OLS slope of logit(p_g) on generation g,  g = 0..G   (log-odds change per generation)

Under directional selection log-odds move linearly at rate s (haploid/clonal model; GrENE-net is
~97% selfing so heterozygosity is negligible and this is the right model). Δp saturates near 0/1;
s does not -> s is the frequency-fair selection statistic. gen0 = founding (p0, SEEDMIX projection),
gen1..G = projected per-gen AF (pool_matrices). freqs clipped to [1e-3,1-1e-3] for logit.

Then the SAME frequency-matched SNP-empirical test as _temporal_sel_drift_maf.py, but on s:
  - p0 deciles (MAF control); within each bin the SNP s-distribution defines the tail.
  - DOWN(purged) tail = s below SNP 5th pctile; UP(favoured) = s above SNP 95th; fold = category
    tail rate / 5%. MAF bands (rare/mid/common) + SV ins/del split. Aggregate over 31 sites.
  - Also raw median |s| and signed s per class (magnitude of selection, descriptive).

CAVEAT unchanged: pool AF is a GLOBAL-mode founder projection -> s of a variant is a projection of
founder-h slopes; SNP-on-same-founders shares it. SNPs = realized-drift null; cannot separate SV-
specific selection from riding selected haplotypes (needs independent local-mode/vg SV AF).

Env: kmate. Writes analysis/grenenet_selection/sv_adaptive/temporal_s_snp_vs_nonsnp.csv (+_summary.csv).
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
MIN_P0, SV_BP, MIN_MAC, NBIN, TAILQ, EPS = 0.02, 50, 12, 10, 0.95, 1e-3
BANDS = [("rare", 0.02, 0.10), ("mid", 0.10, 0.20), ("common", 0.20, 0.501)]


def logit_slope(P, p0):
    """s per variant = OLS slope of logit(freq) on generation. P=[G,nvar] gens 1..G; p0=[nvar] gen0.
    Returns s (nvar,)."""
    traj = np.vstack([p0[None, :], P])                 # (G+1, nvar): gens 0..G
    y = np.log(np.clip(traj, EPS, 1 - EPS) / (1 - np.clip(traj, EPS, 1 - EPS)))
    x = np.arange(traj.shape[0], dtype=float)
    xc = x - x.mean(); Sxx = (xc ** 2).sum()
    return (xc @ y) / Sxx                                # (nvar,)


def down_up_fold(s, p0, snp_dn, snp_up, binf, band):
    maf = np.minimum(p0, 1 - p0)
    sel = (maf >= band[1]) & (maf < band[2])
    if sel.sum() < 20:
        return np.nan, np.nan, int(sel.sum())
    sv = s[sel]; b = binf(p0[sel])
    dn = (sv < snp_dn[b]).mean() / (1 - TAILQ)
    up = (sv > snp_up[b]).mean() / (1 - TAILQ)
    return float(dn), float(up), int(sel.sum())


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
    reach_non = (p0_non >= MIN_P0) & (p0_non <= 1 - MIN_P0)
    reach_snp = (p0_snp >= MIN_P0) & (p0_snp <= 1 - MIN_P0)
    nonsnp_m = (dlen >= 1) & common_non & reach_non
    sv_m = (dlen > SV_BP) & common_non & reach_non
    ind_m = (dlen >= 1) & (dlen <= SV_BP) & common_non & reach_non
    snp_m = common_snp & reach_snp
    print(f"SNP={int(snp_m.sum()):,} nonSNP={int(nonsnp_m.sum()):,} indel={int(ind_m.sum()):,} "
          f"SV={int(sv_m.sum()):,}")

    clim = lib.load_climate()
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                    for s in pd.read_csv(f).site.unique()})
    p0q = np.quantile(p0_snp[snp_m], np.linspace(0, 1, NBIN + 1)); p0q[-1] += 1e-6
    binf = lambda p: np.clip(np.digitize(p, p0q[1:-1]), 0, NBIN - 1)

    rows = []
    for site in sites:
        try:
            _, Fn, _ = site_freq_per_gen(site, "nonsnp"); _, Fs, _ = site_freq_per_gen(site, "snp")
        except Exception as e:
            print(f"  site {site}: skip ({e})"); continue
        G = Fn.shape[0] - 1
        s_non = logit_slope(Fn[1:], p0_non); s_snp = logit_slope(Fs[1:], p0_snp)
        b1 = float(clim["bio1"].reindex([site]).iloc[0]) if site in clim.index else np.nan
        ss = s_snp[snp_m]; bs = binf(p0_snp[snp_m])
        snp_dn = np.array([np.quantile(ss[bs == k], 1 - TAILQ) if (bs == k).any() else -np.inf
                           for k in range(NBIN)])
        snp_up = np.array([np.quantile(ss[bs == k], TAILQ) if (bs == k).any() else np.inf
                           for k in range(NBIN)])
        rec = dict(site=site, bio1=round(b1, 1), gens=G,
                   med_s_snp=round(float(np.median(ss)), 4),
                   med_s_indel=round(float(np.median(s_non[ind_m])), 4),
                   med_s_sv=round(float(np.median(s_non[sv_m])), 4))
        for name, m in (("nonsnp", nonsnp_m), ("indel", ind_m), ("sv", sv_m)):
            for band in BANDS:
                dn, up, n = down_up_fold(s_non[m], p0_non[m], snp_dn, snp_up, binf, band)
                rec[f"{name}_{band[0]}_dn"] = round(dn, 3) if dn == dn else np.nan
                rec[f"{name}_{band[0]}_up"] = round(up, 3) if up == up else np.nan
        for tag, mm in (("svdel", sv_m & is_del), ("svins", sv_m & ~is_del)):
            dn, up, n = down_up_fold(s_non[mm], p0_non[mm], snp_dn, snp_up, binf, ("all", 0.02, 0.501))
            rec[f"{tag}_dn"] = round(dn, 3) if dn == dn else np.nan
        rows.append(rec)
        print(f"  site {site:>2} b1={b1:>5.1f} g={G}: med_s SNP={rec['med_s_snp']:+.3f} "
              f"indel={rec['med_s_indel']:+.3f} SV={rec['med_s_sv']:+.3f} | SV-down mid/com="
              f"{rec.get('sv_mid_dn')}/{rec.get('sv_common_dn')} (del {rec.get('svdel_dn')} "
              f"ins {rec.get('svins_dn')})", flush=True)

    df = pd.DataFrame(rows).sort_values("bio1")
    df.to_csv(f"{lib.GEA}/sv_adaptive/temporal_s_snp_vs_nonsnp.csv", index=False)

    def sgn(col):
        c = df[col].dropna(); n = int((c > 1).sum())
        return round(c.median(), 3), f"{n}/{len(c)}", round(stats.binomtest(n, len(c), 0.5).pvalue, 4)

    print("\n=== SELECTION COEFFICIENT s: is non-SNP under selection more than matched SNPs? ===")
    print("median s per class across sites (log-odds/gen; negative = declining):")
    print(f"  SNP {df.med_s_snp.median():+.4f} | indel {df.med_s_indel.median():+.4f} | "
          f"SV {df.med_s_sv.median():+.4f}")
    print("\nfrequency-matched tail fold (category rate / SNP 5%); DOWN=purged, UP=favoured:")
    for name in ["nonsnp", "indel", "sv"]:
        for band in ["rare", "mid", "common"]:
            md, sd, pd_ = sgn(f"{name}_{band}_dn"); mu, su, pu = sgn(f"{name}_{band}_up")
            print(f"  {name:>6} {band:<7}: DOWN {md}x ({sd}, p={pd_}) | UP {mu}x ({su}, p={pu})")
    for tag in ["svdel", "svins"]:
        m, s, p = sgn(f"{tag}_dn"); print(f"  {tag:>6} DOWN(all): {m}x ({s}, p={p})")
    df.to_csv(f"{lib.GEA}/sv_adaptive/temporal_s_snp_vs_nonsnp.csv", index=False)
    print(f"\n[wrote] temporal_s_snp_vs_nonsnp.csv")


if __name__ == "__main__":
    main()

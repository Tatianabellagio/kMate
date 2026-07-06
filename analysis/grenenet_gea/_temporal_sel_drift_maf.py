#!/usr/bin/env python
"""ROBUSTNESS for the per-site SNP-vs-nonSNP selection comparison (user 2026-07-02): add explicit
(1) MAF-band control and (2) a parametric DRIFT null, on top of the primary SNP-empirical p0-matched
test (_temporal_selection_snp_vs_nonsnp.py). Primary result there: non-SNP & indel ~1.0x (behave
like matched SNPs); SVs show a directional DOWN/purged excess ~1.12x at 27/31 sites. Here we ask:
does that SV-down excess survive MAF stratification and a drift null?

Two nulls, per site:
  A. SNP-EMPIRICAL (realized drift+linked): a variant is "selected DOWN" if Δp < the 5th-percentile
     of the SNP Δp distribution in its p0 decile. SNP down-rate = 5% by construction; fold =
     category down-rate / 5%. Controls MAF (p0 bin) and drift (SNPs embody the actual founder drift).
  B. PARAMETRIC DRIFT (idealized Wright-Fisher): drift SD(p0) = sqrt(p0(1-p0)*(1-(1-1/(2N))**g)) with
     per-site census N (flowers). z = Δp/SD; "selected DOWN" if z < -1.96. We REPORT the SNP down-rate
     under this null: if >>2.5% the census-drift under-disperses (AF is a founder projection -> real
     drift is founder-level, not census) -> null B is mis-specified and null A is the right one. Then
     Ne is RE-CALIBRATED so SNP down-rate == 2.5% (that Ne ~ founder eff-N), and SV-down re-checked.

MAF band = min(p0,1-p0): rare [0.02,0.10) / mid [0.10,0.20) / common [0.20,0.50]. Reports SV, indel,
nonSNP DOWN-fold per band; SV split ins/del. Aggregate across sites (median + sign test).

Env: kmate. Writes results/grenenet_gea/sv_adaptive/temporal_sel_drift_maf.csv (+ _summary.csv).
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
MIN_P0, SV_BP, MIN_MAC, NBIN, TAILQ = 0.02, 50, 12, 10, 0.95
BANDS = [("rare", 0.02, 0.10), ("mid", 0.10, 0.20), ("common", 0.20, 0.501)]
CENSUS = pd.read_csv("results/grenenet_gea/fitness/site_census_N.csv").set_index("site")["N_per_gen"]


def down_fold(dp, p0, snp_dn, binf, band):
    """fraction of focal variants (in MAF band) below the SNP p0-matched 5th-pctile / 0.05."""
    maf = np.minimum(p0, 1 - p0)
    sel = (maf >= band[1]) & (maf < band[2])
    if sel.sum() < 20:
        return np.nan, int(sel.sum())
    d = dp[sel]; b = binf(p0[sel])
    rate = (d < snp_dn[b]).mean()
    return float(rate / (1 - TAILQ)), int(sel.sum())


def calib_Ne(dp_snp, p0_snp, g, target=0.025):
    """Ne so that the fraction of SNPs with drift-z < -1.96 equals target (one-sided down)."""
    from scipy.optimize import brentq
    def frac(N):
        sd = np.sqrt(p0_snp * (1 - p0_snp) * (1 - (1 - 1 / (2 * N)) ** g))
        return (dp_snp / sd < -1.96).mean() - target
    try:
        return float(brentq(frac, 5, 5e5))
    except Exception:
        return np.nan


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
        g = Fn.shape[0] - 1
        dp_non = np.clip(Fn[-1], 0, 1) - p0_non; dp_snp = np.clip(Fs[-1], 0, 1) - p0_snp
        b1 = float(clim["bio1"].reindex([site]).iloc[0]) if site in clim.index else np.nan
        ds = dp_snp[snp_m]; bs = binf(p0_snp[snp_m])
        snp_dn = np.array([np.quantile(ds[bs == k], 1 - TAILQ) if (bs == k).any() else -np.inf
                           for k in range(NBIN)])
        # --- Null A: SNP-empirical, per MAF band ---
        rec = dict(site=site, bio1=round(b1, 1), gens=g)
        for name, m, dp, p0 in (("nonsnp", nonsnp_m, dp_non, p0_non),
                                ("indel", ind_m, dp_non, p0_non), ("sv", sv_m, dp_non, p0_non)):
            for band in BANDS:
                f, n = down_fold(dp[m], p0[m], snp_dn, binf, band)
                rec[f"A_{name}_{band[0]}"] = round(f, 3) if f == f else np.nan
        # SV ins/del down-fold (all MAF)
        for tag, mm in (("svdel", sv_m & is_del), ("svins", sv_m & ~is_del)):
            f, n = down_fold(dp_non[mm], p0_non[mm], snp_dn, binf, ("all", 0.02, 0.501))
            rec[f"A_{tag}"] = round(f, 3) if f == f else np.nan
        # --- Null B: parametric drift (census N, then calibrated Ne) ---
        Ncen = float(CENSUS.reindex([site]).iloc[0]) if site in CENSUS.index else np.nan
        for Ntag, N in (("census", Ncen), ("calib", calib_Ne(ds, p0_snp[snp_m], g))):
            if not (N == N) or N <= 0:
                continue
            sd_snp = np.sqrt(p0_snp[snp_m] * (1 - p0_snp[snp_m]) * (1 - (1 - 1 / (2 * N)) ** g))
            snp_down = float((ds / sd_snp < -1.96).mean())
            sd_sv = np.sqrt(p0_non[sv_m] * (1 - p0_non[sv_m]) * (1 - (1 - 1 / (2 * N)) ** g))
            sv_down = float((dp_non[sv_m] / sd_sv < -1.96).mean())
            rec[f"B_{Ntag}_Ne"] = int(N)
            rec[f"B_{Ntag}_snp_down"] = round(snp_down, 4)
            rec[f"B_{Ntag}_sv_down"] = round(sv_down, 4)
            rec[f"B_{Ntag}_sv_fold"] = round(sv_down / snp_down, 3) if snp_down > 0 else np.nan
        rows.append(rec)
        print(f"  site {site:>2} b1={b1:>5.1f} g={g}: A SV-down rare/mid/com="
              f"{rec.get('A_sv_rare')}/{rec.get('A_sv_mid')}/{rec.get('A_sv_common')} "
              f"(del {rec.get('A_svdel')} ins {rec.get('A_svins')}) | "
              f"B census SNP-down={rec.get('B_census_snp_down')} SVfold={rec.get('B_census_sv_fold')} "
              f"| calib Ne={rec.get('B_calib_Ne')} SVfold={rec.get('B_calib_sv_fold')}", flush=True)

    df = pd.DataFrame(rows).sort_values("bio1")
    df.to_csv(f"{lib.GEA}/sv_adaptive/temporal_sel_drift_maf.csv", index=False)

    def sgn(col):
        c = df[col].dropna(); n = int((c > 1).sum())
        return round(c.median(), 3), f"{n}/{len(c)}", round(stats.binomtest(n, len(c), 0.5).pvalue, 4)

    print("\n=== SV DOWN/purged fold vs matched-SNP null, across sites ===")
    print("NULL A (SNP-empirical, MAF-band):")
    for band in ["rare", "mid", "common"]:
        for name in ["sv", "indel", "nonsnp"]:
            m, s, p = sgn(f"A_{name}_{band}")
            if name == "sv" or band == "rare":
                print(f"  {name:>6} {band:<7}: median {m}x  >1 at {s}  sign-p={p}")
    for tag in ["svdel", "svins"]:
        m, s, p = sgn(f"A_{tag}"); print(f"  {tag:>6} (all) : median {m}x  >1 at {s}  sign-p={p}")
    print("\nNULL B (parametric drift): SNP down-rate should be ~2.5% if the null fits; census N:")
    print(f"  census SNP down-rate median = {df['B_census_snp_down'].median():.3f} "
          f"(>>0.025 => census under-disperses, AF is founder-projection); "
          f"SV fold median = {df['B_census_sv_fold'].median():.2f}x")
    m, s, p = sgn("B_calib_sv_fold")
    print(f"  calibrated-Ne (SNP down=2.5% by design) median Ne={int(df['B_calib_Ne'].median())}; "
          f"SV-down fold median {m}x  >1 at {s}  sign-p={p}")
    df.to_csv(f"{lib.GEA}/sv_adaptive/temporal_sel_drift_maf.csv", index=False)
    print(f"\n[wrote] temporal_sel_drift_maf.csv")


if __name__ == "__main__":
    main()

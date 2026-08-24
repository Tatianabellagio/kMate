#!/usr/bin/env python
"""CLEANER SV-purging test (user 2026-07-02). Two climate-differential views that cancel the founder-
hitchhiking offset and the p0/logit artifact (both per-variant constants across sites):

 (1) PER-VARIANT CLIMATE SLOPE  beta_v = d s_v / d bio1  across all sites (s = plot-replicate logit-
     slope). Constant purging (hitchhiking) -> beta=0; beta<0 only if purging INTENSIFIES with heat.
     Compare beta distribution SV vs frequency-matched SNP -> climate-driven SV purging, confounds
     removed. Saved: beta + mean_s + p0 per variant (all SV/indel; 25k SNP subsample).

 (2) SIGN ENRICHMENT vs bio1 (user's idea): per site, frequency-matched fraction with s<0 (purged)
     for SV / indel / matched-SNP; regress the SV-SNP purging excess on bio1. Saved per site.

Env: kmate. Writes analysis/grenenet_gea/sv_adaptive/s_climate_slope.npz + _sign_by_site.csv .
"""
import os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import importlib.util
_sp = importlib.util.spec_from_file_location(
    "pm", os.path.join(os.path.dirname(os.path.abspath(__file__)), "_temporal_s_plots_snp_vs_nonsnp.py"))
pm = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(pm)

os.chdir("/global/scratch/users/tbellg/kmate")
STORE = lib.AF_STORE; PM = f"{lib.GEA}/pool_matrices"
MIN_P0, SV_BP, MIN_MAC, NBIN, NDRAW = 0.02, 50, 12, 25, 20000
rng = np.random.default_rng(0)


def match_to(tp0, ss, sp0):
    e = np.quantile(np.concatenate([tp0, sp0]), np.linspace(0, 1, NBIN + 1)); e[-1] += 1e-9
    tb = np.digitize(tp0, e[1:-1]); sb = np.digitize(sp0, e[1:-1]); out = []
    for b in range(NBIN):
        k = int(round((tb == b).mean() * NDRAW)); pool = ss[sb == b]
        if k > 0 and pool.size: out.append(rng.choice(pool, k, replace=True))
    return np.concatenate(out)


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
    dl = dlen[cols_non]; sv_k = dl > SV_BP; ind_k = (dl >= 1) & (dl <= SV_BP); isdel = is_del[cols_non]

    # per-site climate (bio1 temp, bio18 dry-summer precip) from pool meta
    cmeta = pd.concat([pd.read_csv(f"{PM}/pool_gen{g}_snp.meta.csv") for g in (1, 2, 3)],
                      ignore_index=True).groupby("site")[["bio1", "bio18"]].mean()
    CVARS = ["bio1", "bio18"]
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                    for s in pd.read_csv(f).site.unique()})

    # accumulators for per-variant climate slope beta_k = cov(s,clim_k)/var(clim_k), per climate var
    ns = np.zeros(cols_snp.size); ssum = np.zeros(cols_snp.size)
    nn = np.zeros(cols_non.size); nsum = np.zeros(cols_non.size)
    ssb = {k: np.zeros(cols_snp.size) for k in CVARS}   # sum s*clim per snp
    nsb = {k: np.zeros(cols_non.size) for k in CVARS}
    B = {k: 0.0 for k in CVARS}; B2 = {k: 0.0 for k in CVARS}; nS = 0.0
    sign = []
    for site in sites:
        s_sn, _, _, npl = pm.plot_replicate_sz(site, "snp", cols_snp, p0s)
        s_no, _, _, _ = pm.plot_replicate_sz(site, "nonsnp", cols_non, p0n)
        if npl < 2 or site not in cmeta.index:
            continue
        cv = {k: float(cmeta.loc[site, k]) for k in CVARS}
        s_sn = s_sn.astype(np.float64); s_no = s_no.astype(np.float64)
        ssum += s_sn; nsum += s_no; ns += 1; nn += 1; nS += 1
        for k in CVARS:
            ssb[k] += s_sn * cv[k]; nsb[k] += s_no * cv[k]; B[k] += cv[k]; B2[k] += cv[k] ** 2
        p_sv = p0n[sv_k]
        msn = match_to(p_sv, s_sn, p0s); mind = match_to(p_sv, s_no[ind_k], p0n[ind_k])
        sign.append(dict(site=site, bio1=round(cv["bio1"], 1), bio18=round(cv["bio18"], 1),
                         sv_neg=float((s_no[sv_k] < 0).mean()),
                         msnp_neg=float((msn < 0).mean()),
                         mindel_neg=float((mind < 0).mean())))
        print(f"  site {site:>2} bio1={cv['bio1']:>5.1f} bio18={cv['bio18']:>5.0f} done", flush=True)

    def beta(sumsc, ssumv, cnt, k):
        return (nS * sumsc[k] - B[k] * ssumv) / (nS * B2[k] - B[k] ** 2)
    i_sn = rng.choice(cols_snp.size, 25000, replace=False)
    out = dict(p0_snp=p0s[i_sn].astype(np.float32), mean_snp=(ssum / ns)[i_sn].astype(np.float32),
               p0_sv=p0n[sv_k].astype(np.float32), mean_sv=(nsum / nn)[sv_k].astype(np.float32),
               sv_isdel=isdel[sv_k], p0_indel=p0n[ind_k].astype(np.float32),
               mean_indel=(nsum / nn)[ind_k].astype(np.float32))
    for k in CVARS:
        bsnp = beta(ssb, ssum, ns, k); bnon = beta(nsb, nsum, nn, k)
        out[f"beta_snp_{k}"] = bsnp[i_sn].astype(np.float32)
        out[f"beta_sv_{k}"] = bnon[sv_k].astype(np.float32)
        out[f"beta_indel_{k}"] = bnon[ind_k].astype(np.float32)
    np.savez_compressed(f"{lib.GEA}/sv_adaptive/s_climate_slope.npz", **out)
    pd.DataFrame(sign).to_csv(f"{lib.GEA}/sv_adaptive/s_climate_slope_sign_by_site.csv", index=False)
    # quick text summary (uses the saved, subsampled arrays)
    from scipy import stats
    for k in CVARS:
        bsv = out[f"beta_sv_{k}"]; bmatch = match_to(out["p0_sv"], out[f"beta_snp_{k}"], out["p0_snp"])
        print(f"\n=== climate slope beta = d s / d {k} ===")
        print(f"  median beta: SNP {np.median(out[f'beta_snp_{k}']):+.5f}  "
              f"indel {np.median(out[f'beta_indel_{k}']):+.5f}  SV {np.median(bsv):+.5f}")
        print(f"  SV vs matched-SNP KS p={stats.ks_2samp(bsv,bmatch).pvalue:.2e}; "
              f"frac beta<0 SV={np.mean(bsv<0):.3f} mSNP={np.mean(bmatch<0):.3f}")
        ins = bsv[~out['sv_isdel']]; dele = bsv[out['sv_isdel']]
        print(f"  insertions median {np.median(ins):+.5f} (KS p={stats.ks_2samp(ins,bmatch).pvalue:.1e}) | "
              f"deletions median {np.median(dele):+.5f} (KS p={stats.ks_2samp(dele,bmatch).pvalue:.1e})")
    print(f"\n[wrote] s_climate_slope.npz + _sign_by_site.csv")


if __name__ == "__main__":
    main()

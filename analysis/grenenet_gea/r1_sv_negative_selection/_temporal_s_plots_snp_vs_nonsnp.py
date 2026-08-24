#!/usr/bin/env python
"""SELECTION COEFFICIENT with PLOTS AS REPLICATES, per variant per site (user 2026-07-02). Gains
accuracy + a real test by using the ~10-12 replicate plots per site instead of the plot-pooled
trajectory.

Per variant v, per site:
  - each plot j (present at the site) gives a trajectory [p0 (shared founding, gen0), AF_{j,g1},
    AF_{j,g2}, ...]; s_{j,v} = OLS slope of logit(freq) on generation = per-plot selection coef.
  - s_v  = mean_j s_{j,v}                    (site selection coefficient)
    se_v = sd_j(s_{j,v}) / sqrt(n_plots)     (AMONG-PLOT SE = empirical drift floor)
    z_v  = s_v / se_v                         (parallel across plots => selection; scattered => drift)
This z is the accuracy gain: it separates real selection (same direction in independent replicate
plots) from drift (plots scatter around 0) WITHOUT any parametric drift model.

Question: are NON-SNPs (indels+SVs) under selection (significant, parallel s) more than SNPs?
  - Absolute: fraction with z < -2 (purged) or z > 2 (favoured), per class.
  - Frequency-matched robustness: SNP z-distribution per p0-decile defines the tail; category fold =
    (category tail rate)/(SNP 5% tail). MAF bands (rare/mid/common) + SV ins/del.
Aggregate over sites (median + sign test).

CAVEAT (unchanged): pool AF is a GLOBAL-mode founder projection, so s_v is a projection of founder-h
slopes and a SNP on the same founders shares it -> cannot separate SV-specific selection from riding
selected haplotypes; needs independent local-mode/vg SV AF. (But plot replicates DO remove the
drift-vs-selection ambiguity, which Δp could not.)

Env: kmate. Writes analysis/grenenet_gea/sv_adaptive/temporal_s_plots_snp_vs_nonsnp.csv (+_summary).
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

os.chdir("/global/scratch/users/tbellg/kmate")
STORE = lib.AF_STORE
PM = f"{lib.GEA}/pool_matrices"
MIN_P0, SV_BP, MIN_MAC, NBIN, TAILQ, EPS = 0.02, 50, 12, 10, 0.95, 1e-3
BANDS = [("rare", 0.02, 0.10), ("mid", 0.10, 0.20), ("common", 0.20, 0.501)]


def plot_replicate_sz(site, kind, cols, p0c):
    """s, se, z per variant (over `cols`) from plot-replicate logit-slopes. gen0 = shared p0c."""
    L0 = np.log(np.clip(p0c, EPS, 1 - EPS) / (1 - np.clip(p0c, EPS, 1 - EPS)))
    plot_pts = {}                                            # plot -> list of (t, logit_af_over_cols)
    for g in (1, 2, 3):
        try:
            meta = pd.read_csv(f"{PM}/pool_gen{g}_{kind}.meta.csv")
        except FileNotFoundError:
            continue
        m4 = meta[meta.site == site]
        if len(m4) == 0:
            continue
        mat = np.load(f"{PM}/pool_gen{g}_{kind}_af.npy", mmap_mode="r")
        for ridx, plot in zip(m4.index.to_numpy(), m4["plot"].to_numpy()):
            af = np.asarray(mat[ridx])[cols].astype(np.float64)
            plot_pts.setdefault(int(plot), []).append(
                (float(g), np.log(np.clip(af, EPS, 1 - EPS) / (1 - np.clip(af, EPS, 1 - EPS)))))
    slopes = []
    for plot, pts in plot_pts.items():
        ts = np.array([0.0] + [t for t, _ in pts])
        if ts.size < 2:
            continue
        Y = np.vstack([L0] + [y for _, y in pts])            # [nt x ncols]
        tc = ts - ts.mean(); sst = float((tc ** 2).sum())
        slopes.append((tc @ Y) / sst)                        # [ncols]
    S = np.vstack(slopes)                                    # [nplots x ncols]
    n = S.shape[0]
    s = S.mean(0)
    se = S.std(0, ddof=1) / np.sqrt(n) if n > 1 else np.full(S.shape[1], np.nan)
    z = np.divide(s, se, out=np.zeros_like(s), where=se > 0)
    return s.astype(np.float32), se.astype(np.float32), z.astype(np.float32), n


def tailfold(z, p0, snp_dn, snp_up, binf, band):
    maf = np.minimum(p0, 1 - p0)
    sel = (maf >= band[1]) & (maf < band[2])
    if sel.sum() < 20:
        return np.nan, np.nan, int(sel.sum())
    zv = z[sel]; b = binf(p0[sel])
    dn = (zv < snp_dn[b]).mean() / (1 - TAILQ)
    up = (zv > snp_up[b]).mean() / (1 - TAILQ)
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
    snp_m = common_snp & reach_snp
    non_keep = (dlen >= 1) & common_non & reach_non          # nonSNP columns we score
    cols_snp = np.where(snp_m)[0]; cols_non = np.where(non_keep)[0]
    p0_snp_k = p0_snp[cols_snp]; p0_non_k = p0_non[cols_non]
    dlen_k = dlen[cols_non]; isdel_k = is_del[cols_non]
    sv_k = dlen_k > SV_BP; ind_k = (dlen_k >= 1) & (dlen_k <= SV_BP)
    print(f"scored cols: SNP={len(cols_snp):,} nonSNP={len(cols_non):,} "
          f"(indel={int(ind_k.sum()):,} SV={int(sv_k.sum()):,})", flush=True)

    clim = lib.load_climate()
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                    for s in pd.read_csv(f).site.unique()})
    p0q = np.quantile(p0_snp_k, np.linspace(0, 1, NBIN + 1)); p0q[-1] += 1e-6
    binf = lambda p: np.clip(np.digitize(p, p0q[1:-1]), 0, NBIN - 1)

    rows = []
    for site in sites:
        s_sn, se_sn, z_sn, n_sn = plot_replicate_sz(site, "snp", cols_snp, p0_snp_k)
        s_no, se_no, z_no, n_no = plot_replicate_sz(site, "nonsnp", cols_non, p0_non_k)
        if n_sn < 2:
            print(f"  site {site}: {n_sn} plots, skip"); continue
        b1 = float(clim["bio1"].reindex([site]).iloc[0]) if site in clim.index else np.nan
        bs = binf(p0_snp_k)
        snp_dn = np.array([np.quantile(z_sn[bs == k], 1 - TAILQ) if (bs == k).any() else -np.inf
                           for k in range(NBIN)])
        snp_up = np.array([np.quantile(z_sn[bs == k], TAILQ) if (bs == k).any() else np.inf
                           for k in range(NBIN)])
        rec = dict(site=site, bio1=round(b1, 1), n_plots=n_sn,
                   med_s_snp=round(float(np.median(s_sn)), 4),
                   med_s_indel=round(float(np.median(s_no[ind_k])), 4),
                   med_s_sv=round(float(np.median(s_no[sv_k])), 4),
                   # absolute parallel-selection rates (z-based), per class
                   snp_dn2=round(float((z_sn < -2).mean()), 4),
                   snp_up2=round(float((z_sn > 2).mean()), 4),
                   ind_dn2=round(float((z_no[ind_k] < -2).mean()), 4),
                   ind_up2=round(float((z_no[ind_k] > 2).mean()), 4),
                   sv_dn2=round(float((z_no[sv_k] < -2).mean()), 4),
                   sv_up2=round(float((z_no[sv_k] > 2).mean()), 4))
        for name, msk in (("nonsnp", non_keep[cols_non] | True), ("indel", ind_k), ("sv", sv_k)):
            zc = z_no if name == "nonsnp" else z_no[msk]
            pc = p0_non_k if name == "nonsnp" else p0_non_k[msk]
            for band in BANDS:
                dn, up, n = tailfold(zc, pc, snp_dn, snp_up, binf, band)
                rec[f"{name}_{band[0]}_dn"] = round(dn, 3) if dn == dn else np.nan
                rec[f"{name}_{band[0]}_up"] = round(up, 3) if up == up else np.nan
        for tag, mm in (("svdel", sv_k & isdel_k), ("svins", sv_k & ~isdel_k)):
            dn, up, n = tailfold(z_no[mm], p0_non_k[mm], snp_dn, snp_up, binf, ("all", 0.02, 0.501))
            rec[f"{tag}_dn"] = round(dn, 3) if dn == dn else np.nan
        rows.append(rec)
        print(f"  site {site:>2} b1={b1:>5.1f} plots={n_sn}: med_s SNP={rec['med_s_snp']:+.3f} "
              f"SV={rec['med_s_sv']:+.3f} | z<-2 SNP={100*rec['snp_dn2']:.1f}% SV={100*rec['sv_dn2']:.1f}% "
              f"| z>2 SNP={100*rec['snp_up2']:.1f}% SV={100*rec['sv_up2']:.1f}% | SVdown fold mid/com="
              f"{rec.get('sv_mid_dn')}/{rec.get('sv_common_dn')}", flush=True)

    df = pd.DataFrame(rows).sort_values("bio1")
    df.to_csv(f"{lib.GEA}/sv_adaptive/temporal_s_plots_snp_vs_nonsnp.csv", index=False)

    def sgn(col, ref=1.0):
        c = df[col].dropna(); n = int((c > ref).sum())
        return round(c.median(), 3), f"{n}/{len(c)}", round(stats.binomtest(n, len(c), 0.5).pvalue, 4)

    print("\n=== PLOT-REPLICATE s: is non-SNP under selection more than SNPs? ===")
    print(f"median s/gen: SNP {df.med_s_snp.median():+.4f} | indel {df.med_s_indel.median():+.4f} | "
          f"SV {df.med_s_sv.median():+.4f}")
    print("\nABSOLUTE parallel-selection rate (fraction with z beyond +/-2 across plots):")
    for d, lab in (("dn2", "purged z<-2"), ("up2", "favoured z>2")):
        print(f"  {lab:<14}: SNP {100*df[f'snp_{d}'].median():.1f}%  indel {100*df[f'ind_{d}'].median():.1f}%"
              f"  SV {100*df[f'sv_{d}'].median():.1f}%")
    print("\nFREQUENCY-MATCHED fold vs SNP z-tail (DOWN=purged, UP=favoured):")
    for name in ["nonsnp", "indel", "sv"]:
        for band in ["rare", "mid", "common"]:
            md, sd, pd_ = sgn(f"{name}_{band}_dn"); mu, su, pu = sgn(f"{name}_{band}_up")
            print(f"  {name:>6} {band:<7}: DOWN {md}x ({sd}, p={pd_}) | UP {mu}x ({su}, p={pu})")
    for tag in ["svdel", "svins"]:
        m, s, p = sgn(f"{tag}_dn"); print(f"  {tag:>6} DOWN(all): {m}x ({s}, p={p})")
    print(f"\n[wrote] temporal_s_plots_snp_vs_nonsnp.csv")


if __name__ == "__main__":
    main()

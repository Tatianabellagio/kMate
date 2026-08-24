#!/usr/bin/env python
"""Two-stage linkage-aware per-block selection: TEMPORAL slope + per-plot SE -> precision-
weighted LD mixed model (default site 4). Upgrades block_ld_lmm.py (which used the
founding->gen3 endpoint mean) to use the full gen-0..3 trajectory AND the replicate structure.

STAGE 1 (per block b, temporal + replicate):
  per plot j: logit-freq trajectory over t=0(founding),1,2,3; per-plot OLS slope beta_{b,j}.
  s_b  = mean_j beta_{b,j}           (selection coefficient, full slope)
  SE_b = sd_j(beta_{b,j})/sqrt(n)    (replicate measurement error)
Per-plot info enters as the EFFECT (slope) and a precision WEIGHT (SE) -- NOT as the null
(that over-precision broke the random-slope Wald). Calibration comes from the LD covariance.

STAGE 2 (across blocks, precision-weighted LD-LMM):
  s = mu*1 + g + e,   g ~ N(0, sg^2 C_LD),   e_b ~ N(0, se^2 + SE_b^2)
  C_LD = haplotype relationship (founder-sharing = LD). g = R u, R = A'/sqrt(nF), u~N(0,sg^2 I).
  V = sg^2 C_LD + diag(se^2 + SE_b^2). REML over (sg^2, se^2) via Woodbury (low-rank, nF=231).
  block-specific residual r_b = s_b - BLUP_b;  z_b = r_b / sqrt(se^2 + SE_b^2).

Writes site<ID>_block_ld_lmm_temporal.csv + _meta.json + manhattan. Env: kmate. SITE via env.
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats, optimize
from scipy.linalg import cho_factor, cho_solve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import build_genotype

H = os.environ.get("HF_DIR", "results/grenenet_gea/hapfreq")  # e.g. .../hapfreq_clq90
SITE = int(os.environ.get("SITE", 4))
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
EPS = 1e-3
Tg = np.array([0.0, 1.0, 2.0, 3.0])


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def main():
    G, founders, reg = build_genotype()
    mat = np.load(f"{H}/hapfreq_matrix.npy")
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    p0 = np.load(f"{H}/hapfreq_p0_seedmix.npy")
    fr = pd.read_csv(f"{H}/hapfreq_registry.csv")
    nhap = mat.shape[1]

    # ---- STAGE 1: per-plot logit-freq trajectory -> per-block slope s_b and SE_b ----
    pt = lib.pool_table()
    s = pt[(pt.site == SITE) & pt.sampleid.astype(str).isin(sidx)].copy()
    # flower-weighted hapfreq per (gen,plot)
    cell = {}
    for (gen, plot), g in s.groupby(["generation", "plot"]):
        rows = g.sampleid.astype(str).map(sidx).to_numpy()
        w = g.flowerscollected.to_numpy(float); w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        cell[(int(gen), int(plot))] = (mat[rows] * w[:, None]).sum(0) / w.sum()
    present = {}
    for (gen, plot) in cell:
        present.setdefault(plot, set()).add(gen)
    plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
    n = len(plots)
    # per-plot slope over t=0..3 (t0 = founding p0, shared); Sigma(t-1.5)^2 = 5
    coef = (Tg - 1.5)
    slopes = np.zeros((n, nhap))
    for j, pl in enumerate(plots):
        y = np.vstack([logit(p0), logit(cell[(1, pl)]), logit(cell[(2, pl)]), logit(cell[(3, pl)])])  # (4,nhap)
        slopes[j] = (coef[:, None] * y).sum(0) / 5.0
    s_all = slopes.mean(0)
    se_all = slopes.std(0, ddof=1) / np.sqrt(n)

    # testable haplotypes, k-1 per block
    fr["unit"] = fr.chrom + ":" + fr.unit_start.astype(str) + "-" + fr.unit_end.astype(str)
    testable = (fr.covered & fr.panel_freq.between(0.05, 0.95)).to_numpy()
    keep = np.zeros(nhap, bool)
    for u, grp in fr[testable].groupby("unit"):
        kept = grp.index.drop(grp.panel_freq.idxmax()) if len(grp) > 1 else grp.index
        keep[kept.to_numpy()] = True
    ki = np.where(keep)[0]; M = len(ki)
    sb = s_all[ki].copy(); sb -= sb.mean()
    D0 = se_all[ki] ** 2                                       # known per-block measurement var
    D0 = np.clip(D0, np.median(D0) * 1e-3, None)               # floor tiny SEs

    # ---- STAGE 2: precision-weighted LD-LMM (Woodbury, low-rank) ----
    A = G[:, ki].astype(np.float64)
    pf = A.mean(0); A = (A - pf) / np.sqrt(pf * (1 - pf) + 1e-9)   # (nF, M)
    nF = A.shape[0]
    one = np.ones(M)

    def pieces(tau, se2):
        Winv = 1.0 / (se2 + D0)                                # (M,)
        AWA = (A * Winv[None, :]) @ A.T / nF                   # (nF,nF) = R'WinvR
        Inner = np.eye(nF) / tau + AWA
        cf = cho_factor(Inner, lower=True)

        def Vinv(x):
            Wx = Winv * x
            t = (A @ Wx) / np.sqrt(nF)
            u = cho_solve(cf, t)
            return Wx - Winv * ((A.T @ u) / np.sqrt(nF))
        logdetV = np.sum(np.log(se2 + D0)) + np.linalg.slogdet(np.eye(nF) + tau * AWA)[1]
        return Vinv, logdetV

    def neg_reml(par):
        tau, se2 = np.exp(par)
        Vinv, logdetV = pieces(tau, se2)
        Vi1 = Vinv(one); Vis = Vinv(sb)
        mu = (one @ Vis) / (one @ Vi1)
        r = sb - mu; quad = r @ Vinv(r)
        return 0.5 * (logdetV + quad + np.log(one @ Vi1))      # REML (X=intercept)

    res = optimize.minimize(neg_reml, x0=np.log([np.var(sb), np.var(sb) / 2]),
                            method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-3})
    tau, se2 = np.exp(res.x)
    Vinv, _ = pieces(tau, se2)
    Vi1 = Vinv(one); mu = (one @ Vinv(sb)) / (one @ Vi1)
    r = sb - mu
    Vir = Vinv(r)
    g_b = (tau / nF) * (A.T @ (A @ Vir))                       # BLUP linked background
    resid = r - g_b
    z = resid / np.sqrt(se2 + D0)
    pval = 2 * stats.norm.sf(np.abs(z))
    lam_gc = np.median(stats.chi2.isf(np.clip(pval, 1e-300, 1), 1)) / stats.chi2.ppf(0.5, 1)
    h2 = tau / (tau + se2)

    d = fr.loc[ki, ["chrom", "unit_start", "unit_end", "panel_freq"]].copy()
    d["s"] = sb; d["se"] = np.sqrt(D0); d["linked_bg"] = g_b; d["resid"] = resid
    d["z"] = z; d["p"] = pval
    m = len(d); o = d.p.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate((d.p.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q"] = np.clip(q, 0, 1)
    bonf = 0.05 / m
    d.sort_values("p").to_csv(f"{H}/site{SITE}_block_ld_lmm_temporal.csv", index=False)
    meta = dict(site=SITE, n_plot=n, M=int(M), nF=int(nF), h2=float(h2),
                sigma_g2=float(tau), sigma_e2=float(se2), lambda_gc=float(lam_gc),
                n_bonf=int((d.p < bonf).sum()), n_fdr05=int((d.q < 0.05).sum()),
                model="two-stage temporal slope + per-plot SE, precision-weighted LD-LMM")
    json.dump(meta, open(f"{H}/site{SITE}_block_ld_lmm_temporal_meta.json", "w"), indent=2)

    print(f"site {SITE}: {n} plots, {M:,} haploblocks (TEMPORAL slope + per-plot SE)")
    print(f"  median SE_b = {np.median(np.sqrt(D0)):.3f}  (replicate measurement error)")
    print(f"  REML h2(linkage) = {h2:.2f}, sigma_e^2 = {se2:.3f}")
    print(f"  residual lambda_GC = {lam_gc:.2f}")
    print(f"  block-specific hits: Bonferroni {int((d.p<bonf).sum())} | FDR q<0.05 {int((d.q<0.05).sum())}")
    print("\n  TOP block-specific blocks:")
    print(d.sort_values("p").head(12)[
        ["chrom", "unit_start", "unit_end", "panel_freq", "s", "se", "linked_bg", "resid", "z", "p", "q"]].to_string(index=False))

    g2 = d.copy(); g2["mid"] = (g2.unit_start + g2.unit_end) / 2
    g2 = g2.sort_values(["chrom", "unit_start"]).reset_index(drop=True)
    off, centers, x = 0.0, [], np.zeros(len(g2))
    for ch in CHROMS:
        mk = (g2.chrom == ch).to_numpy()
        if not mk.any():
            continue
        x[mk] = g2.mid[mk] + off; centers.append(off + g2.mid[mk].max() / 2); off += g2.mid[mk].max() * 1.02
    g2["x"] = x
    fig, ax = plt.subplots(figsize=(11, 4.2))
    nlp = -np.log10(np.clip(g2.p, 1e-300, 1))
    for i, ch in enumerate(CHROMS):
        mk = (g2.chrom == ch).to_numpy()
        ax.scatter(g2.x[mk], nlp[mk], s=5, c=["#3b4cc0", "#7aa0c4"][i % 2], alpha=.5, edgecolors="none", rasterized=True)
    ax.axhline(-np.log10(bonf), color="firebrick", lw=.9, ls="--", label="Bonferroni")
    ax.set_xticks(centers); ax.set_xticklabels(CHROMS)
    ax.set_ylabel(r"$-\log_{10}p$"); ax.set_xlabel("genome position")
    ax.set_title(f"Site {SITE} block-specific selection — temporal slope + per-plot SE, LD-LMM "
                 f"(λ_GC={lam_gc:.2f}, h²_LD={h2:.2f})", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{H}/site{SITE}_block_ld_lmm_temporal.png", dpi=150)
    print(f"\n[done] {H}/site{SITE}_block_ld_lmm_temporal.csv + _meta.json + manhattan")


if __name__ == "__main__":
    main()

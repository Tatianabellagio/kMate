#!/usr/bin/env python
"""Sampling-variance-weighted (beta-binomial-style) temporal LD-LMM (default site 4).

Upgrades block_ld_lmm_temporal.py. The per-block measurement variance D0 in the residual
term e_b ~ N(0, sigma_e^2 + D0_b) was the replicate SE alone, which COLLAPSES for rare /
near-fixed haplotypes (all plots agree near the logit clip -> tiny spread -> huge z ->
inflated tail; 100% of the old Bonferroni hits were panel_freq<=0.10 or >=0.90). Here D0 is
floored at the BINOMIAL SAMPLING variance of the temporal slope implied by the flower census
N_gam = 2*flowerscollected. On the logit scale Var(logit f) ~ 1/(N_gam f(1-f)), so rare/near-
fixed haplotypes get the wide uncertainty they actually have and stop driving the tail --
WITHOUT being deleted by a MAF filter (keeps the biologically interesting rare-sweep case).
sigma_e^2 (REML) absorbs residual drift overdispersion on top -> beta-binomial-style structure.

  D0_b = (coef0/5)^2 * v0_b            (founding p0 sampling var; SHARED across plots, NOT averaged)
       + max( SE_rep_b^2 ,  (1/n^2) sum_j sum_{t>=1} (coef_t/5)^2 v_{t,j,b} )   (evolved: empirical
         replicate spread OR the binomial floor, whichever is larger)
  v_{.,b} = 1 / (N_gam * f(1-f))  on the logit scale.

Founding error enters as a COMMON term (it cancels in the replicate SD, so it must be added
explicitly). Stage 2 (precision-weighted low-rank REML LD-LMM) is byte-identical to the parent.
Writes site<ID>_block_ld_lmm_sampvar.csv + _meta.json + manhattan. Env: kmate. SITE via env.
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

H = "results/grenenet_gea/hapfreq"
SITE = int(os.environ.get("SITE", 4))
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
EPS = 1e-3
Tg = np.array([0.0, 1.0, 2.0, 3.0])


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def logit_var(f, Ngam):
    """Delta-method binomial sampling variance of logit(f): 1/(N f(1-f))."""
    f = np.clip(f, EPS, 1 - EPS)
    return 1.0 / (np.clip(Ngam, 1.0, None) * f * (1 - f))


def main():
    G, founders, reg = build_genotype()
    mat = np.load(f"{H}/hapfreq_matrix.npy")
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    p0 = np.load(f"{H}/hapfreq_p0_seedmix.npy")
    fr = pd.read_csv(f"{H}/hapfreq_registry.csv")
    nhap = mat.shape[1]

    # ---- STAGE 1: per-plot logit-freq trajectory -> slope s_b, replicate SE_b, sampling var ----
    pt = lib.pool_table()
    s = pt[(pt.site == SITE) & pt.sampleid.astype(str).isin(sidx)].copy()
    cell, cellN = {}, {}
    for (gen, plot), g in s.groupby(["generation", "plot"]):
        rows = g.sampleid.astype(str).map(sidx).to_numpy()
        w = g.flowerscollected.to_numpy(float); w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        cell[(int(gen), int(plot))] = (mat[rows] * w[:, None]).sum(0) / w.sum()
        cellN[(int(gen), int(plot))] = float(w.sum())           # total flowers contributing to the cell
    present = {}
    for (gen, plot) in cell:
        present.setdefault(plot, set()).add(gen)
    plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
    n = len(plots)
    coef = (Tg - 1.5)                                           # [-1.5,-0.5,0.5,1.5], sum sq = 5

    Nmed = np.median([cellN[(g, pl)] for (g, pl) in cellN])     # founding effective-flowers proxy
    v0 = logit_var(p0, 2 * Nmed)                                # founding p0 sampling var (shared)
    v0_term = (coef[0] / 5.0) ** 2 * v0                         # founding COMMON term (not averaged)

    slopes = np.zeros((n, nhap))
    evolved_samp = np.zeros(nhap)                              # sum_j sum_{t>=1} (coef_t/5)^2 v_{t,j}
    for j, pl in enumerate(plots):
        cells = [cell[(1, pl)], cell[(2, pl)], cell[(3, pl)]]
        Ns = [cellN[(1, pl)], cellN[(2, pl)], cellN[(3, pl)]]
        y = np.vstack([logit(p0), logit(cells[0]), logit(cells[1]), logit(cells[2])])  # (4,nhap)
        slopes[j] = (coef[:, None] * y).sum(0) / 5.0
        for ti, (c, Nc) in enumerate(zip(cells, Ns), start=1):
            evolved_samp += (coef[ti] / 5.0) ** 2 * logit_var(c, 2 * Nc)
    evolved_samp /= n ** 2                                      # variance of the plot-mean (indep across plots)
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
    rep = se_all[ki] ** 2                                       # empirical replicate variance
    samp_evo = evolved_samp[ki]                                 # binomial floor (evolved cells)
    D0 = v0_term[ki] + np.maximum(rep, samp_evo)                # 2-component point variance
    D0 = np.clip(D0, np.median(D0) * 1e-3, None)

    # ---- STAGE 2: precision-weighted LD-LMM (Woodbury, low-rank) -- identical to parent ----
    A = G[:, ki].astype(np.float64)
    pf = A.mean(0); A = (A - pf) / np.sqrt(pf * (1 - pf) + 1e-9)   # (nF, M)
    nF = A.shape[0]
    one = np.ones(M)

    def pieces(tau, se2):
        Winv = 1.0 / (se2 + D0)
        AWA = (A * Winv[None, :]) @ A.T / nF
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
        return 0.5 * (logdetV + quad + np.log(one @ Vi1))

    res = optimize.minimize(neg_reml, x0=np.log([np.var(sb), np.var(sb) / 2]),
                            method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-3})
    tau, se2 = np.exp(res.x)
    Vinv, _ = pieces(tau, se2)
    Vi1 = Vinv(one); mu = (one @ Vinv(sb)) / (one @ Vi1)
    r = sb - mu
    Vir = Vinv(r)
    g_b = (tau / nF) * (A.T @ (A @ Vir))
    resid = r - g_b
    z = resid / np.sqrt(se2 + D0)
    pval = 2 * stats.norm.sf(np.abs(z))
    chi = stats.chi2.isf(np.clip(pval, 1e-300, 1), 1)
    lam_q = {q: float(np.quantile(chi, q) / stats.chi2.ppf(q, 1)) for q in [0.25, 0.5, 0.75, 0.9, 0.99, 0.999]}
    lam_gc = lam_q[0.5]
    h2 = tau / (tau + se2)

    d = fr.loc[ki, ["chrom", "unit_start", "unit_end", "panel_freq"]].copy()
    d["s"] = sb; d["se"] = np.sqrt(D0); d["linked_bg"] = g_b; d["resid"] = resid
    d["z"] = z; d["p"] = pval
    m = len(d); o = d.p.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate((d.p.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q"] = np.clip(q, 0, 1)
    bonf = 0.05 / m
    d.sort_values("p").to_csv(f"{H}/site{SITE}_block_ld_lmm_sampvar.csv", index=False)
    meta = dict(site=SITE, n_plot=n, M=int(M), nF=int(nF), h2=float(h2),
                sigma_g2=float(tau), sigma_e2=float(se2), lambda_gc=float(lam_gc),
                lambda_quantiles={str(k): v for k, v in lam_q.items()},
                median_D0=float(np.median(D0)), median_rep=float(np.median(rep)),
                median_samp_evo=float(np.median(samp_evo)), Nmed_flowers=float(Nmed),
                n_bonf=int((d.p < bonf).sum()), n_fdr05=int((d.q < 0.05).sum()),
                model="temporal slope + flower-census binomial sampling-variance floor, precision-weighted LD-LMM")
    json.dump(meta, open(f"{H}/site{SITE}_block_ld_lmm_sampvar_meta.json", "w"), indent=2)

    sig = d[d.q < 0.05]
    frac_bound = float(((sig.panel_freq <= 0.10) | (sig.panel_freq >= 0.90)).mean()) if len(sig) else float("nan")
    print(f"site {SITE}: {n} plots, {M:,} haploblocks (TEMPORAL slope + flower-census sampling-var floor)")
    print(f"  median flowers/cell = {Nmed:.0f}  ->  N_gam = {2*Nmed:.0f}")
    print(f"  median D0 = {np.median(D0):.4f}  (replicate-only was {np.median(rep):.4f}; binomial floor {np.median(samp_evo):.4f})")
    print(f"  REML h2(linkage) = {h2:.2f}, sigma_e^2 = {se2:.4f}")
    print("  lambda by quantile:  " + "  ".join(f"q{int(k*1000)/1000}={v:.2f}" for k, v in lam_q.items()))
    print(f"  block-specific hits: Bonferroni {int((d.p<bonf).sum())} | FDR q<0.05 {int((d.q<0.05).sum())}"
          f"   (frac of FDR hits at freq boundary = {frac_bound:.2f})")
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
    ax.set_title(f"Site {SITE} block-specific selection — sampling-var-weighted LD-LMM "
                 f"(λ_GC={lam_gc:.2f}, h²_LD={h2:.2f})", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{H}/site{SITE}_block_ld_lmm_sampvar.png", dpi=150)
    print(f"\n[done] {H}/site{SITE}_block_ld_lmm_sampvar.csv + _meta.json + manhattan")


if __name__ == "__main__":
    main()

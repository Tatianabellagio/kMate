#!/usr/bin/env python
"""LINEAR-scale block LD-LMM + binomial sampling-variance floor (default site 4).

Combines the two fixes the audit converged on:
  * LINEAR scale (block_ld_lmm_linear.py): f_b=Σ_{f∈b}h_f so Δf_b=ΣΔh_f is EXACTLY linear ->
    C_LD absorbs the founder/clade winning cleanly, NO logit-curvature frequency-sign artifact
    (the logit model pushed rare haps to spurious NEGATIVE residuals, near-fixed to positive).
  * SAMPLING-VARIANCE FLOOR (block_ld_lmm_temporal_sampvar.py, here on the LINEAR scale):
    the replicate SE collapses for rare/near-fixed haps (median rep/floor ~0.5 -> spurious huge
    |z|); floor D0 at the flower-census BINOMIAL variance Var(f)=f(1-f)/N_gam, N_gam=2*flowers.
  D0_b = (c0/5)^2 v0_b + max( SE_rep_b^2, (1/n^2) Σ_j Σ_{t≥1} (c_t/5)^2 f_{t,j}(1-f_{t,j})/N_gam )

The block-specific residual r_b = s_b - BLUP_b is now interpretable on the freq scale with the
correct SIGN: r_b > 0 = the block rose MORE than founder-linkage predicts (= block-specific /
convergent selection beyond hitchhiking, the TARGET); r_b < 0 = rose less (decoupling).
Writes site<ID>_block_ld_lmm_linear_sampvar.csv + _meta.json. Env: kmate. SITE via env.
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats, optimize
from scipy.linalg import cho_factor, cho_solve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import build_genotype

H = "results/grenenet_gea/hapfreq"
SITE = int(os.environ.get("SITE", 4))
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
Tg = np.array([0.0, 1.0, 2.0, 3.0])


def main():
    G, founders, reg = build_genotype()
    mat = np.load(f"{H}/hapfreq_matrix.npy")
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    p0 = np.load(f"{H}/hapfreq_p0_seedmix.npy")
    fr = pd.read_csv(f"{H}/hapfreq_registry.csv")
    nhap = mat.shape[1]
    coef = (Tg - 1.5)

    # ---- STAGE 1: per-plot RAW-freq slope + replicate SE + linear binomial sampling floor ----
    pt = lib.pool_table()
    s = pt[(pt.site == SITE) & pt.sampleid.astype(str).isin(sidx)].copy()
    cell, cellN = {}, {}
    for (gen, plot), g in s.groupby(["generation", "plot"]):
        rows = g.sampleid.astype(str).map(sidx).to_numpy()
        w = g.flowerscollected.to_numpy(float); w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        cell[(int(gen), int(plot))] = (mat[rows] * w[:, None]).sum(0) / w.sum()
        cellN[(int(gen), int(plot))] = float(w.sum())
    present = {}
    for (gen, plot) in cell:
        present.setdefault(plot, set()).add(gen)
    plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
    n = len(plots)
    Nmed = np.median([cellN[(g, pl)] for (g, pl) in cellN])

    def bvar(f, N):
        f = np.clip(f, 1e-6, 1 - 1e-6)
        return f * (1 - f) / np.clip(N, 1.0, None)

    slopes = np.zeros((n, nhap)); evo = np.zeros(nhap)
    for j, pl in enumerate(plots):
        cs = [cell[(1, pl)], cell[(2, pl)], cell[(3, pl)]]
        y = np.vstack([p0, cs[0], cs[1], cs[2]])
        slopes[j] = (coef[:, None] * np.clip(y, 0.0, 1.0)).sum(0) / 5.0
        for ti in (1, 2, 3):
            evo += (coef[ti] / 5.0) ** 2 * bvar(cs[ti - 1], 2 * cellN[(ti, pl)])
    evo /= n ** 2
    v0_term = (coef[0] / 5.0) ** 2 * bvar(p0, 2 * Nmed)
    s_all = slopes.mean(0)
    se_all = slopes.std(0, ddof=1) / np.sqrt(n)

    fr["unit"] = fr.chrom + ":" + fr.unit_start.astype(str) + "-" + fr.unit_end.astype(str)
    testable = (fr.covered & fr.panel_freq.between(0.05, 0.95)).to_numpy()
    keep = np.zeros(nhap, bool)
    for u, grp in fr[testable].groupby("unit"):
        kept = grp.index.drop(grp.panel_freq.idxmax()) if len(grp) > 1 else grp.index
        keep[kept.to_numpy()] = True
    ki = np.where(keep)[0]; M = len(ki)
    sb = s_all[ki].copy(); sb -= sb.mean()
    rep = se_all[ki] ** 2
    D0 = v0_term[ki] + np.maximum(rep, evo[ki])
    D0 = np.clip(D0, np.median(D0) * 1e-3, None)

    # ---- STAGE 2: precision-weighted LD-LMM (Woodbury, low-rank) ----
    A = G[:, ki].astype(np.float64)
    pf = A.mean(0); A = (A - pf) / np.sqrt(pf * (1 - pf) + 1e-9)
    nF = A.shape[0]; one = np.ones(M)

    def pieces(tau, se2):
        Winv = 1.0 / (se2 + D0)
        AWA = (A * Winv[None, :]) @ A.T / nF
        cf = cho_factor(np.eye(nF) / tau + AWA, lower=True)

        def Vinv(x):
            Wx = Winv * x
            u = cho_solve(cf, (A @ Wx) / np.sqrt(nF))
            return Wx - Winv * ((A.T @ u) / np.sqrt(nF))
        logdetV = np.sum(np.log(se2 + D0)) + np.linalg.slogdet(np.eye(nF) + tau * AWA)[1]
        return Vinv, logdetV

    def neg_reml(par):
        tau, se2 = np.exp(par)
        Vinv, logdetV = pieces(tau, se2)
        Vi1 = Vinv(one); mu = (one @ Vinv(sb)) / (one @ Vi1)
        r = sb - mu
        return 0.5 * (logdetV + r @ Vinv(r) + np.log(one @ Vi1))

    res = optimize.minimize(neg_reml, x0=np.log([np.var(sb), np.var(sb) / 2]),
                            method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-3})
    tau, se2 = np.exp(res.x)
    Vinv, _ = pieces(tau, se2)
    mu = (one @ Vinv(sb)) / (one @ Vinv(one)); r = sb - mu
    g_b = (tau / nF) * (A.T @ (A @ Vinv(r)))
    resid = r - g_b
    z = resid / np.sqrt(se2 + D0)
    pval = 2 * stats.norm.sf(np.abs(z))
    lam_gc = np.median(stats.chi2.isf(np.clip(pval, 1e-300, 1), 1)) / stats.chi2.ppf(0.5, 1)
    h2 = tau / (tau + se2)

    d = fr.loc[ki, ["chrom", "unit_start", "unit_end", "panel_freq"]].copy()
    d["gid"] = ki                                          # unique global haplotype id (join key)
    d["s"] = sb; d["se"] = np.sqrt(D0); d["linked_bg"] = g_b; d["resid"] = resid
    d["z"] = z; d["p"] = pval
    m = len(d); o = d.p.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate((d.p.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q"] = np.clip(q, 0, 1); bonf = 0.05 / m
    # one-sided p for the TARGET (rose MORE than predicted): resid>0
    d["p_pos"] = stats.norm.sf(z)                          # small when z>>0 (rose-more)
    op = d.p_pos.values.argsort()
    qp = np.empty(m); qp[op] = np.minimum.accumulate((d.p_pos.values[op] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q_pos"] = np.clip(qp, 0, 1)
    d.sort_values("p").to_csv(f"{H}/site{SITE}_block_ld_lmm_linear_sampvar.csv", index=False)
    meta = dict(site=SITE, n_plot=n, M=int(M), nF=int(nF), h2=float(h2),
                sigma_g2=float(tau), sigma_e2=float(se2), lambda_gc=float(lam_gc),
                n_bonf=int((d.p < bonf).sum()), n_fdr05=int((d.q < 0.05).sum()),
                n_pos_fdr05=int((d.q_pos < 0.05).sum()),
                model="linear-scale block LD-LMM + binomial sampling-variance floor")
    json.dump(meta, open(f"{H}/site{SITE}_block_ld_lmm_linear_sampvar_meta.json", "w"), indent=2)

    pos = d[d.resid > 0]; neg = d[d.resid < 0]
    bnd = lambda x: ((x.panel_freq <= 0.10) | (x.panel_freq >= 0.90)).mean()
    print(f"site {SITE} LINEAR + sampling floor: {M:,} blocks, h2={h2:.2f}, lambda_GC={lam_gc:.2f}")
    print(f"  median D0={np.median(D0):.5f}  (replicate-only {np.median(rep):.5f}; binomial floor {np.median(evo[ki]):.5f})")
    print(f"  two-sided hits: Bonferroni {int((d.p<bonf).sum())} | FDR q<0.05 {int((d.q<0.05).sum())}")
    print(f"     of FDR hits: resid>0 (rose-MORE) {int((d.q<0.05).sum() and (d[d.q<0.05].resid>0).sum())}, "
          f"resid<0 (rose-less) {int((d[d.q<0.05].resid<0).sum())}; boundary-freq frac {bnd(d[d.q<0.05]):.2f}")
    print(f"  ONE-SIDED target (rose-MORE, resid>0): FDR q_pos<0.05 = {int((d.q_pos<0.05).sum())}; "
          f"boundary-freq frac {bnd(d[d.q_pos<0.05]):.2f}, median panel_freq {d[d.q_pos<0.05].panel_freq.median():.3f}")
    print("\n  TOP-12 rose-MORE-than-ecotype blocks (one-sided, q_pos):")
    cols = ["chrom", "unit_start", "unit_end", "panel_freq", "s", "linked_bg", "resid", "z", "p_pos", "q_pos"]
    print(d.sort_values("p_pos").head(12)[cols].to_string(index=False))
    print(f"\n[done] {H}/site{SITE}_block_ld_lmm_linear_sampvar.csv + _meta.json")


if __name__ == "__main__":
    main()

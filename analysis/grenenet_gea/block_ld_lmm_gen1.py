#!/usr/bin/env python
"""Founding->GEN1 block LD-LMM (linear scale + binomial sampling floor), per site.

Selection acts strongest in the FIRST generation (founding seedmix -> gen1), and gen1 is available
at far more sites than the full gen1-3 trajectory. So the response here is the one-generation linear
change, available at every site, giving a CONSISTENT per-site statistic for cross-site combination:

  per block b, per plot j:  Delta_{b,j} = f_{b,gen1,j} - p0_b          (linear founding->gen1 change)
  s_b  = mean_j Delta_{b,j}        SE_b = sd_j(Delta_{b,j})/sqrt(n)
  D0_b = p0_b(1-p0_b)/(2 Nmed)  +  max( SE_b^2,  (1/n^2) sum_j f1_{b,j}(1-f1_{b,j})/(2 N1_j) )
         \_ founding sampling (shared)        \_ replicate spread OR gen1 binomial floor

Stage 2 = identical linear LD-LMM (C_LD founder-relationship; BLUP linked background; block-specific
residual r_b = s_b - BLUP; r_b>0 = rose MORE than founder linkage predicts). Writes
site<ID>_block_ld_lmm_gen1.csv (gid,panel_freq,s,se,linked_bg,resid,z,p,q,p_pos,q_pos) + _meta.json.
Env: kmate. SITE via env.
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


def bvar(f, N):
    f = np.clip(f, 1e-6, 1 - 1e-6)
    return f * (1 - f) / np.clip(N, 1.0, None)


def main():
    G, founders, reg = build_genotype()
    mat = np.load(f"{H}/hapfreq_matrix.npy")
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    p0 = np.load(f"{H}/hapfreq_p0_seedmix.npy")
    fr = pd.read_csv(f"{H}/hapfreq_registry.csv")
    nhap = mat.shape[1]

    # ---- STAGE 1: founding->gen1 linear change per plot ----
    pt = lib.pool_table()
    s = pt[(pt.site == SITE) & (pt.generation == 1) & pt.sampleid.astype(str).isin(sidx)].copy()
    cell, cellN = {}, {}
    for plot, g in s.groupby("plot"):
        rows = g.sampleid.astype(str).map(sidx).to_numpy()
        w = g.flowerscollected.to_numpy(float); w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        cell[int(plot)] = (mat[rows] * w[:, None]).sum(0) / w.sum()
        cellN[int(plot)] = float(w.sum())
    plots = sorted(cell)
    n = len(plots)
    Nmed = np.median(list(cellN.values()))

    deltas = np.zeros((n, nhap)); gen1_samp = np.zeros(nhap)
    for j, pl in enumerate(plots):
        f1 = np.clip(cell[pl], 0.0, 1.0)
        deltas[j] = f1 - np.clip(p0, 0.0, 1.0)
        gen1_samp += bvar(f1, 2 * cellN[pl])
    gen1_samp /= n ** 2
    v0_term = bvar(p0, 2 * Nmed)                            # founding sampling var (shared across plots)
    s_all = deltas.mean(0)
    se_all = deltas.std(0, ddof=1) / np.sqrt(n)

    fr["unit"] = fr.chrom + ":" + fr.unit_start.astype(str) + "-" + fr.unit_end.astype(str)
    testable = (fr.covered & fr.panel_freq.between(0.05, 0.95)).to_numpy()
    keep = np.zeros(nhap, bool)
    for u, grp in fr[testable].groupby("unit"):
        kept = grp.index.drop(grp.panel_freq.idxmax()) if len(grp) > 1 else grp.index
        keep[kept.to_numpy()] = True
    ki = np.where(keep)[0]; M = len(ki)
    sb = s_all[ki].copy(); sb -= sb.mean()
    rep = se_all[ki] ** 2
    D0 = v0_term[ki] + np.maximum(rep, gen1_samp[ki])
    D0 = np.clip(D0, np.median(D0) * 1e-3, None)

    # ---- STAGE 2: linear LD-LMM (Woodbury, low-rank) ----
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

    res = optimize.minimize(neg_reml, x0=np.log([np.var(sb) + 1e-12, np.var(sb) / 2 + 1e-12]),
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
    d["gid"] = ki
    d["s"] = sb; d["se"] = np.sqrt(D0); d["linked_bg"] = g_b; d["resid"] = resid; d["z"] = z; d["p"] = pval
    m = len(d); o = d.p.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate((d.p.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q"] = np.clip(q, 0, 1)
    d["p_pos"] = stats.norm.sf(z)
    op = d.p_pos.values.argsort()
    qp = np.empty(m); qp[op] = np.minimum.accumulate((d.p_pos.values[op] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q_pos"] = np.clip(qp, 0, 1)
    d.sort_values("p").to_csv(f"{H}/site{SITE}_block_ld_lmm_gen1.csv", index=False)
    json.dump(dict(site=SITE, n_plot=n, M=int(M), nF=int(nF), h2=float(h2),
                   sigma_g2=float(tau), sigma_e2=float(se2), lambda_gc=float(lam_gc),
                   n_pos_fdr05=int((d.q_pos < 0.05).sum()), Nmed_flowers=float(Nmed)),
              open(f"{H}/site{SITE}_block_ld_lmm_gen1_meta.json", "w"), indent=2)
    print(f"site {SITE} gen1: {n} plots, {M:,} blocks, h2={h2:.2f}, lambda_GC={lam_gc:.2f}, "
          f"rose-more FDR(q_pos<0.05)={int((d.q_pos<0.05).sum())}")


if __name__ == "__main__":
    main()

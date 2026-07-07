#!/usr/bin/env python
"""ONE-SIDED adaptive block-specific selection: blocks that rose MORE than their ecotype
background (residual r_b > 0). Default site 4. NEW, standalone — does NOT touch the other
agent's block_ld_lmm_sampvar / driftnull code, nor block_ld_lmm_temporal / _linear / _binom.

Motivation (audit): the two-sided residual test's significant hits were RISERS that rose LESS
than the founder-BLUP predicted (resid<0) — BLUP-over-prediction artifacts in the WRONG direction.
We want the ADAPTIVE direction: an allele whose carriers' change EXCEEDS what their relatedness
(hitchhiking) predicts -> r_b = s_b - g_b > 0. So we test the UPPER tail only.

Model = the sampling-variance-weighted LD-LMM (logit slope, honest binomial sampling variance with
the FLOWER CENSUS as N_eff, kinship C_LD as the linked background), identical fit; only the test is
one-sided:  p_b = P(Z >= z_b) = Phi(-z_b),  z_b = r_b / sqrt(sigma_e^2 + SE_b^2).
BH-FDR + Bonferroni on the one-sided p. Output site<ID>_block_ld_lmm_oneside.csv + _meta.json +
manhattan (overperformers above 0). Env: kmate. SITE via env.
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
Tg = np.array([0.0, 1.0, 2.0, 3.0])
EPS = 1e-3
SEED_N = 500.0


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

    # ---- STAGE 1: per-(gen,plot) hapfreq + N_eff = flower census; logit slope + sampling var ----
    pt = lib.pool_table()
    s = pt[(pt.site == SITE) & pt.sampleid.astype(str).isin(sidx)].copy()
    cell, neff = {}, {}
    for (gen, plot), g in s.groupby(["generation", "plot"]):
        rows = g.sampleid.astype(str).map(sidx).to_numpy()
        w = g.flowerscollected.to_numpy(float); w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        cell[(int(gen), int(plot))] = (mat[rows] * w[:, None]).sum(0) / w.sum()
        neff[(int(gen), int(plot))] = 2.0 * w.sum()         # census = chromosomes from flowers
    present = {}
    for (gen, plot) in cell:
        present.setdefault(plot, set()).add(gen)
    plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
    n = len(plots)
    coef = (Tg - 1.5)
    slopes = np.zeros((n, nhap)); vsamp = np.zeros(nhap)
    for j, pl in enumerate(plots):
        F = [p0, cell[(1, pl)], cell[(2, pl)], cell[(3, pl)]]
        Ne = [SEED_N, neff[(1, pl)], neff[(2, pl)], neff[(3, pl)]]
        slopes[j] = (coef[:, None] * np.vstack([logit(f) for f in F])).sum(0) / 5.0
        for ti in range(4):
            fc = np.clip(F[ti], EPS, 1 - EPS)
            vsamp += coef[ti] ** 2 / 25.0 / (Ne[ti] * fc * (1 - fc))
    s_all = slopes.mean(0)
    vsamp_b = vsamp / n / n
    vdrift_b = slopes.var(0, ddof=1) / n

    fr["unit"] = fr.chrom + ":" + fr.unit_start.astype(str) + "-" + fr.unit_end.astype(str)
    testable = (fr.covered & fr.panel_freq.between(0.05, 0.95)).to_numpy()
    keep = np.zeros(nhap, bool)
    for u, grp in fr[testable].groupby("unit"):
        kept = grp.index.drop(grp.panel_freq.idxmax()) if len(grp) > 1 else grp.index
        keep[kept.to_numpy()] = True
    ki = np.where(keep)[0]; M = len(ki)
    sb = s_all[ki].copy(); sb -= sb.mean()
    D0 = np.maximum(vdrift_b[ki], vsamp_b[ki]); D0 = np.clip(D0, np.median(D0) * 1e-3, None)

    # ---- STAGE 2: precision-weighted LD-LMM (Woodbury) ----
    A = G[:, ki].astype(np.float64)
    pf = A.mean(0); A = (A - pf) / np.sqrt(pf * (1 - pf) + 1e-9)
    nF = A.shape[0]; one = np.ones(M)

    def pieces(tau, se2):
        Winv = 1.0 / (se2 + D0); AWA = (A * Winv[None, :]) @ A.T / nF
        cf = cho_factor(np.eye(nF) / tau + AWA, lower=True)

        def Vinv(x):
            Wx = Winv * x
            u = cho_solve(cf, (A @ Wx) / np.sqrt(nF))
            return Wx - Winv * ((A.T @ u) / np.sqrt(nF))
        return Vinv, np.sum(np.log(se2 + D0)) + np.linalg.slogdet(np.eye(nF) + tau * AWA)[1]

    def neg_reml(par):
        tau, se2 = np.exp(par); Vinv, ld = pieces(tau, se2)
        Vi1 = Vinv(one); mu = (one @ Vinv(sb)) / (one @ Vi1); r = sb - mu
        return 0.5 * (ld + r @ Vinv(r) + np.log(one @ Vi1))

    tau, se2 = np.exp(optimize.minimize(neg_reml, np.log([np.var(sb), np.var(sb) / 2]),
                                        method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-3}).x)
    Vinv, _ = pieces(tau, se2)
    mu = (one @ Vinv(sb)) / (one @ Vinv(one)); r = sb - mu
    g_b = (tau / nF) * (A.T @ (A @ Vinv(r)))
    resid = r - g_b
    z = resid / np.sqrt(se2 + D0)
    p_one = stats.norm.sf(z)                                 # ONE-SIDED: upper tail (resid > 0 = adaptive)
    h2 = tau / (tau + se2)

    d = fr.loc[ki, ["chrom", "unit_start", "unit_end", "panel_freq"]].copy()
    d["s"] = sb; d["se"] = np.sqrt(D0); d["linked_bg"] = g_b; d["resid"] = resid; d["z"] = z; d["p"] = p_one
    m = len(d); o = d.p.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate((d.p.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q"] = np.clip(q, 0, 1); bonf = 0.05 / m
    # calibration of the one-sided test: under H0 the upper-tail p is U(0,1); lambda via the upper half
    lam_gc = np.median(stats.chi2.isf(np.clip(2 * np.minimum(p_one, 1 - p_one), 1e-300, 1), 1)) / stats.chi2.ppf(0.5, 1)
    d.sort_values("p").to_csv(f"{H}/site{SITE}_block_ld_lmm_oneside.csv", index=False)
    hit = d[d.q < 0.05]
    meta = dict(site=SITE, test="one-sided resid>0 (adaptive: rose MORE than ecotype background)",
                n_plot=n, M=int(M), nF=int(nF), h2=float(h2), sigma_g2=float(tau), sigma_e2=float(se2),
                lambda_gc=float(lam_gc), n_bonf=int((d.p < bonf).sum()), n_fdr05=int((d.q < 0.05).sum()),
                best_q=float(d.q.min()), best_p=float(d.p.min()))
    json.dump(meta, open(f"{H}/site{SITE}_block_ld_lmm_oneside_meta.json", "w"), indent=2)

    print(f"site {SITE}: {M:,} haploblocks — ONE-SIDED adaptive test (resid>0 = rose MORE than ecotype)")
    print(f"  REML h2={h2:.2f}, sigma_e^2={se2:.4g}, two-sided-equiv lambda_GC={lam_gc:.2f}")
    print(f"  adaptive (overperformer) hits: Bonferroni {int((d.p<bonf).sum())} | FDR q<0.05 {int((d.q<0.05).sum())}"
          f"  | best q={d.q.min():.3f} (best p={d.p.min():.2g})")
    print("\n  TOP adaptive overperformers (rose MORE than founder background):")
    print(d.sort_values("p").head(12)[
        ["chrom", "unit_start", "unit_end", "panel_freq", "s", "linked_bg", "resid", "z", "p", "q"]].to_string(index=False))

    # Manhattan — only the adaptive (resid>0) signal; grey = resid<=0 (not tested)
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
    nlp = -np.log10(np.clip(g2.p, 1e-300, 1)); adp = g2.resid.values > 0
    ax.scatter(g2.x[~adp], nlp[~adp], s=4, c="#cccccc", alpha=.4, edgecolors="none", label="resid≤0 (not adaptive)")
    for i, ch in enumerate(CHROMS):
        mk = ((g2.chrom == ch) & adp).to_numpy()
        ax.scatter(g2.x[mk], nlp[mk], s=6, c=["#c0392b", "#e07b39"][i % 2], alpha=.65, edgecolors="none")
    ax.axhline(-np.log10(bonf), color="firebrick", lw=.9, ls="--", label="Bonferroni")
    ax.set_xticks(centers); ax.set_xticklabels(CHROMS)
    ax.set_ylabel(r"$-\log_{10}p$ (one-sided, resid>0)"); ax.set_xlabel("genome position")
    ax.set_title(f"Site {SITE} ADAPTIVE block-specific selection (rose MORE than ecotype) "
                 f"— {int((d.q<0.05).sum())} FDR, best q={d.q.min():.2f}", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False, loc="upper right"); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{H}/site{SITE}_block_ld_lmm_oneside.png", dpi=150)
    print(f"\n[done] {H}/site{SITE}_block_ld_lmm_oneside.csv + _meta.json + manhattan")


if __name__ == "__main__":
    main()

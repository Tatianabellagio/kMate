#!/usr/bin/env python
"""LINEAR-scale block-level Δp + kinship-controlled LD mixed model (default site 4).

NEW, standalone (does NOT touch the logit model block_ld_lmm_temporal.py). This combines the
two views the project converged on:
  * BLOCK-LEVEL Δp (the directly-measured pool allele-frequency change over time), and
  * a GWAS-style KINSHIP control (the founder relationship C_LD),
on the LINEAR (raw-frequency) scale instead of logit.

Why linear: f_b = Σ_{f∈b} h_f, so Δf_b = Σ Δh_f is EXACTLY linear in the founder changes. So the
founder/kinship term C_LD can fully absorb the clade-wide winning (hitchhikers + the "winners are
one clade" confound), and the residual is the allele-specific / convergent signal — with NO
logit-curvature artifact (which made 100% of the logit model's hits frequency-boundary cases).
The residual flags alleles whose carriers won MORE than their relatedness predicts (= the
convergently-shared adaptive allele), resolvable wherever the panel decouples it from clade markers.

STAGE 1 (per block, temporal + replicate, LINEAR scale):
  per plot j: OLS slope beta_{b,j} of RAW freq f over generations t=0(founding),1,2,3;
  s_b = mean_j beta_{b,j} (= block-level Δp rate);  SE_b = sd_j(beta_{b,j})/sqrt(n).
STAGE 2 (precision-weighted LD-LMM):
  s = mu + g + e,  g ~ N(0, sg^2 C_LD),  e_b ~ N(0, se^2 + SE_b^2);  C_LD = A'A/nF (founder-sharing).
  REML (Woodbury, low-rank nF); block-specific residual r_b = s_b - BLUP_b; z_b = r_b/sqrt(se^2+SE_b^2).

Writes site<ID>_block_ld_lmm_linear.csv + _meta.json + manhattan. Env: kmate. SITE via env.
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


def main():
    G, founders, reg = build_genotype()
    mat = np.load(f"{H}/hapfreq_matrix.npy")
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    p0 = np.load(f"{H}/hapfreq_p0_seedmix.npy")
    fr = pd.read_csv(f"{H}/hapfreq_registry.csv")
    nhap = mat.shape[1]

    # ---- STAGE 1: per-plot RAW-frequency trajectory -> per-block linear slope s_b, SE_b ----
    pt = lib.pool_table()
    s = pt[(pt.site == SITE) & pt.sampleid.astype(str).isin(sidx)].copy()
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
    coef = (Tg - 1.5)
    slopes = np.zeros((n, nhap))
    for j, pl in enumerate(plots):
        y = np.vstack([p0, cell[(1, pl)], cell[(2, pl)], cell[(3, pl)]])   # RAW freq (linear)
        slopes[j] = (coef[:, None] * np.clip(y, 0.0, 1.0)).sum(0) / 5.0
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
    D0 = se_all[ki] ** 2
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
    Vir = Vinv(r)
    g_b = (tau / nF) * (A.T @ (A @ Vir))
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
    d["q"] = np.clip(q, 0, 1); bonf = 0.05 / m
    d.sort_values("p").to_csv(f"{H}/site{SITE}_block_ld_lmm_linear.csv", index=False)
    fdrhit = d[d.q < 0.05]
    bdry = float(((fdrhit.panel_freq <= 0.10) | (fdrhit.panel_freq >= 0.90)).mean()) if len(fdrhit) else float("nan")
    meta = dict(site=SITE, scale="linear", n_plot=n, M=int(M), nF=int(nF), h2=float(h2),
                sigma_g2=float(tau), sigma_e2=float(se2), lambda_gc=float(lam_gc),
                n_bonf=int((d.p < bonf).sum()), n_fdr05=int((d.q < 0.05).sum()),
                frac_fdr_boundary=bdry,
                model="block-level Δp (linear/raw freq) + per-plot SE, kinship(C_LD)-controlled LD-LMM")
    json.dump(meta, open(f"{H}/site{SITE}_block_ld_lmm_linear_meta.json", "w"), indent=2)

    print(f"site {SITE}: {n} plots, {M:,} haploblocks (LINEAR scale: raw Δp + kinship-control)")
    print(f"  median SE_b = {np.median(np.sqrt(D0)):.4f}")
    print(f"  REML h2(linkage) = {h2:.2f}, sigma_e^2 = {se2:.4g}")
    print(f"  residual lambda_GC = {lam_gc:.2f}")
    print(f"  block-specific hits: Bonferroni {int((d.p<bonf).sum())} | FDR q<0.05 {int((d.q<0.05).sum())}")
    print(f"  fraction of FDR hits at frequency boundary (pf<=.10 or >=.90): {100*bdry:.0f}% "
          f"(logit model = 100% artifact; background = 28%)")
    print("\n  TOP block-specific blocks (linear):")
    print(d.sort_values("p").head(12)[
        ["chrom", "unit_start", "unit_end", "panel_freq", "s", "se", "linked_bg", "resid", "z", "p", "q"]].to_string(index=False))

    # Manhattan (directional by sign of Δp)
    g2 = d.copy(); g2["mid"] = (g2.unit_start + g2.unit_end) / 2
    g2 = g2.sort_values(["chrom", "unit_start"]).reset_index(drop=True)
    off, centers, x = 0.0, [], np.zeros(len(g2))
    for ch in CHROMS:
        mk = (g2.chrom == ch).to_numpy()
        if not mk.any():
            continue
        x[mk] = g2.mid[mk] + off; centers.append(off + g2.mid[mk].max() / 2); off += g2.mid[mk].max() * 1.02
    g2["x"] = x
    fig, ax = plt.subplots(figsize=(11, 4.4))
    signed = np.sign(g2.s.values) * -np.log10(np.clip(g2.p, 1e-300, 1)); up = g2.s.values > 0
    ax.scatter(g2.x[up], signed[up], s=5, c="#c0392b", alpha=.55, edgecolors="none", label="rose")
    ax.scatter(g2.x[~up], signed[~up], s=5, c="#2471a3", alpha=.55, edgecolors="none", label="fell")
    ax.axhline(0, color="k", lw=.6)
    ax.axhline(-np.log10(bonf), color="firebrick", lw=.9, ls="--", label="Bonferroni"); ax.axhline(np.log10(bonf), color="firebrick", lw=.9, ls="--")
    ax.set_xticks(centers); ax.set_xticklabels(CHROMS)
    ax.set_ylabel(r"signed $-\log_{10}p$ ($\uparrow$ rose, $\downarrow$ fell)"); ax.set_xlabel("genome position")
    ax.set_title(f"Site {SITE} LINEAR block-Δp + kinship LD-LMM (λ_GC={lam_gc:.2f}, h²={h2:.2f}, "
                 f"boundary {100*bdry:.0f}%)", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{H}/site{SITE}_block_ld_lmm_linear.png", dpi=150)
    print(f"\n[done] {H}/site{SITE}_block_ld_lmm_linear.csv + _meta.json + manhattan")


if __name__ == "__main__":
    main()

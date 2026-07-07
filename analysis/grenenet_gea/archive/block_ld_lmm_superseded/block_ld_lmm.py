#!/usr/bin/env python
"""Linkage-aware per-haploblock selection via an LD-covariance mixed model (default site 4).

The locus-level analog of a GWAS LMM. Response = per-block selection coefficient
  s_b = mean over plots of [ logit(f_evolved_b) - logit(f_founding_b) ]   (founding->gen3)
Model:  s = mu*1 + g + e,   g ~ N(0, sg^2 * C_LD),   e ~ N(0, se^2 * I)
where C_LD = the HAPLOTYPE RELATIONSHIP matrix (which haplotypes are carried by the same
founders = the LD among blocks), built from the founder x haplotype indicator G. The
polygenic term g absorbs the LINKED background (s explained by founder-sharing = the
genome-wide ecotype winnowing); linked clusters are down-weighted. The BLUP residual
  r_b = s_b - g_b
is the BLOCK-SPECIFIC selection beyond linkage. REML sets the shrinkage, so the residual
test is calibrated (NOT the among-replicate spread that inflated the earlier tests).

C_LD is low-rank (<=nF founders): we eigendecompose the nF x nF dual matrix. Env: kmate.
SITE via env. Writes block_ld_lmm.csv + manhattan.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats, optimize
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


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def main():
    G, founders, reg = build_genotype()
    mat = np.load(f"{H}/hapfreq_matrix.npy")
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    p0 = np.load(f"{H}/hapfreq_p0_seedmix.npy")                 # founding hapfreq
    nhap = mat.shape[1]

    # per-block selection coefficient s_b = mean_plots [logit(f_evolved)-logit(p0)]
    pt = lib.pool_table()
    s = pt[(pt.site == SITE) & (pt.generation == 3)]
    rows = [sidx[str(x)] for x in s.sampleid if str(x) in sidx]
    dlog = logit(mat[rows]) - logit(p0)[None, :]               # (nplot, nhap)
    s_all = dlog.mean(0)
    n = len(rows)

    # testable haplotypes, k-1 per block
    reg["unit"] = reg.chrom + ":" + reg.start.astype(str) + "-" + reg.end.astype(str)
    fr = pd.read_csv(f"{H}/hapfreq_registry.csv")
    testable = (fr.covered & fr.panel_freq.between(0.05, 0.95)).to_numpy()
    fr["unit"] = fr.chrom + ":" + fr.unit_start.astype(str) + "-" + fr.unit_end.astype(str)
    keep = np.zeros(nhap, bool)
    for u, grp in fr[testable].groupby("unit"):
        kept = grp.index.drop(grp.panel_freq.idxmax()) if len(grp) > 1 else grp.index
        keep[kept.to_numpy()] = True
    ki = np.where(keep)[0]; M = len(ki)
    sb = s_all[ki]; sb = sb - sb.mean()                        # center (profile out mu)

    # haplotype relationship C_LD = A'A/nF (A = founder-standardized indicator columns)
    A = G[:, ki].astype(np.float64)                            # (nF, M)
    p = A.mean(0); A = (A - p) / np.sqrt(p * (1 - p) + 1e-9)   # standardize over founders
    nF = A.shape[0]
    AAt = A @ A.T / nF                                         # (nF, nF) dual
    lam, Wv = np.linalg.eigh(AAt); lam = np.clip(lam, 0, None)
    pos = lam > 1e-8
    lam = lam[pos]; Wv = Wv[:, pos]                            # nonzero modes
    P = (A.T @ Wv) / np.sqrt(lam * nF)                         # (M, r) orthonormal loci-eigvecs
    omega = lam                                               # eigenvalues of C_LD (per mode)
    a = P.T @ sb                                              # (r,) s in eigenbasis
    rss_perp = float(sb @ sb - a @ a)                         # null-space SS (M-r dims)
    r = len(omega)

    # ML over delta = se^2/sg^2
    def negll(logd):
        d = np.exp(logd)
        sg2 = (np.sum(a**2 / (omega + d)) + rss_perp / d) / M
        return 0.5 * (M * np.log(sg2 + 1e-300) + np.sum(np.log(omega + d)) + (M - r) * np.log(d) + M)
    res = optimize.minimize_scalar(negll, bounds=(-12, 12), method="bounded")
    delta = float(np.exp(res.x))
    sg2 = (np.sum(a**2 / (omega + delta)) + rss_perp / delta) / M
    se2 = sg2 * delta
    h2 = sg2 / (sg2 + se2)

    # BLUP linked background + block-specific residual
    shrink = omega / (omega + delta)
    g_b = P @ (shrink * a)                                     # linked background per block
    resid = sb - g_b
    z = resid / np.sqrt(se2)
    p_r = 2 * stats.norm.sf(np.abs(z))
    lam_gc = np.median(stats.chi2.isf(np.clip(p_r, 1e-300, 1), 1)) / stats.chi2.ppf(0.5, 1)

    d = fr.loc[ki, ["chrom", "unit_start", "unit_end", "panel_freq"]].copy()
    d["s"] = sb; d["linked_bg"] = g_b; d["resid"] = resid; d["z"] = z; d["p"] = p_r
    m = len(d); o = d.p.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate((d.p.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q"] = np.clip(q, 0, 1)
    bonf = 0.05 / m
    d.sort_values("p").to_csv(f"{H}/site{SITE}_block_ld_lmm.csv", index=False)

    print(f"site {SITE}: {M:,} haploblocks, {nF} founders, C_LD rank {r}")
    print(f"  REML: h2(linked/founder-explained) = {h2:.2f}, delta = {delta:.3f}")
    print(f"  s var explained by LD background = {1 - np.var(resid)/np.var(sb):.2f}")
    print(f"  residual lambda_GC = {lam_gc:.2f}  (target ~1)")
    print(f"  block-specific hits: Bonferroni(p<{bonf:.1e}) {int((d.p<bonf).sum())}; "
          f"FDR q<0.05 {int((d.q<0.05).sum())}; q<0.10 {int((d.q<0.10).sum())}")
    print("\n  TOP block-specific blocks (s beyond the linked background):")
    print(d.sort_values("p").head(12)[
        ["chrom", "unit_start", "unit_end", "panel_freq", "s", "linked_bg", "resid", "z", "p", "q"]].to_string(index=False))

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
    ax.set_title(f"Site {SITE} block-specific selection (LD-LMM residual, linkage down-weighted) "
                 f"λ_GC={lam_gc:.2f}, h²_LD={h2:.2f}", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{H}/site{SITE}_block_ld_lmm.png", dpi=150)
    print(f"\n[done] {H}/site{SITE}_block_ld_lmm.csv + manhattan")


if __name__ == "__main__":
    main()

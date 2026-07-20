#!/usr/bin/env python
"""Analytic Omega-conditioned per-haploblock contrast (our fast model; Omega FROM BayPass).

BayPass's slow MCMC is only needed to estimate Omega (the npop x npop coancestry covariance)
+ per-locus ancestral freqs. Given Omega, the per-locus XtX / C2 contrast statistics are
closed-form GLS/Mahalanobis quantities — computed here for ALL loci at once in numpy (no
per-locus MCMC). This keeps the temporal contrast (founding vs evolved) Omega-conditioned and
per-block, but runs in seconds.

Model (Gunther & Coop 2013): standardized allele freq z_i = (p_i - pi_i)/sqrt(pi_i(1-pi_i)),
Cov(z_i) = Omega under neutrality. Per locus:
  pi_i  = GLS ancestral freq = (1' W p_i)/(1' W 1),  W = Omega^-1
  XtX_i = z_i' W z_i                      ~ chi2_(npop)   (overall differentiation)
  C2_i  = (c' W z_i)^2 / (c' W c)         ~ chi2_1        (contrast c = founding -1 / evolved +1)

Inputs: OUTDIR/{geno.txt, contrast.txt, loci.csv, core_mat_omega.out}. Writes omega_contrast.csv
+ manhattan. Env: kmate. OUTDIR via env.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = "analysis/grenenet_gea/archive/window_hapfreq_retired/hapfreq"
SITE = int(os.environ.get("SITE", 4))
OUTDIR = os.environ.get("OUTDIR", f"{H}/baypass_temporal_site{SITE}")
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
EPS = 1e-6


def main():
    G = np.loadtxt(f"{OUTDIR}/geno.txt", dtype=float)          # (M, 2*npop)
    c = np.loadtxt(f"{OUTDIR}/contrast.txt", dtype=float)      # (npop,)
    omega = np.loadtxt(os.environ.get("OMEGA", f"{OUTDIR}/core_mat_omega.out"))  # (npop, npop)
    loci = pd.read_csv(f"{OUTDIR}/loci.csv")
    npop = len(c); M = G.shape[0]
    n1 = G[:, 0::2]; n2 = G[:, 1::2]                           # (M, npop)
    F = (n1 / np.clip(n1 + n2, 1, None)).T                    # (npop, M) allele freq per pop

    # ridge for robustness (MoM Omega is rank npop-1 from mean-centering; BayPass Omega is
    # full-rank so ridge≈0 is a no-op there). RIDGE = fraction of mean diagonal.
    rf = float(os.environ.get("RIDGE", 0.0))
    if rf > 0:
        omega = omega + np.eye(npop) * rf * np.diag(omega).mean()
    W = np.linalg.inv(omega)
    wsum = W.sum(1)                                            # = 1' W (symmetric)
    denom1 = W.sum()
    pi = (wsum @ F) / denom1                                  # (M,) GLS ancestral freq
    pi = np.clip(pi, EPS, 1 - EPS)
    z = (F - pi[None, :]) / np.sqrt(pi * (1 - pi))[None, :]    # (npop, M) standardized

    Wz = W @ z                                                # (npop, M)
    xtx = (z * Wz).sum(0)                                     # ~ chi2_npop
    cw = W @ c                                                # (npop,)
    c2 = (cw @ z) ** 2 / (c @ cw)                             # (M,) ~ chi2_1

    p_c2 = stats.chi2.sf(c2, 1)
    p_xtx = stats.chi2.sf(xtx, npop)
    lam = np.median(c2) / stats.chi2.ppf(0.5, 1)

    d = loci.copy()
    d["xtx"] = xtx; d["c2"] = c2; d["p_c2"] = p_c2; d["p_xtx"] = p_xtx
    d["contrast_sign"] = np.sign(cw @ z)                       # + = up in evolved
    m = len(d); o = d.p_c2.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate((d.p_c2.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q_c2"] = np.clip(q, 0, 1)
    bonf = 0.05 / m
    d.sort_values("p_c2").to_csv(f"{OUTDIR}/omega_contrast.csv", index=False)

    print(f"site {SITE}: {M:,} haploblock loci, {npop} pops")
    print(f"  C2 lambda_GC = {lam:.2f}  (target ~1 if Omega calibrates)")
    print(f"  C2 hits: Bonferroni(p<{bonf:.1e}) {int((d.p_c2<bonf).sum())}; FDR q<0.05 {int((d.q_c2<0.05).sum())}; "
          f"q<0.10 {int((d.q_c2<0.10).sum())}")
    print("\n  TOP blocks (founding -> evolved, Omega-conditioned):")
    print(d.sort_values("p_c2").head(12)[
        ["chrom", "unit_start", "unit_end", "c2", "contrast_sign", "p_c2", "q_c2"]].to_string(index=False))

    # Manhattan
    g = d.copy(); g["mid"] = (g.unit_start + g.unit_end) / 2
    g = g.sort_values(["chrom", "unit_start"]).reset_index(drop=True)
    off, centers, x = 0.0, [], np.zeros(len(g))
    for ch in CHROMS:
        mk = (g.chrom == ch).to_numpy()
        if not mk.any():
            continue
        x[mk] = g.mid[mk] + off; centers.append(off + g.mid[mk].max() / 2); off += g.mid[mk].max() * 1.02
    g["x"] = x
    fig, ax = plt.subplots(figsize=(11, 4.2))
    nlp = -np.log10(np.clip(g.p_c2, 1e-300, 1))
    for i, ch in enumerate(CHROMS):
        mk = (g.chrom == ch).to_numpy()
        ax.scatter(g.x[mk], nlp[mk], s=5, c=["#3b4cc0", "#7aa0c4"][i % 2], alpha=.5, edgecolors="none", rasterized=True)
    ax.axhline(-np.log10(bonf), color="firebrick", lw=.9, ls="--", label="Bonferroni")
    ax.set_xticks(centers); ax.set_xticklabels(CHROMS)
    ax.set_ylabel(r"$-\log_{10}p$ (C2)"); ax.set_xlabel("genome position")
    ax.set_title(f"Site {SITE} Omega-conditioned temporal C2 (founding→evolved haploblocks) "
                 f"λ_GC={lam:.2f}", fontsize=11, loc="left")
    ax.legend(fontsize=8, frameon=False); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{OUTDIR}/omega_contrast_manhattan.png", dpi=150)
    print(f"\n[done] {OUTDIR}/omega_contrast.csv + manhattan")


if __name__ == "__main__":
    main()

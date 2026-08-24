#!/usr/bin/env python
"""Regime-correct per-haploblock TEMPORAL selection test (default site 4).

Our regime = E&R from a common founding: selection makes replicate plots CONVERGE (parallel),
and is genome-wide (whole winning ecotypes hitchhike). So (unlike BayPass/Omega, which is built
for DIVERGED natural populations) the signal is CONSISTENCY across replicates beyond drift, and
the structure to remove is the genome-wide ecotype winnowing.

Per haplotype (block u, cluster c) we form the HITCHHIKING-CONDITIONED residual:
    r_{j} = f_obs_{j}  -  f_exp_{j}
    f_obs = local block haplotype freq (hapfreq_matrix = sum_{f in c} h_local(u,f,j))
    f_exp = SAME membership projected through the GENOME-WIDE founder freqs
            = global_h(j) @ G   (what the block would be if it just rode its ecotypes)
We subtract the founding residual (seedmix) so r is the EVOLVED block-specific deviation
(≈0 where the block tracks its ecotype; ≠0 only where recombination decoupled it).

Block-specific selection = r consistent across the 11 plots beyond drift:
    z = mean_j(r) / (sd_j(r)/sqrt(nplot)),  p ~ t(nplot-1);  + sign-concordance.
Conditioning out the genome-wide winnowing should make this SPARSE and calibrated.

Env: kmate. SITE via env. Writes block_residual_consistency.csv + manhattan.
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import build_genotype
from ecotype_selection_site import genome_h

H = "results/grenenet_gea/hapfreq"
WIN = "results/grenenet_kmate_window"
SEED = "results/grenenet_kmate_window_seedmix"
SITE = int(os.environ.get("SITE", 4))
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]


def main():
    G, founders, reg = build_genotype()                        # (231, nhap) founder x haplotype
    mat = np.load(f"{H}/hapfreq_matrix.npy")                    # (2168, nhap) f_obs (local)
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    nhap = mat.shape[1]
    fr = pd.read_csv(f"{H}/hapfreq_registry.csv")

    # founding residual: seedmix local hapfreq - seedmix global projected
    p0_local = np.load(f"{H}/hapfreq_p0_seedmix.npy")           # (nhap,) local seedmix hapfreq
    seeds = sorted({p.split("/")[-1].split("_Chr")[0] for p in glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz")})
    seed_gh = np.mean([genome_h(s, SEED) for s in seeds], 0)    # (231,) seedmix global_h
    p0_exp = seed_gh @ G                                        # (nhap,) global-projected seedmix
    r0 = p0_local - p0_exp                                      # founding residual (~0)

    # site gen-3 plots: per-plot local hapfreq (f_obs) and global projection (f_exp)
    pt = lib.pool_table()
    s = pt[(pt.site == SITE) & (pt.generation == 3)].copy()
    s = s[s.sampleid.astype(str).isin(sidx)]
    R = []                                                      # per-plot evolved residual change
    for plot, g in s.groupby("plot"):
        samp = str(g.sampleid.iloc[0])
        gh = genome_h(samp, WIN)
        if gh is None:
            continue
        f_obs = mat[sidx[samp]]
        f_exp = gh @ G
        R.append((f_obs - f_exp) - r0)                         # subtract founding residual
    R = np.vstack(R)                                           # (nplot, nhap)
    n = R.shape[0]

    rbar = R.mean(0); sd = R.std(0, ddof=1)
    z = rbar / np.clip(sd / np.sqrt(n), 1e-9, None)
    p = 2 * stats.t.sf(np.abs(z), df=n - 1)
    nup = (R > 0).sum(0); conc = np.maximum(nup, n - nup) / n

    # testable haplotypes: covered, panel_freq 0.05-0.95, k-1 per block (drop ref)
    fr["unit"] = fr.chrom + ":" + fr.unit_start.astype(str) + "-" + fr.unit_end.astype(str)
    testable = (fr.covered & fr.panel_freq.between(0.05, 0.95)).to_numpy()
    keep = np.zeros(nhap, bool)
    for u, grp in fr[testable].groupby("unit"):
        kept = grp.index.drop(grp.panel_freq.idxmax()) if len(grp) > 1 else grp.index
        keep[kept.to_numpy()] = True
    ki = np.where(keep)[0]

    d = fr.loc[ki, ["chrom", "unit_start", "unit_end", "panel_freq"]].copy()
    d["resid_mean"] = rbar[ki]; d["z"] = z[ki]; d["p"] = p[ki]; d["concordance"] = conc[ki]
    m = len(d); o = d.p.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate((d.p.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q"] = np.clip(q, 0, 1)
    bonf = 0.05 / m
    lam = np.median(stats.chi2.isf(np.clip(d.p, 1e-300, 1), 1)) / stats.chi2.ppf(0.5, 1)
    d.sort_values("p").to_csv(f"{H}/site{SITE}_block_residual_consistency.csv", index=False)

    print(f"site {SITE}: {n} plots, {m:,} testable haploblocks")
    print(f"  residual |r| median {np.median(np.abs(d.resid_mean)):.4f} (≈0 = tracks ecotype)")
    print(f"  lambda_GC = {lam:.2f}  (hitchhiking conditioned out → expect ~1 + sparse tail)")
    print(f"  hits: Bonferroni(p<{bonf:.1e}) {int((d.p<bonf).sum())}; FDR q<0.05 {int((d.q<0.05).sum())}; "
          f"q<0.10 {int((d.q<0.10).sum())}")
    print("\n  TOP block-specific blocks (decoupled + consistent):")
    print(d.sort_values("p").head(12)[
        ["chrom", "unit_start", "unit_end", "panel_freq", "resid_mean", "z", "concordance", "p", "q"]].to_string(index=False))

    # Manhattan
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
    ax.set_title(f"Site {SITE} block-specific temporal selection (hitchhiking-conditioned residual, "
                 f"consistency vs drift) λ_GC={lam:.2f}", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{H}/site{SITE}_block_residual_consistency.png", dpi=150)
    print(f"\n[done] {H}/site{SITE}_block_residual_consistency.csv + manhattan")


if __name__ == "__main__":
    main()

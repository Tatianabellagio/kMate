#!/usr/bin/env python
"""Why do the two site-4 analyses (point) differ? Show that analysis-2 == analysis-1 + kinship.

Analysis 1 = per-haploblock selection: slope of the pooled hap-frequency over time (no
             structure control; frequency-weighted; smeared by hitchhiking).
Analysis 2 = founder GWAS: effect of carrying a haplotype on per-founder selection s_f,
             with vs without kinship.

For every kept marker we compute, on the SAME haplotypes:
  s1     = analysis-1 logit slope of pooled hap freq  (Hbar @ G over time)
  b_naive= analysis-2 effect, NO kinship
  b_kin  = analysis-2 effect, WITH kinship
Prediction: corr(s1, b_naive) HIGH (both uncontrolled -> agree), corr(s1, b_kin) LOW
(kinship removes the shared structure). That is the answer: they agree until kinship.

Env: kmate. SITE via env.
"""
import os, sys, glob
import numpy as np
import pandas as pd
from scipy import stats, optimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import build_genotype, emma_reml_delta
from ecotype_selection_site import genome_h

WIN = "results/grenenet_kmate_window"
SEED = "results/grenenet_kmate_window_seedmix"
OUT = "results/grenenet_gea/hapfreq"
SITE = int(os.environ.get("SITE", 4))
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
EPS, MAC = 1e-3, 3
Tg = np.array([0.0, 1.0, 2.0, 3.0])


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def main():
    G, founders, reg = build_genotype()                  # (231, nhap), founder x haplotype
    nF = len(founders)

    # ---- founder pooled trajectory at the site (analysis-1 ingredient) ----
    seeds = sorted({p.split("/")[-1].split("_Chr")[0]
                    for p in glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz")})
    p0 = np.vstack([genome_h(s, SEED) for s in seeds]).mean(0)
    pt = lib.pool_table()
    s = pt[pt.site == SITE].copy()
    s = s[s.sampleid.astype(str).apply(
        lambda x: os.path.exists(f"{WIN}/{x}_Chr1.h_blocks_per_chrom.npz"))]
    cell = {}
    for (gen, plot), g in s.groupby(["generation", "plot"]):
        hs = [genome_h(str(r.sampleid), WIN) for _, r in g.iterrows()]
        hs = [h for h in hs if h is not None]
        if hs:
            cell[(int(gen), int(plot))] = np.mean(hs, 0)
    present = {}
    for (gen, plot) in cell:
        present.setdefault(plot, set()).add(gen)
    plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
    Hbar = np.zeros((4, nF)); Hbar[0] = p0
    for gen in (1, 2, 3):
        Hbar[gen] = np.mean([cell[(gen, pl)] for pl in plots], 0)   # pooled over plots

    # phenotype (per-founder selection s_bar) for present founders
    ph = pd.read_csv(f"{OUT}/site{SITE}_ecotype_selection.csv")
    ph["founder"] = ph.founder.astype(str)
    pmap = dict(zip(ph.founder, ph.s_bar))
    fi = np.array([i for i, f in enumerate(founders) if f in pmap])
    y = np.array([pmap[founders[i]] for i in fi])
    n = len(y); yq = stats.norm.ppf((stats.rankdata(y) - 0.5) / n)
    G194 = G[fi]

    # kept markers: drop per-block reference + MAC (same as the GWAS)
    cnt = G194.sum(0); blk = reg.block.to_numpy()
    order = np.lexsort((-cnt, blk)); sb = blk[order]
    is_ref = np.zeros(len(cnt), bool)
    is_ref[order[np.concatenate([[True], sb[1:] != sb[:-1]])]] = True
    poly = (~is_ref) & (cnt >= MAC) & (cnt <= n - MAC)
    Gp = G194[:, poly]; M = Gp.shape[1]

    # ---- analysis 1: pooled hap-freq slope per marker (uses ALL founders for the pool) ----
    PQ = Hbar @ G[:, poly]                                # (4, M) hap pool freq over time
    s1 = ((Tg - 1.5)[:, None] * logit(PQ)).sum(0) / 5.0   # logit slope

    # ---- analysis 2: founder-GWAS effect, naive vs kinship ----
    yc = yq - yq.mean(); Gc = Gp - Gp.mean(0)
    b_naive = (yc @ Gc) / np.clip((Gc**2).sum(0), 1e-12, None)     # OLS slope (no kinship)

    p = Gp.mean(0); Z = (Gp - p) / np.sqrt(p * (1 - p) + 1e-9); K = Z @ Z.T / M
    K = (K + K.T) / 2
    lam, U = np.linalg.eigh(K); lam = np.clip(lam, 1e-9, None)
    yr = U.T @ yq; Xr0 = U.T @ np.ones((n, 1))
    delta = emma_reml_delta(yr, Xr0, lam); w = 1.0 / (lam + delta)
    Br = U.T @ Gp; a = Xr0[:, 0]
    Saa = (w * a * a).sum(); Say = (w * a * yr).sum()
    Sab = ((w * a)[None, :] @ Br).ravel(); Sbb = (w[:, None] * Br**2).sum(0)
    Sby = ((w * yr)[None, :] @ Br).ravel()
    detm = Saa * Sbb - Sab**2
    b_kin = (Saa * Sby - Sab * Say) / np.clip(detm, 1e-30, None)

    r_naive = np.corrcoef(s1, b_naive)[0, 1]
    r_kin = np.corrcoef(s1, b_kin)[0, 1]
    print(f"site {SITE}: {M:,} markers, {n} founders")
    print(f"  corr( analysis-1 block selection , analysis-2 effect ):")
    print(f"     NO kinship : r = {r_naive:.3f}   <- they AGREE (both structure-confounded)")
    print(f"     kinship    : r = {r_kin:.3f}   <- agreement BREAKS (kinship removed the shared structure)")
    print(f"  variance of analysis-1 explained by structure = {r_naive**2 - r_kin**2:.2f} "
          f"(r²_naive {r_naive**2:.2f} - r²_kin {r_kin**2:.2f})")

    fig, ax = plt.subplots(1, 2, figsize=(12, 5.5), sharey=False)
    for a_, b, tit, rr in [(ax[0], b_naive, "analysis-2 NO kinship", r_naive),
                           (ax[1], b_kin, "analysis-2 WITH kinship", r_kin)]:
        a_.scatter(s1, b, s=4, alpha=.25, edgecolors="none", c="#34495e", rasterized=True)
        a_.set_xlabel("analysis 1: block selection $s_1$ (pooled hap-freq slope)")
        a_.set_ylabel(f"{tit}: founder-GWAS effect")
        a_.set_title(f"{tit}\nr = {rr:.3f}", loc="left", fontsize=11)
        a_.axhline(0, color="k", lw=.5); a_.axvline(0, color="k", lw=.5)
        a_.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Site {SITE}: the two analyses agree UNTIL you add kinship "
                 f"(which is what 'controlling structure' means)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = f"{OUT}/site{SITE}_compare_analyses.png"
    fig.savefig(out, dpi=150)
    print(f"\n[done] {out}")


if __name__ == "__main__":
    main()

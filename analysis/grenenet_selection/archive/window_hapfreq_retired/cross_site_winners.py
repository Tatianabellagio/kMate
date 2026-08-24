#!/usr/bin/env python
"""Do different CLADES win at different sites, and does it track CLIMATE? (cross-site winners).

Decides whether cross-site / GxE can even break the clade-vs-allele collinearity: if the SAME
clade wins everywhere (generalist), cross-site cannot help; if DIFFERENT clades win in different
climates (local adaptation), there is clade×climate structure to exploit.

Per site (with persistent gen 1/2/3 plots): per-founder selection s_f = mean-over-plots logit-slope
of its genome-wide frequency over gens 0-3. Then:
  (1) cross-site correlation of the per-founder selection vectors -> generalist (high r) vs local (low/structured),
  (2) founder kinship PCA; per-site 'winning-clade' = selection-weighted founder-PC centroid;
      correlate the winning-clade PC with site climate (bio1) -> does climate drive which clade wins.
Writes cross_site_winners.{csv,png}. Env: kmate.
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

H = "analysis/grenenet_gea/archive/window_hapfreq_retired/hapfreq"
WIN = lib.OUT; SEED = lib.SEEDMIX
Tg = np.array([0.0, 1.0, 2.0, 3.0]); EPS = 1e-3


def logit(p):
    p = np.clip(p, EPS, 1 - EPS); return np.log(p / (1 - p))


def main():
    G, founders, reg = build_genotype(); nF = len(founders)
    seeds = sorted({p.split("/")[-1].split("_Chr")[0] for p in glob.glob(f"{SEED}/*_Chr1.h_per_chrom.npz")})
    p0 = np.mean([genome_h(s, SEED) for s in seeds], 0)
    clim = lib.load_climate()                                  # site-indexed bio1..19
    pt = lib.pool_table()
    pt = pt[pt.sampleid.astype(str).apply(lambda x: os.path.exists(f"{WIN}/{x}_Chr1.h_per_chrom.npz"))]

    coef = (Tg - 1.5); gh_cache = {}
    def gh(samp):
        if samp not in gh_cache:
            gh_cache[samp] = genome_h(samp, WIN)
        return gh_cache[samp]

    S = {}                                                     # site -> per-founder selection (nF,)
    for site, sd in pt.groupby("site"):
        cell = {}
        for (gen, plot), g in sd.groupby(["generation", "plot"]):
            hs = [gh(str(x)) for x in g.sampleid]; hs = [h for h in hs if h is not None]
            if hs:
                cell[(int(gen), int(plot))] = np.mean(hs, 0)
        present = {}
        for (gen, plot) in cell:
            present.setdefault(plot, set()).add(gen)
        plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
        if len(plots) < 2:
            continue
        sl = []
        for pl in plots:
            y = np.vstack([logit(p0), logit(cell[(1, pl)]), logit(cell[(2, pl)]), logit(cell[(3, pl)])])
            sl.append((coef[:, None] * y).sum(0) / 5.0)
        S[int(site)] = np.mean(sl, 0)
    sites = sorted(S)
    Smat = np.vstack([S[s] for s in sites])                   # (nsite, nF)
    bio1 = clim["bio1"].reindex(sites).to_numpy(float)
    print(f"{len(sites)} sites with persistent plots | bio1 range {np.nanmin(bio1):.1f}–{np.nanmax(bio1):.1f}°C")

    # (1) cross-site agreement on winners
    R = np.corrcoef(Smat)
    iu = np.triu_indices(len(sites), 1)
    print(f"\n(1) cross-site correlation of per-founder selection:")
    print(f"    mean pairwise r = {R[iu].mean():.2f}  (high ~ generalist/same winners; low ~ local/different winners)")
    # does cross-site similarity decay with climate distance?
    cd = np.abs(bio1[:, None] - bio1[None, :])[iu]
    ok = np.isfinite(cd)
    rcd = np.corrcoef(R[iu][ok], cd[ok])[0, 1]
    print(f"    corr(site-similarity, climate distance) = {rcd:.2f}  (negative ~ climate structures the winners)")

    # (2) founder kinship PCA + winning-clade centroid per site vs climate
    p = G.mean(0); Z = (G - p) / np.sqrt(p * (1 - p) + 1e-9); K = (Z @ Z.T) / G.shape[1]
    lam, V = np.linalg.eigh(K); pc = V[:, ::-1][:, :2] * np.sqrt(np.clip(lam[::-1][:2], 0, None))
    # winning-clade PC = selection-weighted founder centroid (top-decile winners per site)
    cen = np.zeros((len(sites), 2))
    for i, s in enumerate(sites):
        w = np.clip(Smat[i], 0, None); thr = np.quantile(Smat[i], 0.9); m = Smat[i] >= thr
        cen[i] = (pc[m] * w[m, None]).sum(0) / w[m].sum() if w[m].sum() > 0 else pc[m].mean(0)
    for k in (0, 1):
        r = np.corrcoef(cen[np.isfinite(bio1), k], bio1[np.isfinite(bio1)])[0, 1]
        print(f"\n(2) winning-clade PC{k+1} vs bio1: r = {r:.2f}  (|r| large ~ climate drives which clade wins)")

    out = pd.DataFrame({"site": sites, "bio1": bio1, "win_pc1": cen[:, 0], "win_pc2": cen[:, 1],
                        "mean_r_to_others": [np.delete(R[i], i).mean() for i in range(len(sites))]})
    out.to_csv(f"{H}/cross_site_winners.csv", index=False)

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))
    ax[0].scatter(pc[:, 0], pc[:, 1], s=15, c="#cccccc", edgecolors="none", zorder=1)
    sc = ax[0].scatter(cen[:, 0], cen[:, 1], c=bio1, cmap="coolwarm", s=130, edgecolors="k", lw=.5, zorder=3)
    plt.colorbar(sc, ax=ax[0], label="site bio1 (°C)")
    ax[0].set_xlabel("founder PC1 (kinship)"); ax[0].set_ylabel("PC2")
    ax[0].set_title("winning-clade centroid per site, colored by climate\n(spread+climate-ordered ⇒ local; one spot ⇒ generalist)", fontsize=9, loc="left")
    ax[0].spines[["top", "right"]].set_visible(False)
    rr = np.corrcoef(cen[np.isfinite(bio1), 0], bio1[np.isfinite(bio1)])[0, 1]
    ax[1].scatter(bio1, cen[:, 0], s=60, c=bio1, cmap="coolwarm", edgecolors="k", lw=.4)
    ax[1].set_xlabel("site bio1 (°C)"); ax[1].set_ylabel("winning-clade PC1")
    ax[1].set_title(f"does climate drive which clade wins?  r={rr:.2f}", fontsize=10, loc="left")
    ax[1].spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Cross-site winners ({len(sites)} sites): clade×climate structure?  "
                 f"mean cross-site r={R[iu].mean():.2f}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(f"{H}/cross_site_winners.png", dpi=150)
    print(f"\n[done] {H}/cross_site_winners.csv + png")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Clade-level local adaptation on the SAME 30 sites + SAME trait as the multisite founder GWAS.

The original cross_site_winners.py used a stricter inclusion (>=2 plots surviving all of gen 1+2+3)
and a logit fixed-4-gen slope -> only 17 (mostly cold, long-surviving) sites. That is NOT the site
set the multisite GWAS / permutation null run on (30 sites: gen1 anchor + full available trajectory,
linear OLS slope + QN, exclude site 33). This recomputes the winning-clade-PC1-vs-bio1 signal on the
EXACT 30-site set and trait, by importing build_trait_site/EXCLUDE_SITES from founder_gwas_multisite,
so the clade panel in the notebook is apples-to-apples with the Manhattan + FWER null.

Writes cross_site_winners_30.{csv,png}. Env: kmate. N_PERM (10000), PERM_SEED (0).
"""
import os, sys, glob, json
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
from founder_gwas_multisite import build_trait_site, EXCLUDE_SITES, TRAIT_GENS

H = "results/grenenet_gea/hapfreq"
WIN = "results/grenenet_kmate_window"; SEED = "results/grenenet_kmate_window_seedmix"
N_PERM = int(os.environ.get("N_PERM", 10000))
PERM_SEED = int(os.environ.get("PERM_SEED", 0))
SUFFIX = os.environ.get("OUT_SUFFIX", "")   # "_clq90" pairs with the clq0.9 GWAS/membership run


def main():
    G, founders, reg = build_genotype(); nF = len(founders)
    seeds = sorted({p.split("/")[-1].split("_Chr")[0] for p in glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz")})
    p0 = np.mean([genome_h(s, SEED) for s in seeds], 0)
    clim = lib.load_climate()
    pt = lib.pool_table()
    pt = pt[pt.sampleid.astype(str).apply(lambda x: os.path.exists(f"{WIN}/{x}_Chr1.h_blocks_per_chrom.npz"))]
    gh_cache = {}
    def gh(s):
        if s not in gh_cache: gh_cache[s] = genome_h(s, WIN)
        return gh_cache[s]

    # ---- replicate the multisite GWAS 30-site selection + trait EXACTLY ----
    sites, Srows, present_by_site = [], [], []
    for site, sd in pt.groupby("site"):
        if int(site) in EXCLUDE_SITES:
            continue
        b1 = float(clim["bio1"].reindex([int(site)]).iloc[0]) if int(site) in clim.index else np.nan
        if not np.isfinite(b1):
            continue
        s_f, present = build_trait_site(sd, p0, gh, nF)          # linear OLS slope, full available trajectory
        if s_f is None:
            continue
        sites.append(int(site)); Srows.append(s_f); present_by_site.append(present)
    sites = np.array(sites); Smat = np.vstack(Srows)             # (S, nF) per-founder selection
    bio1 = clim["bio1"].reindex(sites).to_numpy(float)
    S = len(sites)
    print(f"{S} sites (same set as multisite GWAS) | trait gens {list(TRAIT_GENS)} | "
          f"bio1 {np.nanmin(bio1):.1f}-{np.nanmax(bio1):.1f}C | excluded {sorted(EXCLUDE_SITES)}")

    # ---- founder kinship PCA + winning-clade centroid per site (identical to cross_site_winners.py) ----
    p = G.mean(0); Zg = (G - p) / np.sqrt(p * (1 - p) + 1e-9); K = (Zg @ Zg.T) / G.shape[1]
    lam, V = np.linalg.eigh(K); pc = V[:, ::-1][:, :2] * np.sqrt(np.clip(lam[::-1][:2], 0, None))
    cen = np.zeros((S, 2))
    for i in range(S):
        w = np.clip(Smat[i], 0, None); thr = np.quantile(Smat[i], 0.9); mk = Smat[i] >= thr
        cen[i] = (pc[mk] * w[mk, None]).sum(0) / w[mk].sum() if w[mk].sum() > 0 else pc[mk].mean(0)

    R = np.corrcoef(Smat); iu = np.triu_indices(S, 1)
    fin = np.isfinite(bio1)
    out = pd.DataFrame({"site": sites, "bio1": bio1, "win_pc1": cen[:, 0], "win_pc2": cen[:, 1],
                        "n_gen": [len(g) for g in present_by_site],
                        "mean_r_to_others": [np.delete(R[i], i).mean() for i in range(S)]})
    out.to_csv(f"{H}/cross_site_winners_30{SUFFIX}.csv", index=False)

    # ---- winning-clade PC1 vs bio1 + permutation p (both PCs; sign-free) ----
    res = {}
    for k in (0, 1):
        x, y = bio1[fin], cen[fin, k]
        r = float(np.corrcoef(x, y)[0, 1])
        rng = np.random.default_rng(PERM_SEED)
        nd = np.array([abs(np.corrcoef(rng.permutation(x), y)[0, 1]) for _ in range(N_PERM)])
        pp = float((nd >= abs(r)).mean())
        res[f"pc{k+1}"] = dict(r=r, perm_p=pp)
        print(f"(2) winning-clade PC{k+1} vs bio1: r = {r:.2f}, perm p = {pp:.4f}  ({fin.sum()} sites)")
    print(f"(1) mean cross-site r of per-founder selection = {R[iu].mean():.2f}")
    json.dump(dict(n_sites=int(S), n_perm=N_PERM, mean_crosssite_r=float(R[iu].mean()),
                   trait_gens=list(TRAIT_GENS), excluded=sorted(EXCLUDE_SITES), **res),
              open(f"{H}/cross_site_winners_30{SUFFIX}.json", "w"), indent=2)

    # ---- figure ----
    r1, p1 = res["pc1"]["r"], res["pc1"]["perm_p"]
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
    ax[0].scatter(pc[:, 0], pc[:, 1], s=15, c="#cccccc", edgecolors="none", zorder=1)
    sc = ax[0].scatter(cen[:, 0], cen[:, 1], c=bio1, cmap="coolwarm", s=130, edgecolors="k", lw=.5, zorder=3)
    plt.colorbar(sc, ax=ax[0], label="site bio1 (°C)")
    ax[0].set_xlabel("founder PC1 (kinship)"); ax[0].set_ylabel("PC2")
    ax[0].set_title("winning-clade centroid per site, colored by climate\n(climate-ordered spread ⇒ local adaptation)",
                    fontsize=9, loc="left"); ax[0].spines[["top", "right"]].set_visible(False)
    ax[1].scatter(bio1, cen[:, 0], s=70, c=bio1, cmap="coolwarm", edgecolors="k", lw=.4, zorder=3)
    bb1, bb0 = np.polyfit(bio1[fin], cen[fin, 0], 1); xs = np.array([np.nanmin(bio1), np.nanmax(bio1)])
    ax[1].plot(xs, bb0 + bb1 * xs, color="firebrick", lw=1.5, zorder=2)
    ax[1].set_xlabel("site bio1 (°C)"); ax[1].set_ylabel("winning-clade PC1")
    ax[1].set_title(f"does climate drive which clade wins?  r={r1:.2f}, perm p={p1:.3f}", fontsize=10, loc="left")
    ax[1].spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Cross-site winners — {S} sites (matched to multisite GWAS): clade×climate structure  "
                 f"(mean cross-site r={R[iu].mean():.2f})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(f"{H}/cross_site_winners_30{SUFFIX}.png", dpi=150)
    print(f"\n[done] {H}/cross_site_winners_30{SUFFIX}.csv + .json + .png")


if __name__ == "__main__":
    main()

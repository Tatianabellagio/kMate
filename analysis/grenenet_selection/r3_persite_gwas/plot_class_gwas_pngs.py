#!/usr/bin/env python
"""Render ALL Manhattan + QQ plots for the class-split GWAS to disk as PNGs (not embedded in a
notebook, for browsing on demand): multitrait JOINT/GLOBAL/CLIMATE (all 19 bioclim vars + PC1) and
per-site individual scans, for SNP/non-SNP/SV-only. Companion to class_gwas_multitrait.ipynb /
class_gwas_persite.ipynb (which show only the peak-overlap TABLES, no plots, per user request).

Output -> analysis/grenenet_selection/r3_persite_gwas/results/varexp/gwas_plots/{multitrait,persite}_{class}_{name}_{manhattan,qq}.png
Env: basic (matplotlib hangs in `plotting`).
"""
import os, sys, time
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 100, "font.size": 9})

sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection")
import lib

OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r3_persite_gwas/results/varexp"
PLOTDIR = f"{OUT}/gwas_plots"
os.makedirs(PLOTDIR, exist_ok=True)
CHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5"]
CLASSES = ["snp", "nonsnp", "sv"]
COLORS = {"snp": ["#3b4cc0", "#7aa0c4"], "nonsnp": ["#D55E00", "#f2a679"], "sv": ["#009E73", "#66c2a5"]}


def load_class(name):
    return dict(np.load(f"{OUT}/class_gwas_{name}.npz", allow_pickle=True))


def genome_x_df(chrom, pos, p):
    df = pd.DataFrame({"chrom": chrom, "pos": pos, "p": p})
    df = df.sort_values(["chrom", "pos"]).reset_index(drop=True)
    off, centers, x = 0.0, [], np.zeros(len(df))
    for ch in CHROMS:
        mk = (df.chrom == ch).to_numpy()
        if not mk.any():
            continue
        x[mk] = df.pos[mk].to_numpy() + off
        centers.append(off + df.pos[mk].max() / 2)
        off += df.pos[mk].max() * 1.02
    df["x"] = x
    return df, centers


def block_thresholds(blk, p, q=0.05):
    """BLOCK-level Bonferroni + BH-FDR significance thresholds, matching the clq0.9-block overlap
    tables in the notebooks (NOT marker-level). Collapse markers to one lead (min-p) per clq0.9
    block, then:
      - Bonferroni p-threshold = 0.05 / n_blocks
      - BH-FDR p-threshold = the block-level BH critical p (largest lead-p_(k) with p_(k) <= k*q/m),
        or None if no block passes.
    Both are hard p-value cutoffs on the block-lead p-values, and since a block's lead is its
    min-p marker, drawing a horizontal line at either threshold on the MARKER manhattan is exact:
    every marker above the line sits in a significant block, and every significant block has its
    lead above the line. Returns (bonf_p, fdr_p, n_bonf_blocks, n_fdr_blocks)."""
    keep = blk != ""
    nlp = -np.log10(np.clip(p[keep], 1e-300, 1))
    df = pd.DataFrame({"block": blk[keep], "nlp": nlp})
    lead = lib.collapse_to_blocks(df, stat="nlp", block_col="block")
    lp = np.sort(10.0 ** (-lead["nlp"].to_numpy()))
    m = len(lp)
    bonf_p = 0.05 / m
    crit = lp <= (np.arange(1, m + 1) * q / m)
    fdr_p = float(lp[crit][-1]) if crit.any() else None
    return bonf_p, fdr_p, int((lp < bonf_p).sum()), int((lp <= (fdr_p if fdr_p else -1)).sum())


def manhattan(df, centers, colors, title, bonf_p, fdr_p, n_bonf, n_fdr, lam, fname):
    fig, ax = plt.subplots(figsize=(11, 3.6))
    nlp = -np.log10(np.clip(df["p"], 1e-300, 1))
    for i, ch in enumerate(CHROMS):
        mk = (df.chrom == ch).to_numpy()
        ax.scatter(df.x[mk], nlp[mk], s=4, c=colors[i % 2], alpha=.5, edgecolors="none", rasterized=True)
    ax.axhline(-np.log10(bonf_p), color="firebrick", lw=.9, ls="--",
               label=f"Bonferroni 0.05 ({n_bonf} blocks)")
    if fdr_p is not None:
        ax.axhline(-np.log10(fdr_p), color="steelblue", lw=.9, ls=":",
                   label=f"BH FDR q<0.05 ({n_fdr} blocks)")
    else:
        ax.plot([], [], " ", label="BH FDR q<0.05 (0 blocks)")
    ax.set_xticks(centers); ax.set_xticklabels([c.replace("chr", "Chr") for c in CHROMS])
    ax.set_ylabel(r"$-\log_{10}p$")
    ax.set_title(f"{title}  (lambda={lam:.3f})", loc="left", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(fname, dpi=100, bbox_inches="tight")
    plt.close(fig)


def qqplot(p, lam, title, fname):
    p = np.sort(np.clip(p[np.isfinite(p)], 1e-300, 1))
    n = len(p)
    if n == 0:
        return
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(p)
    if n > 100_000:
        idx = np.unique(np.linspace(0, n - 1, 100_000).astype(int))
        exp, obs = exp[idx], obs[idx]
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.scatter(exp, obs, s=4, alpha=.4, edgecolors="none", rasterized=True)
    mx = max(exp.max(), obs.max())
    ax.plot([0, mx], [0, mx], "r--", lw=1)
    ax.set_xlabel("expected -log10(p)"); ax.set_ylabel("observed -log10(p)")
    ax.set_title(f"{title}  (lambda={lam:.3f})", loc="left", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(fname, dpi=100, bbox_inches="tight")
    plt.close(fig)


def climate_p(Z, cache_key, cache, bvec):
    S = Z.shape[1]
    if cache_key not in cache:
        cache[cache_key] = np.corrcoef(Z.T)
    C = cache[cache_key]
    Cinv = np.linalg.pinv(C); one = np.ones(S)
    dg = float(one @ Cinv @ one)
    c0 = (bvec - bvec.mean()) / bvec.std()
    c = c0 - (float(one @ Cinv @ c0) / dg) * one
    dc = float(c @ Cinv @ c)
    z_clim = (Z @ Cinv @ c) / np.sqrt(dc)
    return 2 * stats.norm.sf(np.abs(z_clim))


def block_ids(chrom_lower, pos):
    chrom_cap = np.array([c.replace("chr", "Chr") for c in chrom_lower])
    return lib.assign_clq_blocks(chrom_cap, pos, r2=0.9)


def main():
    t0 = time.time()
    data = {c: load_class(c) for c in CLASSES}
    clim_all = lib.load_climate()
    Ccache = {}
    # clq0.9 block ids per class (same across all contrasts for a class) -- for block-level
    # Bonferroni/FDR threshold lines that MATCH the overlap tables in the notebooks
    blk = {c: block_ids(data[c]["chrom"], data[c]["pos"]) for c in CLASSES}
    print("block ids computed", flush=True)

    print("=== multitrait: JOINT / GLOBAL / CLIMATE x (19 bio vars + PC1) ===", flush=True)
    for cls in CLASSES:
        d = data[cls]
        chrom, pos, sites = d["chrom"], d["pos"], d["sites"]
        clim = clim_all.reindex(sites)[lib.BIO_COLS]
        Bz = (clim - clim.mean()) / clim.std()
        Uu, Ss, _ = np.linalg.svd(Bz.to_numpy(), full_matrices=False)
        pc1 = Uu[:, 0] * Ss[0]
        axes = {f"bio{i}": clim[f"bio{i}"].to_numpy(float) for i in range(1, 20)}
        axes["PC1_allbio"] = pc1

        contrasts = [("JOINT", d["p_joint"]), ("GLOBAL", d["p_global"])]
        for name, bvec in axes.items():
            contrasts.append((f"CLIMATE_{name}", climate_p(d["Z"], cls, Ccache, bvec)))

        for name, p in contrasts:
            lam = lib.lamgc(p)
            bonf_p, fdr_p, n_bonf, n_fdr = block_thresholds(blk[cls], p)
            dfx, centers = genome_x_df(chrom, pos, p)
            manhattan(dfx, centers, COLORS[cls], f"{cls} multitrait {name}",
                      bonf_p, fdr_p, n_bonf, n_fdr, lam,
                      f"{PLOTDIR}/multitrait_{cls}_{name}_manhattan.png")
            qqplot(p, lam, f"{cls} multitrait {name}", f"{PLOTDIR}/multitrait_{cls}_{name}_qq.png")
        print(f"  {cls}: {len(contrasts)} contrasts done ({time.time()-t0:.0f}s)", flush=True)

    print("=== per-site: 31 individual gardens ===", flush=True)
    for cls in CLASSES:
        d = data[cls]
        chrom, pos, Z, sites, bio1 = d["chrom"], d["pos"], d["Z"], d["sites"], d["bio1"]
        for si, site_id in enumerate(sites):
            z = Z[:, si]
            p = 2 * stats.norm.sf(np.abs(z))
            lam = float(np.nanmedian(z ** 2) / stats.chi2.ppf(0.5, 1))
            bonf_p, fdr_p, n_bonf, n_fdr = block_thresholds(blk[cls], p)
            dfx, centers = genome_x_df(chrom, pos, p)
            title = f"{cls} site{int(site_id)} (bio1={bio1[si]:.1f}C)"
            manhattan(dfx, centers, COLORS[cls], title, bonf_p, fdr_p, n_bonf, n_fdr, lam,
                      f"{PLOTDIR}/persite_{cls}_site{int(site_id)}_manhattan.png")
            qqplot(p, lam, title, f"{PLOTDIR}/persite_{cls}_site{int(site_id)}_qq.png")
        print(f"  {cls}: {len(sites)} sites done ({time.time()-t0:.0f}s)", flush=True)

    n_png = len(os.listdir(PLOTDIR))
    print(f"\nwrote {n_png} PNGs -> {PLOTDIR}  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()

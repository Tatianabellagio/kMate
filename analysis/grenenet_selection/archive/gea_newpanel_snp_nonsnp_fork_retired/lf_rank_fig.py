#!/usr/bin/env python
"""Render QQ (raw vs K16-adjusted) + Manhattan (SNP/nonSNP) for LF-rank GEA.
Run in `basic` env. Reads fig_inputs.npz + fig_gif.npz."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/results/lf_rank"
Z = np.load(f"{OUT}/fig_inputs.npz", allow_pickle=True)
G = np.load(f"{OUT}/fig_gif.npz", allow_pickle=True)
CHR = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
CMAP = {c: col for c, col in zip(CHR, ["#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677"])}


def qq_pts(p, n_max=60000):
    p = p[np.isfinite(p)]
    p = np.clip(p, 1e-300, 1)
    obs = -np.log10(np.sort(p))
    exp = -np.log10((np.arange(1, p.size + 1) - 0.5) / p.size)
    if p.size > n_max:                       # thin for plotting
        idx = np.unique(np.linspace(0, p.size - 1, n_max).astype(int))
        obs, exp = obs[idx], exp[idx]
    return exp, obs


fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.32, wspace=0.22)

# ---- QQ panels ----
for j, cls in enumerate(["snp", "nonsnp"]):
    ax = fig.add_subplot(gs[0, j])
    ex_u, ob_u = qq_pts(Z[f"{cls}_punadj"])
    ex_a, ob_a = qq_pts(Z[f"{cls}_p"])
    gu, ga = float(G[f"{cls}_unadj"]), float(G[f"{cls}_adj"])
    ax.scatter(ex_u, ob_u, s=6, c="#BBBBBB", label=f"raw (K=0)  GIF={gu:.2f}", rasterized=True)
    ax.scatter(ex_a, ob_a, s=6, c="#EE6677" if cls == "snp" else "#228833",
               label=f"K=16 adj  GIF={ga:.2f}", rasterized=True)
    lim = max(ex_u.max(), ob_u.max())
    ax.plot([0, lim], [0, lim], "k--", lw=0.8)
    ax.set_title(f"{cls.upper()}  QQ (partial Spearman vs bio1)", fontsize=11)
    ax.set_xlabel("expected -log10 p"); ax.set_ylabel("observed -log10 p")
    ax.legend(fontsize=8, loc="upper left", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

# ---- Manhattan panels (K16-adjusted) ----
def manhattan(ax, cls, title):
    chrom = Z[f"{cls}_chrom"].astype(str)
    pos = Z[f"{cls}_pos"].astype(float)
    p = np.clip(Z[f"{cls}_p"].astype(float), 1e-300, 1)
    q = Z[f"{cls}_q"].astype(float)
    m = int(np.isfinite(Z[f"{cls}_p"]).sum())
    bonf = 0.05 / m
    off = 0.0; ticks = []
    for c in CHR:
        sel = chrom == c
        if not sel.any():
            continue
        x = pos[sel] + off
        nlp = -np.log10(p[sel])
        # thin non-significant points for size
        sig = nlp > 2
        keep = sig | (np.random.default_rng(0).random(sel.sum()) < 0.06)
        ax.scatter(x[keep], nlp[keep], s=3, c=CMAP[c], rasterized=True)
        ticks.append((off + pos[sel].max() / 2, c))
        off += pos[sel].max() + 5e6
    ax.axhline(-np.log10(bonf), color="k", ls="--", lw=0.8,
               label=f"Bonferroni {bonf:.1e}")
    # FDR q<0.05 threshold = largest p with q<0.05
    fdr_p = Z[f"{cls}_p"][np.isfinite(q) & (q < 0.05)]
    if fdr_p.size:
        ax.axhline(-np.log10(fdr_p.max()), color="grey", ls=":", lw=0.8,
                   label="FDR q<0.05")
    ga = float(G[f"{cls}_adj"])
    ax.set_title(f"{title}   (K=16, GIF={ga:.2f}, Bonf hits={int((p<bonf).sum())})",
                 fontsize=11)
    ax.set_xticks([t[0] for t in ticks]); ax.set_xticklabels([t[1] for t in ticks])
    ax.set_ylabel("-log10 p (adj)")
    ax.legend(fontsize=8, loc="upper right", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

manhattan(fig.add_subplot(gs[1, 0]), "snp", "SNP")
manhattan(fig.add_subplot(gs[1, 1]), "nonsnp", "non-SNP (indel+SV)")

rho = float(Z["block_rho"])
fig.suptitle(f"Latent-factor (K=16) adjusted partial-rank climate-GEA vs bio1  |  "
             f"SNP-vs-nonSNP block-peak Spearman rho={rho:.2f}", fontsize=12)
fig.savefig(f"{OUT}/lf_rank_qq_manhattan.png", dpi=150, bbox_inches="tight")
print(f"saved {OUT}/lf_rank_qq_manhattan.png")

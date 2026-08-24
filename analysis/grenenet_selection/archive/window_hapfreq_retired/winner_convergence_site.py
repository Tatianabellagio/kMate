#!/usr/bin/env python
"""Are site-4 winners one clade or convergently diverse? + convergent-allele test (kinship-corrected).

The adaptive signal = alleles the WINNING ecotypes convergently share, resolvable wherever the
founder PANEL decouples them from hitchhikers (panel LD), NOT experimental recombination. This:
  (1) WINNER DIVERSITY: founder GRM (kinship) PCA colored by selection response; mean pairwise
      kinship among winners vs background. If winners are ONE clade, everything they share is
      confounded with kinship (can't fine-map). If they SPAN the tree (convergent), shared
      adaptive alleles separate from clade-wide hitchhikers.
  (2) CONVERGENT-ALLELE test: kinship-corrected association of each haplotype with the founder
      selection response (EMMAX). An allele surviving kinship = shared by winners beyond their
      relatedness = the convergent candidate. For the top survivors, measure CARRIER SPREAD in
      kinship space (do the winning-allele carriers come from different clades?).

Env: kmate. SITE via env. Reuses founder_genotype.build_genotype / emma_reml_delta and the
per-founder selection from site<ID>_ecotype_selection.csv.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from founder_genotype import build_genotype, emma_reml_delta

H = "analysis/grenenet_gea/archive/window_hapfreq_retired/hapfreq"
SITE = int(os.environ.get("SITE", 4))
MAC = 3


def main():
    G, founders, reg = build_genotype()
    ph = pd.read_csv(f"{H}/site{SITE}_ecotype_selection.csv")
    ph["founder"] = ph.founder.astype(str)
    smap = dict(zip(ph.founder, ph.s_bar)); dmap = dict(zip(ph.founder, ph.d_freq))
    qmap = dict(zip(ph.founder, ph.q))
    fi = np.array([i for i, f in enumerate(founders) if f in smap])
    y = np.array([smap[founders[i]] for i in fi])              # selection coeff per founder
    dfreq = np.array([dmap[founders[i]] for i in fi])
    qv = np.array([qmap[founders[i]] for i in fi])
    fnames = np.array([founders[i] for i in fi])
    G = G[fi]; n = len(y)
    yq = stats.norm.ppf((stats.rankdata(y) - 0.5) / n)         # quantile-normalised

    # markers + GRM (kinship)
    cnt = G.sum(0); poly = (cnt >= MAC) & (cnt <= n - MAC)
    Gp = G[:, poly]; M = Gp.shape[1]
    p = Gp.mean(0); Z = (Gp - p) / np.sqrt(p * (1 - p) + 1e-9)
    K = (Z @ Z.T) / M; K = (K + K.T) / 2

    # ---- (1) winner diversity: PCA of K, colored by selection ----
    lam, V = np.linalg.eigh(K)
    pc = V[:, ::-1][:, :2] * np.sqrt(np.clip(lam[::-1][:2], 0, None))   # PC scores
    winners = (qv < 0.05) & (dfreq > 0)                       # significant risers
    losers = (qv < 0.05) & (dfreq < 0)
    # mean pairwise kinship among winners vs random
    def meank(idx):
        if idx.sum() < 2:
            return np.nan
        sub = K[np.ix_(idx, idx)]; iu = np.triu_indices(idx.sum(), 1)
        return sub[iu].mean()
    rng = np.random.RandomState(0)
    rand_k = np.mean([meank(np.isin(np.arange(n), rng.choice(n, winners.sum(), replace=False)))
                      for _ in range(200)])
    print(f"site {SITE}: {n} founders, {M:,} markers | winners(q<0.05,rose)={int(winners.sum())}, "
          f"losers={int(losers.sum())}")
    print(f"  mean pairwise kinship among winners = {meank(winners):.3f}  vs random sets = {rand_k:.3f}")
    print(f"  -> winners {'CLUSTER (one clade, hard to fine-map)' if meank(winners) > 2*rand_k else 'SPREAD across the tree (convergent, fine-mappable)'}")

    # ---- (2) convergent-allele test: kinship-corrected EMMAX on s_bar ----
    lamK, U = np.linalg.eigh(K); lamK = np.clip(lamK, 1e-9, None)
    yr = U.T @ yq; Xr = U.T @ np.ones((n, 1))
    delta = emma_reml_delta(yr, Xr, lamK); w = 1.0 / (lamK + delta)
    Br = U.T @ Gp; a = Xr[:, 0]
    Saa = (w * a * a).sum(); Say = (w * a * yr).sum()
    Sab = ((w * a)[None, :] @ Br).ravel(); Sbb = (w[:, None] * Br**2).sum(0)
    Sby = ((w * yr)[None, :] @ Br).ravel()
    det = Saa * Sbb - Sab**2
    beta = (Saa * Sby - Sab * Say) / np.clip(det, 1e-30, None)
    alpha = (Sbb * Say - Sab * Sby) / np.clip(det, 1e-30, None)
    rss = (w * yr * yr).sum() - alpha * Say - beta * Sby
    s2 = np.clip(rss, 1e-30, None) / (n - 2)
    se = np.sqrt(np.clip(s2 * Saa / np.clip(det, 1e-30, None), 1e-30, None))
    t = beta / se; pk = 2 * stats.t.sf(np.abs(t), df=n - 2)
    lam_gc = np.median(stats.chi2.isf(np.clip(pk, 1e-300, 1), 1)) / stats.chi2.ppf(0.5, 1)

    regp = reg.iloc[np.where(poly)[0]][["chrom", "start", "end"]].reset_index(drop=True)
    regp["beta"] = beta; regp["p_kin"] = pk; regp["mac"] = cnt[poly]
    m = len(regp); o = regp.p_kin.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate((regp.p_kin.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    regp["q"] = np.clip(q, 0, 1)

    # carrier SPREAD: for each top allele, are winning-allele carriers spread in kinship/PC space?
    top = regp.sort_values("p_kin").head(15).copy()
    spreads = []
    for idx, r in top.iterrows():
        gcol = Gp[:, idx]                                      # carriers of this haplotype
        carr = gcol > 0
        # spread of carriers in PC1 (sd) relative to all; and mean kinship among carriers
        sp = pc[carr, 0].std() if carr.sum() > 1 else np.nan
        mk = meank(carr)
        spreads.append((sp, mk, int(carr.sum())))
    top["carrier_pc1_sd"] = [s[0] for s in spreads]
    top["carrier_meank"] = [s[1] for s in spreads]
    top["n_carry"] = [s[2] for s in spreads]
    top["clade_spread"] = top.carrier_meank < 2 * rand_k       # carriers NOT one clade = convergent

    print(f"\n  kinship-corrected EMMAX on s_bar: lambda_GC={lam_gc:.2f}, "
          f"Bonferroni hits {int((regp.p_kin<0.05/m).sum())}, FDR q<0.05 {int((regp.q<0.05).sum())}")
    print(f"  overall PC1 sd (all founders) = {pc[:,0].std():.2f}  (carrier sd >~ this = convergent)")
    print("\n  TOP kinship-corrected alleles (do winning carriers span clades?):")
    print(top[["chrom", "start", "end", "mac", "beta", "p_kin", "q",
               "n_carry", "carrier_meank", "clade_spread"]].to_string(index=False))

    # ---- figure: founder PCA colored by selection, winners marked ----
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    sc = ax[0].scatter(pc[:, 0], pc[:, 1], c=dfreq, cmap="coolwarm",
                       vmin=-np.abs(dfreq).max(), vmax=np.abs(dfreq).max(), s=35, edgecolors="k", lw=.3)
    ax[0].scatter(pc[winners, 0], pc[winners, 1], s=120, facecolors="none", edgecolors="red", lw=1.5, label="winner (q<.05)")
    plt.colorbar(sc, ax=ax[0], label="Δfreq (winning)")
    ax[0].set_xlabel("founder PC1 (kinship)"); ax[0].set_ylabel("PC2")
    ax[0].set_title(f"Site {SITE}: are winners one clade or convergent?", fontsize=10, loc="left")
    ax[0].legend(fontsize=8, frameon=False); ax[0].spines[["top", "right"]].set_visible(False)
    # QQ of the kinship-corrected convergent test
    pv = np.sort(regp.p_kin.values); exp = -np.log10((np.arange(1, m + 1) - .5) / m)
    ax[1].scatter(exp, -np.log10(np.clip(pv, 1e-300, 1)), s=5, c="#34495e", alpha=.5, edgecolors="none")
    lim = exp.max(); ax[1].plot([0, lim], [0, lim], color="firebrick", lw=1)
    ax[1].set_xlabel(r"expected $-\log_{10}p$"); ax[1].set_ylabel(r"observed")
    ax[1].set_title(f"convergent-allele test (kinship-corrected) λ={lam_gc:.2f}", fontsize=10, loc="left")
    ax[1].spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{H}/site{SITE}_winner_convergence.png", dpi=150)
    regp.sort_values("p_kin").to_csv(f"{H}/site{SITE}_convergent_alleles.csv", index=False)
    print(f"\n[done] {H}/site{SITE}_winner_convergence.png + convergent_alleles.csv")


if __name__ == "__main__":
    main()

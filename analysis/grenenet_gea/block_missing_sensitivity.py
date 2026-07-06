#!/usr/bin/env python
"""Did missing-data imputation distort the haplotype clusters? (Chr1, CLQ0.9 + gate2)

For each multi-hap block (n_eff>2), cluster the 231 founders with the SAME algorithm
(average-linkage agglomerative, cut at k = the production xmeans k) under TWO missing
treatments, isolating the imputation effect:
  - IMPUTED distance : Hamming on geno=(raw>=0.5), i.e. missing -> major allele
    (what the production clusters use).
  - MASKED distance  : Hamming computed ONLY over co-called sites (var_called),
    using the real allele, ignoring missing.
Compare assignments via adjusted Rand index (ARI). ARI~1 => imputation irrelevant;
ARI drops in missing-heavy blocks => imputation reshaped the clusters there.
Also report ARI(xmeans vs agglom-imputed) for algorithm-sensitivity context.

Run in kmate env. Output: chr1_missing_sensitivity.csv
"""
import os, sys, argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import adjusted_rand_score
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from block_haplotype_counts import hap_counts
from block_cluster_pc1ve import cluster_founders   # xmeans (shim applied on import)

BR = "results/grenenet_gea/blocks_recompute"
MAF, MINCF, GATE = 0.05, 0.5, 2.0


def load(pref):
    vp = sp.load_npz(f"{pref}.var_pa.npz").toarray().astype(np.float64)
    vc = sp.load_npz(f"{pref}.var_called.npz").toarray().astype(np.float64)
    pos = np.load(f"{pref}.meta.npz", allow_pickle=True)["pos"].astype(np.int64)
    F = vp.shape[0]
    n_called = vc.sum(0); n_alt = vp.sum(0)
    af = np.divide(n_alt, n_called, out=np.full_like(n_alt, np.nan), where=n_called > 0)
    keep = (af > MAF) & (af < 1 - MAF) & (n_called >= MINCF * F)
    order = np.where(keep)[0]
    seen = set(); uniq = []
    for c in order:
        p = int(pos[c])
        if p not in seen:
            seen.add(p); uniq.append(c)
    uniq = np.array(uniq)
    vp_u = vp[:, uniq]; vc_u = vc[:, uniq]; afc = af[uniq]
    raw = vp_u.copy(); miss = vc_u == 0
    raw[miss] = np.repeat(afc[None, :], F, axis=0)[miss]
    geno_imp = (raw >= 0.5).astype(np.int8)          # missing -> major (production)
    return geno_imp, vp_u.astype(np.int8), vc_u.astype(np.int8), np.array([int(p) for p in pos[uniq]])


def dist_imputed(G):
    """Hamming distance matrix on imputed binary genotypes (F x m)."""
    m = G.shape[1]
    X = G.astype(np.float64)
    diff = X @ (1 - X).T + (1 - X) @ X.T
    return diff / max(m, 1)


def dist_masked(Vp, Vc):
    """Hamming over CO-CALLED sites only. Vp,Vc: F x m (0/1)."""
    X = (Vp * Vc).astype(np.float64)        # called & alt
    M = Vc.astype(np.float64)               # called
    C = M @ M.T                             # co-called counts
    diff = X @ (M - X).T + (M - X) @ X.T    # differing among co-called
    with np.errstate(invalid="ignore", divide="ignore"):
        D = diff / C
    D[C == 0] = 1.0                         # no shared info -> maximally distant
    np.fill_diagonal(D, 0.0)
    return 0.5 * (D + D.T)                  # symmetrize tiny fp asymmetry


def agglom(D, k):
    if k <= 1 or D.shape[0] <= k:
        return np.zeros(D.shape[0], dtype=int) if k <= 1 else np.arange(D.shape[0])
    Z = linkage(squareform(D, checks=False), method="average")
    return fcluster(Z, t=k, criterion="maxclust") - 1


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--chrom", default="Chr1")
    a = ap.parse_args(); chrlc = a.chrom.lower()
    print(f"[{a.chrom}] loading ...", flush=True)
    Gimp, Vp, Vc, positions = load(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}")
    bl = pd.read_csv(f"{BR}/{chrlc}_clq0.9_blocks_clq0.9.tsv", sep="\t")

    rows = []
    for bi, (s, e) in enumerate(zip(bl.start_pos, bl.end_pos)):
        lo = int(np.searchsorted(positions, s)); hi = int(np.searchsorted(positions, e, side="right"))
        if hi - lo < 2:
            continue
        gi = Gimp[:, lo:hi]; vp = Vp[:, lo:hi]; vc = Vc[:, lo:hi]
        nd, neff, ng = hap_counts(gi)
        if neff <= GATE:
            continue
        try:
            xm = cluster_founders(gi)[0]                 # production xmeans labels
        except Exception:
            continue
        k = len(np.unique(xm))
        if k < 2:
            continue
        li = agglom(dist_imputed(gi), k)
        lm = agglom(dist_masked(vp, vc), k)
        rows.append(dict(block=bi, n_eff=round(float(neff), 3), nvar=hi - lo, k=k,
                         mean_callrate=round(float(vc.mean()), 4),
                         ari_imp_vs_mask=round(adjusted_rand_score(li, lm), 4),
                         ari_xm_vs_imp=round(adjusted_rand_score(xm, li), 4)))
        if (bi + 1) % 4000 == 0:
            print(f"  ...{bi+1}/{len(bl)}", flush=True)

    df = pd.DataFrame(rows)
    out = f"{BR}/{chrlc}_missing_sensitivity.csv"
    df.to_csv(out, index=False)
    print(f"\nwrote {len(df)} multi-hap blocks -> {out}")
    print(f"ARI imputed-vs-masked: median {df.ari_imp_vs_mask.median():.3f}  "
          f"mean {df.ari_imp_vs_mask.mean():.3f}  %>=0.8 {100*(df.ari_imp_vs_mask>=0.8).mean():.0f}%")
    print("by block mean call-rate:")
    for lo, hi in [(0, .8), (.8, .9), (.9, .95), (.95, 1.01)]:
        b = df[(df.mean_callrate >= lo) & (df.mean_callrate < hi)]
        if len(b):
            print(f"  call-rate {lo:.2f}-{hi:.2f}: n={len(b):5d}  median ARI {b.ari_imp_vs_mask.median():.3f}")


if __name__ == "__main__":
    main()

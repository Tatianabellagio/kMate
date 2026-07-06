#!/usr/bin/env python3
"""Selection patterns BEYOND the linear climate gradient -- generalizes the founder-GWAS
CLIMATE contrast (continuous bioPC1 gradient, found null: best q=0.25) two ways, to test
whether SV-carrying haploblocks are selected in a pattern a straight-line gradient test is
structurally blind to.

(1) CATEGORICAL: k-means-cluster the 30 gardens on standardized bio1-19 (K by silhouette,
    min cluster size 3). EMPIRICAL FINDING (see the printed silhouette sweep): this cohort's
    30 gardens do NOT support >=3 well-separated climate archetypes -- K>=3 always peels off
    a size-1 outlier garden rather than splitting the main group, so K=2 is the only
    well-supported split. At K=2 the "between-cluster" and "each cluster vs the rest" tests
    are the SAME 1-df contrast (z_cluster1 = -z_cluster0 exactly, since cluster1 = 1-cluster0
    after GLOBAL-orthogonalization) -- reported once, not as independent looks. Its relation
    to the existing linear CLIMATE contrast is reported via their cross-marker correlation
    (this split turns out to load more on precipitation than temperature, so it is NOT simply
    a coarsened version of the bio1/bioPC1 gradient -- see the printed corr).
(2) QUADRATIC (the literal test of "selected up in BOTH hot and cold, flat in the middle"):
    a 1-df U-shape contrast on the SAME climate axis as CLIMATE (squared, then C-orthogonalized
    against both GLOBAL and the linear CLIMATE contrast) -- this is what actually targets
    non-monotonic-in-temperature selection; a k-means split cannot, by construction, isolate
    "both extremes vs the middle" along one axis as cleanly as a polynomial term can.

Both are exact extensions of the same orthogonal chi-square decomposition already used for
GLOBAL/CLIMATE in founder_gwas_multisite.py (JOINT = GLOBAL + CLIMATE + ... + residual), so
they nest cleanly and reuse Z/C from multisite_founder_gwas_clq90_pc1.npz -- no new GWAS run.

Then: size-matched permutation SV enrichment (same design as sv_enrichment.py) on top blocks
for BETWEEN_CLUSTER and QUADRATIC, plus pairwise Jaccard overlap of the top blocks flagged by
JOINT / CLIMATE(linear) / BETWEEN_CLUSTER / QUADRATIC -- do these surface the same loci or
genuinely different ones?

Run in `basic` env from the kmate repo root. Deterministic (seed=0).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import lib

OUT = Path("results/grenenet_gea/sv_adaptive")
NPERM = 10000
TOP = [0.005, 0.01, 0.02]
EDGES = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 250, 10**9]
K_RANGE = range(2, 7)
MIN_CLUSTER_N = 3

raw_all = lib.multisite_gwas_raw("clq90_pc1")
print(f"loaded {raw_all['M']:,} markers x {raw_all['S']} sites "
      f"(self-check OK, JOINT lambda={lib.lamgc(raw_all['p_joint']):.4f})")
MIN_MAF = 0.01
raw = lib.maf_filter(raw_all, MIN_MAF)
Z, Cinv, sites = raw["Z"], raw["Cinv"], raw["sites"]
one, dg, S, M = raw["one"], raw["dg"], raw["S"], raw["M"]
z_global, chi2_joint, p_joint = raw["z_global"], raw["chi2_joint"], raw["p_joint"]
z_clim = raw["z_clim"]                            # linear CLIMATE, axis = bioPC1 for this tag
print(f"MAF filter (founder MAF>={MIN_MAF:.0%}): kept {M:,}/{raw_all['M']:,} markers, "
      f"new JOINT lambda={lib.lamgc(p_joint):.4f}")

# ---- 1. climate clusters (k-means on standardized bio1-19) ----
clim19 = lib.load_climate().reindex(sites)[lib.BIO_COLS]
assert clim19.notna().all().all(), "missing bio1-19 for a cohort site"
Xc = StandardScaler().fit_transform(clim19.to_numpy())

sil, labels_by_k = {}, {}
for K in K_RANGE:
    lab = KMeans(K, n_init=10, random_state=0).fit_predict(Xc)
    if np.bincount(lab).min() < MIN_CLUSTER_N:
        continue
    labels_by_k[K] = lab
    sil[K] = silhouette_score(Xc, lab)
K_best = max(sil, key=sil.get)
labels = labels_by_k[K_best]
print("silhouette by K (K's with a <3-site cluster are dropped):",
      {k: round(v, 3) for k, v in sil.items()}, "-> chosen K =", K_best)
for k in range(K_best):
    s = sites[labels == k]
    b1, b12 = clim19.loc[s, "bio1"], clim19.loc[s, "bio12"]
    print(f"  cluster {k}: {len(s)} sites {list(s)}  bio1 {b1.mean():.1f}+-{b1.std():.1f}C  "
          f"bio12(precip) {b12.mean():.0f}+-{b12.std():.0f}mm")

pd.DataFrame({"site": sites, "cluster": labels, **{c: clim19[c].to_numpy() for c in lib.BIO_COLS}}
            ).to_csv(OUT / "climate_clusters.csv", index=False)

# ---- 2. BETWEEN-CLUSTER (K-1 df) + per-cluster one-vs-rest contrasts ----
D = np.zeros((S, K_best))
for k in range(K_best):
    D[:, k] = (labels == k).astype(float)
A = D.T @ Cinv @ D
U = Z @ Cinv @ D
Q_K = np.einsum("mk,kl,ml->m", U, np.linalg.pinv(A), U)
df_between = K_best - 1
Q_between = np.clip(Q_K - z_global ** 2, 0, None)
p_between = stats.chi2.sf(Q_between, df_between)
nest_ok = (Q_between <= chi2_joint + 1e-6).mean()
assert nest_ok > 0.999, f"between-cluster stat exceeds JOINT for {1 - nest_ok:.1%} of markers"
print(f"BETWEEN_CLUSTER ({df_between} df): lambda={lib.lamgc(p_between):.3f}, "
      f"q<0.05 n={int((lib.bh(p_between) < 0.05).sum()):,}")

zk = np.zeros((M, K_best))
for k in range(K_best):
    dk = D[:, k]
    dk_orth = dk - (float(one @ Cinv @ dk) / dg) * one
    zk[:, k] = (Z @ Cinv @ dk_orth) / np.sqrt(float(dk_orth @ Cinv @ dk_orth))
pk = 2 * stats.norm.sf(np.abs(zk))
if K_best == 2:
    corr_c0_clim = np.corrcoef(zk[:, 0], z_clim)[0, 1]
    print(f"  NOTE K=2: cluster-0 and cluster-1 one-vs-rest are the SAME test as "
          f"BETWEEN_CLUSTER (z_c1=-z_c0 exactly) -- 1 independent look, not 3.")
    print(f"  corr(cluster-0 contrast, linear CLIMATE) = {corr_c0_clim:+.3f} "
          f"({'loads mostly on the same axis as CLIMATE' if abs(corr_c0_clim) > 0.7 else 'materially different axis from CLIMATE (see bio12 spread above)'})")
else:
    for k in range(K_best):
        print(f"  cluster {k} one-vs-rest: lambda={lib.lamgc(pk[:, k]):.3f}, "
              f"q<0.05 n={int((lib.bh(pk[:, k]) < 0.05).sum()):,}")

# ---- 3. QUADRATIC ("U-shape") contrast on the SAME axis as CLIMATE -- orthogonal to GLOBAL
# and to the linear CLIMATE contrast, so it isolates "selected at the climate extremes
# (either/both), not on a straight line" -- the literal test of the up-in-hot-AND-cold example.
bio1 = raw["bio1"]                                # bioPC1 composite for this tag
c0 = (bio1 - bio1.mean()) / bio1.std()
c = c0 - (float(one @ Cinv @ c0) / dg) * one
dc = float(c @ Cinv @ c)
assert np.allclose((Z @ Cinv @ c) / np.sqrt(dc), z_clim), "c reconstruction != saved z_clim"
q_raw = c0 ** 2
q1 = q_raw - (float(one @ Cinv @ q_raw) / dg) * one
q1 = q1 - (float(c @ Cinv @ q1) / dc) * c
dqv = float(q1 @ Cinv @ q1)
z_quad = (Z @ Cinv @ q1) / np.sqrt(dqv)
p_quad = 2 * stats.norm.sf(np.abs(z_quad))
p_clim = 2 * stats.norm.sf(np.abs(z_clim))
print(f"\nQUADRATIC (U-shape, 1 df, orthogonal to GLOBAL+CLIMATE): "
      f"lambda={lib.lamgc(p_quad):.3f}, q<0.05 n={int((lib.bh(p_quad) < 0.05).sum()):,}, "
      f"corr with linear CLIMATE z = {np.corrcoef(z_quad, z_clim)[0, 1]:+.3f} (sanity ~0)")

# ---- 4. collapse to block level (min p per unit, matches sv_enrichment.py convention) ----
mk = pd.DataFrame({"unit": raw["unit"], "p_joint": p_joint, "p_clim": p_clim,
                   "p_between": p_between, "p_quad": p_quad})
for k in range(K_best):
    mk[f"p_c{k}"] = pk[:, k]
gc = mk.groupby("unit", as_index=False).agg({c: "min" for c in mk.columns if c != "unit"})

L = pd.read_csv(OUT / "sv_landscape_clq0.9.csv")
df = gc.merge(L, left_on="unit", right_on="block_id", how="inner").reset_index(drop=True)
df["bin"] = pd.cut(df.n_kept, EDGES, right=False, labels=False)
df = df[df.bin.notna()].reset_index(drop=True)
df["bin"] = df["bin"].astype(int)
print(f"\n{len(df):,} blocks joined to SV landscape ({df.has_sv.sum():,} carry an SV)")

# ---- 5. size-matched SV enrichment: BETWEEN_CLUSTER + per-cluster + QUADRATIC ----
bins = df["bin"].to_numpy()
rows = []
contrasts = ([("BETWEEN_CLUSTER", "p_between")]
            + ([] if K_best == 2 else [(f"CLUSTER_{k}", f"p_c{k}") for k in range(K_best)])
            + [("QUADRATIC", "p_quad")])
for name, pc in contrasts:
    order = df[pc].to_numpy().argsort()
    for frac in TOP:
        n = int(round(frac * len(df)))
        idx = order[:n]
        for metric in ["has_sv", "sv_frac"]:
            vals = df[metric].to_numpy()
            obs, null, ratio, p = lib.matched_perm_test(vals, bins, idx, NPERM)
            rows.append(dict(contrast=name, top=f"{frac:.1%}", n_top=n, metric=metric,
                             observed=round(obs, 4), null=round(null, 4),
                             enrich=round(ratio, 3), p_perm=round(p, 4)))
res = pd.DataFrame(rows)
res.to_csv(OUT / "climate_cluster_enrichment.csv", index=False)
print("\n== size-matched SV enrichment (BETWEEN_CLUSTER + QUADRATIC" +
      ("" if K_best == 2 else " + per-cluster") + ") ==")
pd.set_option("display.width", 160, "display.max_rows", 200)
print(res.to_string(index=False))

# ---- 6. specificity: do JOINT / CLIMATE / BETWEEN_CLUSTER / QUADRATIC flag the SAME loci? ----
def top_units(pc, frac, pool):
    n = max(int(round(frac * len(pool))), 1)
    return set(pool.nsmallest(n, pc).unit)

tests = ["p_joint", "p_clim", "p_between", "p_quad"]
jac_rows = []
for i in range(len(tests)):
    for j in range(i + 1, len(tests)):
        for tag, pool, frac in [("all_blocks", df, 0.01), ("sv_blocks", df[df.has_sv == 1], 0.10)]:
            a, b = top_units(tests[i], frac, pool), top_units(tests[j], frac, pool)
            u = a | b
            jac_rows.append(dict(test_a=tests[i], test_b=tests[j], subset=tag,
                                 n_a=len(a), n_b=len(b), n_overlap=len(a & b),
                                 jaccard=len(a & b) / len(u) if u else np.nan))
jac = pd.DataFrame(jac_rows)
jac.to_csv(OUT / "climate_cluster_specificity.csv", index=False)
print("\n== pairwise top-block overlap across tests (low jaccard = flags different loci) ==")
print(jac.to_string(index=False))

# per-block cluster/quadratic loadings (for plotting / gene lookups) -- the block's own most-
# significant-anywhere marker, so z_* values reported together are self-consistent
lead_idx = mk.groupby("unit")["p_joint"].idxmin().to_numpy()
zk_df = pd.DataFrame({"unit": raw["unit"], "z_quad": z_quad,
                      **{f"z_c{k}": zk[:, k] for k in range(K_best)}})
lead = zk_df.loc[lead_idx].merge(
    df[["unit", "chrom", "start_pos", "end_pos", "has_sv", "n_kept"]], on="unit")
lead.to_csv(OUT / "climate_cluster_blocks.csv", index=False)

print(f"\n[done] K={K_best} -> {OUT}/climate_clusters.csv, climate_cluster_enrichment.csv, "
      f"climate_cluster_specificity.csv, climate_cluster_blocks.csv")

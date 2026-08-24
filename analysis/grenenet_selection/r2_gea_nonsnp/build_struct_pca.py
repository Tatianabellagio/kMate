#!/usr/bin/env python
"""Estimate population-structure axes (for LFMM K) from SNPs, per generation.

Their old K=3 came from Tracy-Widom on a 231-sample LD-pruned matrix — not our
gen-3 pools and TW flags every PC as significant anyway. So we recompute K for the
ACTUAL gen-3 pools, from SNPs (denser / more neutral than the SVs we test):
  - flower-weighted gen-`g` pool allele frequency for a MAF>=0.05, position-thinned
    SNP subset (LD-pruning proxy),
  - PCA of the [pools x SNP] frequency matrix -> eigenvalues / variance explained
    (scree elbow = candidate K),
  - correlation of each PC with bio1 (collinearity warning: a climate-tracking PC
    means structure correction would erase real signal).

Outputs (--out): struct_pca_gen{g}.npz  (eigvals, var_explained, pc_scores [pools x 20],
pc_bio1_corr, pools) and prints the scree.
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", type=int, default=3)
    ap.add_argument("--n-snps", type=int, default=150000, help="thinned subset size")
    ap.add_argument("--maf", type=float, default=0.05)
    ap.add_argument("--out", default=f"{lib.GEA}/structure")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    store = lib.AF_STORE
    p0 = np.load(f"{store}/p0_snp.npy")
    n_snp = len(p0)
    maf = np.minimum(p0, 1 - p0)
    common = np.where(maf >= args.maf)[0]
    # position-thin to ~n_snps (even stride over the MAF-passing SNPs)
    stride = max(1, len(common) // args.n_snps)
    sel = common[::stride]
    print(f"SNPs: {n_snp:,} | MAF>={args.maf}: {len(common):,} | thinned subset: {len(sel):,}",
          flush=True)

    pt = lib.pool_table(store)
    pt = pt[pt.generation == args.gen]
    pools = sorted(pt.pool.unique(),
                   key=lambda p: tuple(int(x) for x in p.split("_")))
    # flower-weighted pool freq over the SNP subset
    F = np.zeros((len(pools), len(sel)))
    for i, pool in enumerate(pools):
        mem = pt[pt.pool == pool]
        acc = np.zeros(len(sel)); wsum = 0.0
        for s, w in zip(mem.sampleid, mem.flowerscollected.astype(float)):
            v = lib.decode_af(np.load(f"{store}/af_snp/{s}.npy")[sel])
            ok = np.isfinite(v)
            acc[ok] += w * v[ok]; wsum += w
        F[i] = acc / wsum
        if (i + 1) % 50 == 0:
            print(f"  pooled {i+1}/{len(pools)}", flush=True)

    # PCA on standardized SNP frequencies (pools x SNPs)
    X = F - F.mean(0)
    sd = X.std(0); sd[sd == 0] = 1
    X = X / sd
    # SVD: eigenvalues of the pool covariance
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    eig = S ** 2
    ve = eig / eig.sum() * 100
    scores = U * S                                  # pool PC scores [pools x npc]
    meta = pd.read_csv(f"{lib.GEA}/pool_matrices/pool_gen{args.gen}_nonsnp.meta.csv")
    meta = meta.set_index("pool").loc[pools]
    bio1 = meta.bio1.to_numpy()
    npc = min(20, scores.shape[1])
    pc_bio1 = np.array([np.corrcoef(scores[:, k], bio1)[0, 1] for k in range(npc)])

    np.savez(f"{args.out}/struct_pca_gen{args.gen}.npz",
             eigvals=eig, var_explained=ve, pc_scores=scores[:, :npc],
             pc_bio1_corr=pc_bio1, pools=np.array(pools),
             n_snps=len(sel))
    print(f"\nscree (var explained %): {np.round(ve[:12], 2)}")
    print(f"cumulative %:            {np.round(np.cumsum(ve[:12]), 1)}")
    print(f"|corr(PC, bio1)| top:    {np.round(np.abs(pc_bio1[:8]), 2)}")
    bad = np.where(np.abs(pc_bio1[:8]) > 0.5)[0]
    if len(bad):
        print(f"  WARNING: PC{(bad+1).tolist()} track climate (|r|>0.5) — "
              "correcting on them would remove real signal")
    print(f"-> {args.out}/struct_pca_gen{args.gen}.npz")


if __name__ == "__main__":
    main()

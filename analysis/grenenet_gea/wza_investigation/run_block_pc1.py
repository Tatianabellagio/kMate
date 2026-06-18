#!/usr/bin/env python
"""Density-invariant block GEA: PC1 of each LD block's AF matrix vs climate.

For each LD block we have an allele-frequency matrix [pools x block-SNPs]. Because
within-block SNPs are in LD they covary across pools, so almost all across-pool
variation collapses onto one axis (the dominant haplotype-frequency axis). We take
that axis (PC1, one score per pool) and Kendall-correlate it with climate -> ONE
test per block, invariant to how many SNPs tag the block. This sidesteps WZA's
#SNPs-as-LD-proxy entirely: a 16-SNP block and a 1168-SNP block each yield exactly
one test, so a block cannot win just by having more (correlated) SNPs.

Choices (see discussion): NaN imputed with per-SNP across-pool mean; SNPs scaled to
unit variance (correlation-PCA); PC1 variance-explained reported per block.

Output: pc1_{class}_gen{g}_{climate}.csv  per block:
  block, chrom, pos, n_snps, pc1_var_explained, tau, pval   (-> BH downstream)
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
from scipy.stats import kendalltau
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib


def block_pc1(M):
    """M: [pools x snps] AF (may contain NaN). Returns (pc1_scores[pools], var_explained)."""
    M = np.asarray(M, float)
    col_mean = np.nanmean(M, axis=0)
    col_mean = np.where(np.isfinite(col_mean), col_mean, 0.0)
    idx = np.where(~np.isfinite(M))
    M[idx] = np.take(col_mean, idx[1])                 # impute NaN with per-SNP mean
    M = M - M.mean(axis=0, keepdims=True)              # center
    sd = M.std(axis=0, ddof=0)
    keep = sd > 0
    if keep.sum() == 0:
        return None, 0.0
    M = M[:, keep] / sd[keep]                          # scale -> correlation PCA
    # PCA via SVD; scores = U*S (pool coordinates on PCs)
    U, S, _ = np.linalg.svd(M, full_matrices=False)
    pc1 = U[:, 0] * S[0]
    var_exp = float(S[0] ** 2 / np.sum(S ** 2))
    return pc1, var_exp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True, choices=["snp", "sv", "smallindel"])
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--cmdir", default=f"{lib.GEA}/phase1_replication/class_matrices")  # input
    ap.add_argument("--out", default=f"{lib.GEA}/wza_investigation/pc1")                 # output (relocated)
    ap.add_argument("--min-snps", type=int, default=2, help="blocks with >= this many SNPs")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    recs = pd.read_csv(f"{args.cmdir}/{args.cls}_gen{args.gen}.records.csv")
    pools = pd.read_csv(f"{args.cmdir}/gen{args.gen}.pools.csv")
    af = np.load(f"{args.cmdir}/{args.cls}_gen{args.gen}_af.npy")        # [pools x rec]
    clim = pools[args.climate].to_numpy(float)
    recs["block"] = recs["block"].astype(str).str.replace(r"\.0$", "", regex=True)
    recs["_row"] = np.arange(len(recs))
    print(f"{args.cls} gen{args.gen}: {af.shape[0]} pools x {af.shape[1]:,} records, "
          f"{recs.block.nunique():,} blocks", flush=True)

    rows = []
    for blk, g in recs.groupby("block"):
        if blk in ("", "None") or len(g) < args.min_snps:
            continue
        cols = g["_row"].to_numpy()
        pc1, ve = block_pc1(af[:, cols])
        if pc1 is None:
            continue
        ok = np.isfinite(clim) & np.isfinite(pc1)
        if ok.sum() < 3 or np.ptp(pc1[ok]) == 0:
            continue
        tau, pv = kendalltau(pc1[ok], clim[ok])
        rows.append(dict(block=blk, chrom=g.chrom.iloc[0], pos=float(g.pos.mean()),
                         n_snps=len(g), pc1_var_explained=round(ve, 4), tau=tau, pval=pv))
    out = pd.DataFrame(rows)
    # BH q
    p = out["pval"].to_numpy(); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = (p[o] * n) / (np.arange(n) + 1); q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    out["q"] = np.clip(q, 0, 1)
    f = f"{args.out}/pc1_{args.cls}_gen{args.gen}_{args.climate}.csv"
    out.sort_values("pval").to_csv(f, index=False)
    nbh = int((out.q < 0.05).sum()); nbonf = int((p < 0.05/n).sum())
    print(f"  {n:,} blocks tested | median PC1 var-explained {out.pc1_var_explained.median():.2f} "
          f"| BH q<0.05: {nbh} | Bonferroni: {nbonf}\n  -> {f}", flush=True)


if __name__ == "__main__":
    main()

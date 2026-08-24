#!/usr/bin/env python
"""Pairwise founder LD (r^2) among variants in a window, for the LD-triangle panel.
r^2 = squared Pearson correlation of ALT dosage (0/1) across the 231 founders
(missing mean-imputed). Emits a long table: i,j,pos_i,pos_j,r2 (i<j) + the
diagonal (r2=1) so the apex row is complete.

Usage: compute_ld.py OUT_TAG [LO HI]   (default = 3' region 11533810-11534333)
"""
import sys, numpy as np, pandas as pd, pysam

HERE = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/cam5_replication/pang_cam5"
TAG = sys.argv[1] if len(sys.argv) > 1 else "region"
LO = int(sys.argv[2]) if len(sys.argv) > 3 else 11533810
HI = int(sys.argv[3]) if len(sys.argv) > 3 else 11534333

vcf = pysam.VariantFile(f"{HERE}/cam5_region.vcf.gz")
samp = list(vcf.header.samples)
pos, G = [], []
for r in vcf.fetch("Chr2", LO, HI):
    pos.append(r.pos)
    G.append([(r.samples[s].get("GT", (None,))[0]) for s in samp])
pos = np.array(pos)
G = np.array([[np.nan if g is None else float(g) for g in row] for row in G])  # variants x founders
# mean-impute missing per variant
col_means = np.nanmean(G, axis=1, keepdims=True)
G = np.where(np.isnan(G), col_means, G)
# drop monomorphic (no variance -> undefined corr)
keep = G.std(axis=1) > 0
G, pos = G[keep], pos[keep]
R = np.corrcoef(G)            # variants x variants Pearson r
R2 = R ** 2
n = len(pos)

rows = []
for i in range(n):
    for j in range(i, n):
        rows.append((i, j, int(pos[i]), int(pos[j]), float(R2[i, j])))
df = pd.DataFrame(rows, columns=["i", "j", "pos_i", "pos_j", "r2"])
df.to_csv(f"{HERE}/gg_ld_{TAG}.tsv", sep="\t", index=False)
print(f"[{TAG}] window {LO}-{HI}: {n} polymorphic variants, {len(df)} pairs")
print(f"  r2 quantiles: 50%={np.nanmedian(R2[np.triu_indices(n,1)]):.3f} "
      f"max(offdiag)={np.nanmax(R2[np.triu_indices(n,1)]):.3f}")
# quick check: LD among the 3 Bonferroni-significant sites
sig = [11533967, 11534031, 11534063]
idx = [int(np.where(pos == p)[0][0]) for p in sig if p in pos]
if len(idx) > 1:
    print("  r2 among significant sites:",
          [round(R2[a, b], 2) for a in idx for b in idx if a < b])

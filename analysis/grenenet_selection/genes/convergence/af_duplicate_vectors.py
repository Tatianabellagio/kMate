#!/usr/bin/env python
"""How many non-SNP records share a byte-identical pool-AF vector with a DISTANT record?

RLM3 (Chr4:9.56 Mb) and AT4G02540 (Chr4:1.11 Mb) turned out to have identical gen9
pool-AF values in all 352 pools and identical p0, despite different founder carrier sets in
the panel VCF. A third record at Chr4:11.87 Mb shares the same vector. Identical vectors at
the SAME locus are expected (split records of one event, perfect local LD); identical
vectors 10 Mb apart with different carriers mean the pool frequencies of those records are
not independent estimates -- plausibly allele k-mers shared across repeated / TE-derived
sequence -- and any climate signal on them is not locus-specific.

Hashes every record's AF vector (MAF > 0.05 only, i.e. the GEA-testable set) and reports
groups spanning > 50 kb or several chromosomes. Writes results/af_duplicate_groups.csv.
env: kmate. Compute node, a few minutes.
"""
import os, sys, hashlib
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                       # noqa: E402

M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
ch, pos = idx["chrom"].astype(str), idx["pos"]
# the GEA filters on the gen9 pool-mean MAF, not the seed-mix p0: RLM3 and AT4G02540 sit at
# p0 0.033 and were missed by a p0-only filter. Take either definition.
mu = np.empty(M.shape[1], np.float32)
for s0 in range(0, M.shape[1], 200_000):
    mu[s0:s0 + 200_000] = np.nanmean(np.asarray(M[:, s0:s0 + 200_000], np.float32), 0)
test = np.where((np.minimum(p0, 1 - p0) > 0.05) | (np.minimum(mu, 1 - mu) > 0.05))[0]
h = np.empty(len(test), dtype=object)
for s in range(0, len(test), 100_000):
    cols = test[s:s + 100_000]
    X = np.ascontiguousarray(np.asarray(M[:, cols], np.float32).T)
    X = np.round(np.nan_to_num(X, nan=-1.0), 6)
    h[s:s + len(cols)] = [hashlib.md5(r.tobytes()).hexdigest() for r in X]
D = pd.DataFrame({"row": test, "chrom": ch[test], "pos": pos[test], "h": h,
                  "ref_len": idx["ref_len"][test], "alt_len": idx["alt_len"][test]})
g = D.groupby("h")
span = g.pos.agg(lambda x: x.max() - x.min())
nchr = g.chrom.nunique()
bad = span[(span > 50_000) | (nchr > 1)].index
B = D[D.h.isin(bad)].sort_values(["h", "chrom", "pos"])
B.to_csv(f"{HERE}/results/af_duplicate_groups.csv", index=False)
print(f"testable non-SNP records: {len(D):,}")
print(f"records in a vector group spanning >50 kb or >1 chromosome: {len(B):,} "
      f"({100*len(B)/len(D):.3f}%) in {B.h.nunique():,} groups")
sv = (B.alt_len - B.ref_len).abs() > 50
print(f"  of which SV (>50 bp): {int(sv.sum()):,}   -- SV share among all testable: "
      f"{100*((D.alt_len-D.ref_len).abs()>50).mean():.1f}%")
print("largest groups:", B.groupby("h").size().sort_values(ascending=False).head(5).tolist())

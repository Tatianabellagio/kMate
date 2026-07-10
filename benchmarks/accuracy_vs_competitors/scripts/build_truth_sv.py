#!/usr/bin/env python3
"""Truth AF for SV records of a g0 pool, projected from var_pa founder genotypes.

For g0 pools (each individual = one whole founder) the pool AF of ANY variant is
exactly (founder frequencies) . (founder genotype). We compute that for the SV
records (|alt_len - ref_len| >= SVLEN) directly from the var_pa founder x variant
matrix + pool_weights. Emits both conventions (see build_truth.py):
  truth_af_phys = Sum_f w_f*varpa[f]      / Sum_f w_f          (missing -> REF; physical)
  truth_af      = Sum_f w_f*varpa[f]      / Sum_f w_f*called[f] (MAR; called only)
"""
import argparse, sys
import numpy as np, pandas as pd
import scipy.sparse as sp

ap = argparse.ArgumentParser()
ap.add_argument("--var-pa", required=True)
ap.add_argument("--var-called", required=True)
ap.add_argument("--meta", required=True)
ap.add_argument("--weights", required=True)
ap.add_argument("--svlen", type=int, default=50, help="min |alt_len-ref_len| to call SV")
ap.add_argument("--out", required=True)
a = ap.parse_args()

m = np.load(a.meta, allow_pickle=True)
founders = np.array([str(x) for x in m["founders"]])
chrom = np.array([str(x) for x in m["chrom"]]); pos = m["pos"]
rl = m["ref_len"].astype(int); al = m["alt_len"].astype(int)
vp = sp.load_npz(a.var_pa).tocsr().astype(np.float64)        # F x N
vc = sp.load_npz(a.var_called).tocsr().astype(np.float64)    # F x N

pw = pd.read_csv(a.weights, sep="\t", dtype={"founder": str})
wmap = dict(zip(pw["founder"], pw["weight"]))
w = np.array([wmap.get(f, 0.0) for f in founders], dtype=np.float64)
print(f"[sv-truth] {len(founders)} founders, {int((w>0).sum())} weighted, Sw={w.sum():.4f}", file=sys.stderr)

sv = (np.abs(al - rl) >= a.svlen)
print(f"[sv-truth] {sv.sum():,} SV records (|indel|>={a.svlen}bp) of {len(rl):,}", file=sys.stderr)

num = w @ vp                 # 1 x N  (weighted alt-carrier mass)
wcl = w @ vc                 # 1 x N  (weighted called mass)
num = np.asarray(num).ravel(); wcl = np.asarray(wcl).ravel()
ncall = np.asarray(vc.sum(axis=0)).ravel().astype(int)
wall = w.sum()

df = pd.DataFrame({
    "chrom": chrom[sv], "pos": pos[sv], "ref_len": rl[sv], "alt_len": al[sv],
    "truth_af_phys": num[sv] / wall,
    "truth_af": np.where(wcl[sv] > 0, num[sv] / np.where(wcl[sv] > 0, wcl[sv], 1), np.nan),
    "n_called": ncall[sv],
})
df.to_csv(a.out, sep="\t", index=False)
print(f"[sv-truth] wrote {len(df):,} SV rows -> {a.out}", file=sys.stderr)

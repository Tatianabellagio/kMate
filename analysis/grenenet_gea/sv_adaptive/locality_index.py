#!/usr/bin/env python3
"""Locality index: how much of a haploblock's "selected somewhere" signal (JOINT, S=30
sites) is a linear bio1 gradient, vs a non-monotonic / climate-type-specific / idiosyncratic
residual?

Exact orthogonal chi-square decomposition already implicit in founder_gwas_multisite.py:
    chi2_joint (S df) = z_global^2 (1 df, generalist) + z_clim^2 (1 df, linear bio1 gradient)
                        + RESIDUAL (S-2 df, everything NOT a generalist shift or straight-
                          line climate gradient -- local/non-monotonic response)
locality = RESIDUAL / chi2_joint, in [0,1] by construction (undefined only at chi2_joint==0).
This is the direct block-level test of "is SV selection more LOCAL (higher residual share)
than SNP selection, even though neither tracks the bio1 gradient (CLIMATE is null for both)?"

Reuses lib.multisite_gwas_raw (Z/C from multisite_founder_gwas_clq90_pc1.npz) -- no new GWAS
run. Tests has_sv==1 blocks against a size-matched null of same-n_kept-bin controls, both
genome-wide and restricted to top-JOINT ("actually selected somewhere") subsets, since the
ratio is just noise for blocks with ~0 signal.
Run in `basic` env from the kmate repo root. Deterministic (seed=0).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path
import numpy as np, pandas as pd
import lib

OUT = Path("results/grenenet_gea/sv_adaptive")
NPERM = 10000
EDGES = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 250, 10**9]
TOP_JOINT = [0.01, 0.05, 0.10]     # "actually selected somewhere" subsets, on top of all_blocks

raw_all = lib.multisite_gwas_raw("clq90_pc1")
MIN_MAF = 0.01
raw = lib.maf_filter(raw_all, MIN_MAF)
print(f"MAF filter (founder MAF>={MIN_MAF:.0%}): kept {raw['M']:,}/{raw_all['M']:,} markers")
chi2_joint, z_global, z_clim = raw["chi2_joint"], raw["z_global"], raw["z_clim"]
residual = np.clip(chi2_joint - z_global ** 2 - z_clim ** 2, 0, None)
locality = np.divide(residual, chi2_joint, out=np.full_like(chi2_joint, np.nan), where=chi2_joint > 0)
assert np.nanmax(locality) <= 1 + 1e-6, "locality index exceeds 1 -- decomposition is not orthogonal"
print(f"{raw['M']:,} markers | residual df = {raw['S'] - 2} | "
      f"median locality (all markers) = {np.nanmedian(locality):.3f}")

mk = pd.DataFrame({"unit": raw["unit"], "p_joint": raw["p_joint"],
                   "chi2_joint": chi2_joint, "locality": locality})
# lead marker per block = the block's own most-significant-anywhere marker (self-consistent
# triple: don't mix chi2_joint from one marker with locality computed from another)
lead_idx = mk.groupby("unit")["p_joint"].idxmin().to_numpy()
lead = mk.loc[lead_idx].reset_index(drop=True)

L = pd.read_csv(OUT / "sv_landscape_clq0.9.csv")
df = lead.merge(L, left_on="unit", right_on="block_id", how="inner").reset_index(drop=True)
df["bin"] = pd.cut(df.n_kept, EDGES, right=False, labels=False)
df = df[df.bin.notna()].reset_index(drop=True)
df["bin"] = df["bin"].astype(int)
df.to_csv(OUT / "locality_index.csv", index=False)
print(f"{len(df):,} blocks | median locality (block-lead) = {df.locality.median():.3f} | "
      f"{int(df.has_sv.sum())} carry an SV")

# ---- is SV-bearing more LOCAL than size-matched controls? genome-wide + top-JOINT subsets ----
rows = []
subsets = [("all_blocks", df)] + [
    (f"top{int(f * 100)}pct_joint", df.nsmallest(int(round(f * len(df))), "p_joint"))
    for f in TOP_JOINT]
for tag, pool in subsets:
    pool = pool.reset_index(drop=True)
    bins, vals = pool["bin"].to_numpy(), pool["locality"].to_numpy()
    idx = np.where(pool.has_sv.to_numpy() == 1)[0]
    if len(idx) < 5:
        print(f"  {tag}: too few SV blocks ({len(idx)}) to test, skipping")
        continue
    obs, null, ratio, p = lib.matched_perm_test(vals, bins, idx, NPERM)
    rows.append(dict(subset=tag, n_sv_blocks=len(idx), n_pool=len(pool),
                     obs_locality=round(obs, 4), null_locality=round(null, 4),
                     ratio=round(ratio, 3), p_perm=round(p, 4)))
res = pd.DataFrame(rows)
res.to_csv(OUT / "locality_enrichment.csv", index=False)
print("\n== is SV-bearing more LOCAL (higher residual share) than size-matched controls? ==")
print(res.to_string(index=False))
print(f"\n[done] -> {OUT}/locality_index.csv, locality_enrichment.csv")

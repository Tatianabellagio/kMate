#!/usr/bin/env python
"""Assemble the [pools x blocks] PC1 matrix for the structure-corrected block test.

Each block -> its dominant haplotype axis (PC1 across pools; same block_pc1 as
run_block_pc1.py). Stacking blocks gives a [pools x n_blocks] matrix whose columns
are density-invariant block summaries. Feeding THAT to LFMM (K=16) then gives a
per-block test that is BOTH density-invariant (one column/block) AND structure-
corrected (latent factors) -> should pull genomic inflation lambda down to ~1.

Writes (for class CLS, gen G), reusing the run_lfmm_lastgen.R binary contract:
  pc1lfmm_{CLS}_gen{G}_Y.f64     float64 C-order [n_pools x n_blocks]
  pc1lfmm_{CLS}_gen{G}_dims.txt   "n_pools n_blocks"
  pc1lfmm_{CLS}_gen{G}_env.csv    z-scored bio (pool order)
  pc1lfmm_{CLS}_gen{G}_blocks.csv block,chrom,pos,n_snps,pc1_var_explained (col order)
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_block_pc1 import block_pc1
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True, choices=["snp", "sv", "smallindel"])
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--cmdir", default=f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices")  # input
    ap.add_argument("--out", default=f"{lib.GEA}/wza/investigation/results/pc1")                 # output (relocated)
    ap.add_argument("--min-snps", type=int, default=2)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    recs = pd.read_csv(f"{args.cmdir}/{args.cls}_gen{args.gen}.records.csv")
    pools = pd.read_csv(f"{args.cmdir}/gen{args.gen}.pools.csv")
    af = np.load(f"{args.cmdir}/{args.cls}_gen{args.gen}_af.npy")          # [pools x rec]
    recs["block"] = recs["block"].astype(str).str.replace(r"\.0$", "", regex=True)
    recs["_row"] = np.arange(len(recs))
    n_pools = af.shape[0]

    cols, meta = [], []
    for blk, g in recs.groupby("block"):
        if blk in ("", "None") or len(g) < args.min_snps:
            continue
        pc1, ve = block_pc1(af[:, g["_row"].to_numpy()])
        if pc1 is None or not np.all(np.isfinite(pc1)):
            continue
        cols.append(pc1.astype(np.float64))
        meta.append((blk, g.chrom.iloc[0], float(g.pos.mean()), len(g), round(ve, 4)))
    Y = np.column_stack(cols)                                              # [pools x blocks]
    print(f"{args.cls} gen{args.gen}: PC1 matrix [{Y.shape[0]} x {Y.shape[1]:,} blocks]", flush=True)

    stem = f"{args.out}/pc1lfmm_{args.cls}_gen{args.gen}"
    Y.astype(np.float64).tofile(f"{stem}_Y.f64")
    open(f"{stem}_dims.txt", "w").write(f"{Y.shape[0]} {Y.shape[1]}\n")
    x = pools[args.climate].to_numpy(float); z = (x - x.mean()) / x.std(ddof=0)
    pd.DataFrame({args.climate: z}).to_csv(f"{stem}_env.csv", index=False)
    pd.DataFrame(meta, columns=["block", "chrom", "pos", "n_snps", "pc1_var_explained"]).to_csv(
        f"{stem}_blocks.csv", index=False)
    print(f"  -> {stem}_Y.f64 / _dims.txt / _env.csv / _blocks.csv", flush=True)


if __name__ == "__main__":
    main()

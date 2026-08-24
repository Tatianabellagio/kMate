#!/usr/bin/env python
"""Build LFMM inputs for the phase-1 replication (per class, per gen).

Mirrors phase-1's LFMM exactly: Y = Δp (pool AF - founding p0) per record, env =
standardized bioclim, lfmm_ridge K=16, lfmm_test calibrate="gif" (the R step).
Phase-1's hapFIRE Δp had no missing values; kMate Δp does (records finite in
>=50% of pools by the class filter), so we impute each record's NaN with that
record's across-pool mean Δp (standard LFMM locus-mean imputation).

Writes (for class CLS, gen G):
  lfmm_{CLS}_gen{G}_Y.f64        float64 binary, C-order [n_pools x n_records]
  lfmm_{CLS}_gen{G}_dims.txt      "n_pools n_records"
  lfmm_{CLS}_gen{G}_env.csv       one col, z-scored bio (pool order matches Y rows)
  (records come from the class matrix; the WZA join reads {CLS}_gen{G}.records.csv)
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

_P0 = {"snp": "p0_snp.npy", "sv": "p0_nonsnp.npy", "smallindel": "p0_nonsnp.npy",
       "nonsnp": "p0_nonsnp.npy"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True,
                    choices=["snp", "sv", "smallindel", "nonsnp"])
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--cmdir", default=f"{lib.GEA}/phase1_replication/results/class_matrices")
    ap.add_argument("--out", default=f"{lib.GEA}/phase1_replication/results/lfmm")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    recs = pd.read_csv(f"{args.cmdir}/{args.cls}_gen{args.gen}.records.csv")
    pools = pd.read_csv(f"{args.cmdir}/gen{args.gen}.pools.csv")
    af = np.load(f"{args.cmdir}/{args.cls}_gen{args.gen}_af.npy")          # [pools x rec], NaN=missing
    p0 = np.load(f"{lib.AF_STORE}/{_P0[args.cls]}")[recs['col'].to_numpy()]  # founding AF per kept record
    n_pools, n_rec = af.shape
    assert n_rec == len(recs) and n_pools == len(pools), (af.shape, len(recs), len(pools))

    dp = af - p0[None, :]                                  # Δp [pools x rec]
    # impute NaN with each record's across-pool mean Δp (column mean)
    with np.errstate(invalid="ignore"):
        colmean = np.nanmean(dp, axis=0)
    colmean = np.where(np.isfinite(colmean), colmean, 0.0)
    miss = ~np.isfinite(dp)
    dp[miss] = np.take(colmean, np.where(miss)[1])
    print(f"{args.cls} gen{args.gen}: Δp [{n_pools} x {n_rec:,}] | imputed "
          f"{miss.sum():,} cells ({100*miss.mean():.1f}%)", flush=True)

    stem = f"{args.out}/lfmm_{args.cls}_gen{args.gen}"
    dp.astype(np.float64).tofile(f"{stem}_Y.f64")          # C-order [pools x rec]
    open(f"{stem}_dims.txt", "w").write(f"{n_pools} {n_rec}\n")
    x = pools[args.climate].to_numpy(float)
    z = (x - x.mean()) / x.std(ddof=0)                     # standardize env (LFMM convention)
    pd.DataFrame({args.climate: z}).to_csv(f"{stem}_env.csv", index=False)
    print(f"  -> {stem}_Y.f64 / _dims.txt / _env.csv (env z-scored {args.climate})", flush=True)


if __name__ == "__main__":
    main()

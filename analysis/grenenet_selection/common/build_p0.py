#!/usr/bin/env python
"""Build founding (gen-0) p0 from the 8 SEEDMIX reps, aligned to the AF store.

p0[record] = NaN-aware mean of alt_freq over the 8 SEEDMIX kMate reps (the
founding seed mix; phase-1's average_seedmix_p0), split SNP / non-SNP to match
the store's index.

Outputs (in the store dir):
  p0_nonsnp.npy  float32 [n_nonsnp]   founding AF for the SV-GEA substrate
  p0_snp.npy     float32 [n_snp]

These are the Delta-p / log(p1/p0) reference for every pool matrix.
"""
from __future__ import annotations
import glob, json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib


def main():
    store = lib.AF_STORE
    snp = np.load(f"{store}/snp_mask.npy")
    meta = json.load(open(f"{store}/meta.json"))
    n_full = meta["n_full"]
    fs = sorted(f for f in glob.glob(f"{lib.SEEDMIX}/SEEDMIX_S*.tsv")
                if "_Chr" not in os.path.basename(f))
    if not fs:
        raise SystemExit(f"no SEEDMIX TSVs in {lib.SEEDMIX}")
    ssum = np.zeros(n_full); scnt = np.zeros(n_full)
    for f in fs:
        a = pd.read_csv(f, sep="\t", usecols=["alt_freq"]).alt_freq.to_numpy(float)
        if len(a) != n_full:
            raise SystemExit(f"{f}: {len(a):,} rows != panel {n_full:,}")
        ok = np.isfinite(a)
        ssum[ok] += a[ok]; scnt += ok
        print(f"  {os.path.basename(f)}: merged", flush=True)
    p0 = (ssum / np.where(scnt > 0, scnt, np.nan)).astype(np.float32)
    np.save(f"{store}/p0_nonsnp.npy", p0[~snp])
    np.save(f"{store}/p0_snp.npy", p0[snp])
    print(f"p0 from {len(fs)} SEEDMIX reps -> p0_nonsnp.npy ({(~snp).sum():,}) "
          f"+ p0_snp.npy ({snp.sum():,}); nonsnp mean={np.nanmean(p0[~snp]):.4f} "
          f"NaN={np.isnan(p0[~snp]).mean()*100:.3f}%")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""LFMM input for the clq0.9 HAPLOBLOCK last-gen matrix (structure-corrected GEA).

Uses the SAME response as the Kendall run -- the raw gen9 haploblock frequencies
(hap_gen9_af.npy) -- so LFMM's K latent factors isolate exactly the founder/clade
structure that inflated the raw Kendall scan (lambda=8.95). No Delta-p / p0 subtraction
(that is the variant-pipeline convention); no imputation (hap freq has no NaN).

Writes (consumed by run_lfmm_lastgen.R):
  lfmm_hap_gen9_Y.f64    float64 C-order [n_pools x n_hap]
  lfmm_hap_gen9_dims.txt "n_pools n_hap"
  lfmm_hap_gen9_env.csv  z-scored bio (pool order matches Y rows)
Env: kmate.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--cmdir", default=f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices")
    ap.add_argument("--out", default=f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/lfmm")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    af = np.load(f"{a.cmdir}/hap_gen9_af.npy")                 # [pools x hap], no NaN
    pools = pd.read_csv(f"{a.cmdir}/gen9.pools.csv")
    recs = pd.read_csv(f"{a.cmdir}/hap_gen9.records.csv")
    np_, nh = af.shape
    assert np_ == len(pools) and nh == len(recs), (af.shape, len(pools), len(recs))
    assert np.isfinite(af).all(), "hap freq matrix has non-finite entries"
    print(f"hap gen9: Y [{np_} pools x {nh:,} haploblocks] (raw freq, matches Kendall)", flush=True)

    stem = f"{a.out}/lfmm_hap_gen9"
    af.astype(np.float64).tofile(f"{stem}_Y.f64")             # C-order [pools x hap]
    open(f"{stem}_dims.txt", "w").write(f"{np_} {nh}\n")
    x = pools[a.climate].to_numpy(float)
    z = (x - x.mean()) / x.std(ddof=0)
    pd.DataFrame({a.climate: z}).to_csv(f"{stem}_env.csv", index=False)
    print(f"  -> {stem}_Y.f64 / _dims.txt / _env.csv (env z-scored {a.climate})", flush=True)


if __name__ == "__main__":
    main()

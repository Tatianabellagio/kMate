#!/usr/bin/env python
"""Cache per-sample genome-wide founder h (231-vector) for every sequenced timepoint sample.

One-time Lustre-heavy pass (genome_h loads 5 window npz per sample). Reused by the
Price/fitness analyses (plot composition <-> realized census fitness). Output:
analysis/grenenet_gea/fitness/sample_genome_h.npz  ->  samples (str[N]), H (N x 231), founders.
Env: kmate. Run via run_sample_h_cache.sbatch.
"""
import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from lib import genome_h, CHROMS

WIN = lib.OUT
OUT = "analysis/grenenet_gea/fitness"


def main():
    os.makedirs(OUT, exist_ok=True)
    founders = np.load(glob.glob(f"{WIN}/*_Chr1.h_per_chrom.npz")[0],
                       allow_pickle=True)["founders"].astype(str)
    pt = lib.pool_table()
    samples = sorted(pt.sampleid.astype(str).unique())
    print(f"{len(samples)} unique timepoint samples | {len(founders)} founders", flush=True)

    keep, H = [], []
    miss = 0
    for i, s in enumerate(samples):
        h = genome_h(s, WIN)
        if h is None:
            miss += 1
            continue
        keep.append(s); H.append(h)
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(samples)} done ({miss} missing)", flush=True)
    H = np.vstack(H)
    np.savez(f"{OUT}/sample_genome_h.npz",
             samples=np.array(keep), H=H, founders=founders)
    print(f"[done] cached {len(keep)} samples ({miss} missing) -> {OUT}/sample_genome_h.npz "
          f"shape {H.shape}", flush=True)


if __name__ == "__main__":
    main()

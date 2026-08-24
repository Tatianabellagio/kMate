#!/usr/bin/env python
"""Build per-GENERATION, flower-weighted site_gen_plot POOL matrices.

The 2,168 cohort samples are individual collection TIMEPOINTS. The GrENE-Net
analysis unit is the `site_gen_plot` POOL: within a (site, generation, plot) the
timepoints are merged into one allele frequency, weighted by flowers collected:

    p_pool[record] = sum_t flowers_t * p_t[record] / sum_t flowers_t

(verified vs the phase-1 merged_hapFIRE columns to ~1e-15). NaN timepoints are
dropped per-record and the weights renormalised over the finite ones. The 745
pools reproduce the phase-1 merge units exactly; plots are kept SEPARATE (reps)
so the mixed models keep site/ecotype random structure.

Per generation `g` and `kind` in {nonsnp, snp}:
  pool_gen{g}_{kind}_af.npy   float32 [n_pools_g x n_records]  (NaN where no data)
  pool_gen{g}_{kind}.meta.csv pool, site, plot, generation, n_timepoints,
                              total_flowers, mean_coverage, bio1..bio19
                              (row-aligned to the matrix)

Records share the panel order; the index is the store's index_{kind}.npz.
Pooled means are floats -> stored float32 (gen1 nonsnp ~ 326 x 2.25M ~ 2.9 GB).

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY analysis/grenenet_selection/common/build_pool_matrix.py --kind nonsnp --gens 1 2 3 \
      --out analysis/grenenet_selection/pool_matrices
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib


def build(gen: int, kind: str, store: str, out: str):
    pt = lib.pool_table(store)
    pt = pt[pt.generation == gen]
    clim = lib.load_climate()
    pools = sorted(pt.pool.unique(),
                   key=lambda p: tuple(int(x) for x in p.split("_")))
    n_rec = json.load(open(f"{store}/meta.json"))[f"n_{kind}"]
    os.makedirs(out, exist_ok=True)
    mat = np.lib.format.open_memmap(
        f"{out}/pool_gen{gen}_{kind}_af.npy", mode="w+",
        dtype=np.float32, shape=(len(pools), n_rec))
    rows = []
    for i, pool in enumerate(pools):
        members = pt[pt.pool == pool]
        acc = np.zeros(n_rec, dtype=np.float64)      # flower-weighted AF sum
        wsum = np.zeros(n_rec, dtype=np.float64)     # finite-weight sum per record
        for s, w in zip(members.sampleid, members.flowerscollected.astype(float)):
            f = lib.decode_af(np.load(f"{store}/af_{kind}/{s}.npy"))
            ok = np.isfinite(f)
            acc[ok] += w * f[ok]
            wsum[ok] += w
        mat[i] = np.where(wsum > 0, acc / np.where(wsum > 0, wsum, 1), np.nan)
        site = int(members.site.iloc[0])
        rows.append(dict(pool=pool, site=site, plot=int(members["plot"].iloc[0]),
                         generation=gen, n_timepoints=len(members),
                         total_flowers=float(members.flowerscollected.sum()),
                         mean_coverage=float(members.coverage.mean())))
        if (i + 1) % 50 == 0:
            print(f"  gen{gen} {kind}: {i+1}/{len(pools)} pools", flush=True)
    mat.flush()
    meta = pd.DataFrame(rows).join(clim, on="site")
    meta.to_csv(f"{out}/pool_gen{gen}_{kind}.meta.csv", index=False)
    print(f"gen{gen} {kind}: {len(pools)} pools x {n_rec:,} records "
          f"({mat.nbytes/1e9:.1f} GB) -> {out}/pool_gen{gen}_{kind}_af.npy")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default=lib.AF_STORE)
    ap.add_argument("--out", default=f"{lib.GEA}/pool_matrices")
    ap.add_argument("--kind", default="nonsnp", choices=["nonsnp", "snp"])
    ap.add_argument("--gens", type=int, nargs="+", default=[1, 2, 3])
    args = ap.parse_args()
    for g in args.gens:
        build(g, args.kind, args.store, args.out)


if __name__ == "__main__":
    main()

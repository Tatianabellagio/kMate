#!/usr/bin/env python
"""Phase-1 Kendall-tau GEA on a filtered class matrix (one class, one gen, bio1).

Mirrors the phase-1 `kendall_tau/.../kendall.py`: for each kept record, Kendall's
tau-b between the record's per-pool allele frequency and the pool's climate value
(bio1) across the generation's site_gen_plot pools. Pools within a site share a
climate value (ties; tau-b handles them) and are pseudoreplicated — exactly as in
phase-1; structure correction is the job of LFMM / WZA, not this step.

Reads the outputs of build_class_matrices.py:
  {class}_gen{g}_af.npy        [pools x records]
  {class}_gen{g}.records.csv   record meta (incl. maf, block)
  gen{g}.pools.csv             pool rowmeta with bio1..bio19

Output (--out):
  kendall_{class}_gen{g}_{climate}.csv   per-record: chrom,pos,ref_len,alt_len,
                                         maf,block,tau,pval  (-> WZA input)
Columns named so wza_script.py consumes it directly: needs `MAF`, a p-value
summary stat, and a window column (`block`).

Usage:
  PY=.../kmate/bin/python
  $PY run_kendall.py --class sv --gen 1 --climate bio1
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy.stats import kendalltau
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

_X = None      # climate vector per pool
_AFT = None    # AF transposed [records x pools]


def _chunk(rng):
    a, b = rng
    tau = np.full(b - a, np.nan); pv = np.full(b - a, np.nan)
    for i in range(a, b):
        y = _AFT[i]
        ok = np.isfinite(y)
        if ok.sum() >= 3 and np.ptp(y[ok]) > 0:        # need variation + ≥3 finite pools
            t, p = kendalltau(_X[ok], y[ok])
            tau[i - a] = t; pv[i - a] = p
    return a, tau, pv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True,
                    choices=["snp", "sv", "smallindel"])
    ap.add_argument("--gen", type=int, required=True)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--cmdir", default=f"{lib.GEA}/phase1_replication/class_matrices")
    ap.add_argument("--out", default=f"{lib.GEA}/phase1_replication/kendall")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    recs = pd.read_csv(f"{args.cmdir}/{args.cls}_gen{args.gen}.records.csv")
    pools = pd.read_csv(f"{args.cmdir}/gen{args.gen}.pools.csv")
    x = pools[args.climate].to_numpy(float)
    af = np.load(f"{args.cmdir}/{args.cls}_gen{args.gen}_af.npy")   # [pools x rec]
    assert af.shape[0] == len(x), (af.shape, len(x))
    assert af.shape[1] == len(recs), (af.shape, len(recs))
    n_rec = af.shape[1]
    print(f"{args.cls} gen{args.gen}: {af.shape[0]} pools x {n_rec:,} records vs "
          f"{args.climate} (range {np.nanmin(x):.1f}-{np.nanmax(x):.1f})", flush=True)

    global _X, _AFT
    _X = x
    _AFT = np.ascontiguousarray(af.T)      # [records x pools]
    del af

    step = 20000
    ranges = [(a, min(a + step, n_rec)) for a in range(0, n_rec, step)]
    tau = np.empty(n_rec, np.float64); pv = np.empty(n_rec, np.float64)
    with Pool(args.threads) as pool:
        for a, t, p in pool.imap_unordered(_chunk, ranges):
            tau[a:a + len(t)] = t; pv[a:a + len(p)] = p

    recs = recs.assign(tau=tau, pval=pv)
    recs = recs.rename(columns={"maf": "MAF"})
    out = f"{args.out}/kendall_{args.cls}_gen{args.gen}_{args.climate}.csv"
    recs[["chrom", "pos", "ref_len", "alt_len", "MAF", "block",
          "tau", "pval"]].to_csv(out, index=False)
    fin = np.isfinite(pv)
    print(f"  tested {int(fin.sum()):,}/{n_rec:,} | p<0.05: {int((pv[fin]<0.05).sum()):,} "
          f"| min p={np.nanmin(pv):.2e}\n  -> {out}", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""LFMM input files for one class x climate axis, for both arms of the GEA (steps 2 and 4 of _submit_gea.sh).

  Y    delta p = pool frequency - founding frequency p0, [pools x tested units].
       Record arm (snp, smallindel): built here from --cmdir's {cls}_gen{g}_af.npy and af_store's p0.
       Event arm (sv, nonsnp): --dp, the delta p of Truvari events built by event_gea/build_event_matrix.py.
  env  the climate axis from --cmdir's gen{g}.pools.csv, z-scored over pools.

Writes {out}/lfmm_{cls}_gen{g}_Y.f64 (float64, C order), _dims.txt ("n_pools n_units") and _env.csv, which
multiaxis/run_lfmm_nogif.R reads. Rows are the pools of gen{g}.pools.csv; columns of Y are the rows of the
records / units table, in order, so the p-values join back by position.

Environment: kmate.
Usage: build_lfmm_input.py --class snp --climate bio1 --cmdir <class matrices> --out <work dir> [--dp event_dp.npy]
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

P0 = {"snp": "p0_snp.npy", "smallindel": "p0_nonsnp.npy", "sv": "p0_nonsnp.npy", "nonsnp": "p0_nonsnp.npy"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True, choices=list(P0))
    ap.add_argument("--climate", required=True, help="a column of gen{g}.pools.csv: bio1-19, pc1-3")
    ap.add_argument("--cmdir", required=True, help="build_class_matrices.py output (pools; records for the record arm)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--dp", help="event arm: a ready delta p matrix [pools x units] (.npy)")
    a = ap.parse_args()
    lib.check_version_dir(a.cmdir, "--cmdir")               # never mix panel versions (KMATE_PANEL)
    if lib.VSUF and lib.VSUF.lstrip("_") not in os.path.realpath(a.out):
        raise SystemExit(f"--out {a.out} is not under a {lib.PANEL} directory; refusing")

    pools = pd.read_csv(f"{a.cmdir}/gen{a.gen}.pools.csv")
    assert not pools.site.astype(int).isin(lib.EXCLUDED_SITES).any(), "an excluded site reached the GEA input"
    if a.dp:
        dp = np.load(a.dp)
    else:
        col = pd.read_csv(f"{a.cmdir}/{a.cls}_gen{a.gen}.records.csv", usecols=["col"])["col"].to_numpy()
        af = np.load(f"{a.cmdir}/{a.cls}_gen{a.gen}_af.npy")
        assert af.shape[1] == len(col)
        dp = af - np.load(f"{lib.AF_STORE}/{P0[a.cls]}")[col][None, :]
    assert dp.shape[0] == len(pools), (dp.shape, len(pools))
    assert np.isfinite(dp).all(), "delta p has missing values; LFMM needs a complete matrix"

    os.makedirs(a.out, exist_ok=True)
    stem = f"{a.out}/lfmm_{a.cls}_gen{a.gen}"
    dp.astype(np.float64).tofile(f"{stem}_Y.f64")
    with open(f"{stem}_dims.txt", "w") as fh:
        fh.write(f"{dp.shape[0]} {dp.shape[1]}\n")
    x = pools[a.climate].to_numpy(float)
    pd.DataFrame({a.climate: (x - x.mean()) / x.std(ddof=0)}).to_csv(f"{stem}_env.csv", index=False)
    print(f"{a.cls} gen{a.gen} {a.climate}: Y [{dp.shape[0]} x {dp.shape[1]:,}] -> {stem}_*", flush=True)


if __name__ == "__main__":
    main()

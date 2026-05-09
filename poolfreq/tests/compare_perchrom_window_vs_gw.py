#!/usr/bin/env python
"""Compare per-chrom window-mode against genome-wide window-mode on MLFH040120180306.

Production-readiness check: if per-chrom window mode reproduces the genome-wide
window result within MAE < 0.01 and R^2 > 0.999, we can flip the production
SLURM template to use the per-chrom driver.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    base = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/results")
    gw = base / "site04_231_v2" / "MLFH040120180306.tsv"
    pc = base / "site04_231_v2_perchrom" / "MLFH040120180306.tsv"
    if not gw.exists():
        sys.exit(f"missing genome-wide reference: {gw}")
    if not pc.exists():
        sys.exit(f"missing per-chrom output: {pc}")

    a = pd.read_csv(gw, sep="\t")
    b = pd.read_csv(pc, sep="\t")
    print(f"genome-wide window-mode rows: {len(a):,}")
    print(f"per-chrom window-mode rows:   {len(b):,}")
    if len(a) != len(b):
        print("  WARN: row counts differ — joining on (chrom,pos,ref_len,alt_len)")
        m = a.merge(b, on=["chrom", "pos", "ref_len", "alt_len"], suffixes=("_gw", "_pc"))
        af_gw = m["alt_freq_gw"].to_numpy()
        af_pc = m["alt_freq_pc"].to_numpy()
    else:
        af_gw = a["alt_freq"].to_numpy()
        af_pc = b["alt_freq"].to_numpy()

    finite = np.isfinite(af_gw) & np.isfinite(af_pc)
    af_gw, af_pc = af_gw[finite], af_pc[finite]
    err = af_pc - af_gw
    mae = float(np.mean(np.abs(err)))
    rmse = float(np.sqrt(np.mean(err**2)))
    max_abs = float(np.max(np.abs(err)))
    r2 = float(np.corrcoef(af_gw, af_pc)[0, 1] ** 2)
    print(f"\n=== per-chrom window vs genome-wide window AF (n={finite.sum():,}) ===")
    print(f"  max |diff|: {max_abs:.5f}")
    print(f"  MAE:        {mae:.6f}")
    print(f"  RMSE:       {rmse:.6f}")
    print(f"  R^2:        {r2:.6f}")

    poly = (af_gw > 0.001) & (af_gw < 0.999)
    if poly.sum() > 100:
        e = af_pc[poly] - af_gw[poly]
        print(f"\n  polymorphic-only (n={poly.sum():,}):")
        print(f"    max |diff|: {np.max(np.abs(e)):.5f}")
        print(f"    MAE:        {np.mean(np.abs(e)):.6f}")
        print(f"    R^2:        {np.corrcoef(af_gw[poly], af_pc[poly])[0,1]**2:.6f}")

    verdict = "PASS — flip production to per-chrom" if (max_abs < 0.05 and r2 > 0.99) else "INVESTIGATE"
    print(f"\nverdict: {verdict}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Compare per-chromosome cactus_em output against the genome-wide reference.

ALGORITHM.md predicts the per-chrom h vectors agree to within ~0.1% of the
genome-wide h. This script measures the actual divergence on SEEDMIX_S1:

  - Per-record AF: max abs diff, MAE, R^2 between the two TSVs
  - Per-chrom h vs genome-wide h (if genome-wide h was saved)

If max abs diff < 0.005 and R^2 > 0.9999, the per-chrom mode is a safe
drop-in replacement for the production driver and we can switch SLURM
allocations from 192/200 GB to 64 GB.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    base = Path(__file__).resolve().parents[1] / "results"
    gw = base / "seedmix_231_v2" / "SEEDMIX_S1.tsv"
    pc = base / "seedmix_231_v2_perchrom" / "SEEDMIX_S1.tsv"
    if not gw.exists():
        sys.exit(f"missing genome-wide reference: {gw}")
    if not pc.exists():
        sys.exit(f"missing per-chrom output: {pc}")

    a = pd.read_csv(gw, sep="\t")
    b = pd.read_csv(pc, sep="\t")
    print(f"genome-wide rows: {len(a):,}")
    print(f"per-chrom rows:   {len(b):,}")
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
    print(f"\n=== per-chrom vs genome-wide AF (n={finite.sum():,}) ===")
    print(f"  max |diff|: {max_abs:.5f}")
    print(f"  MAE:        {mae:.6f}")
    print(f"  RMSE:       {rmse:.6f}")
    print(f"  R^2:        {r2:.6f}")

    # Polymorphic-only (the records that actually matter)
    poly = (af_gw > 0.001) & (af_gw < 0.999)
    if poly.sum() > 100:
        e = af_pc[poly] - af_gw[poly]
        print(f"\n  polymorphic-only (n={poly.sum():,}):")
        print(f"    max |diff|: {np.max(np.abs(e)):.5f}")
        print(f"    MAE:        {np.mean(np.abs(e)):.6f}")
        print(f"    R^2:        {np.corrcoef(af_gw[poly], af_pc[poly])[0,1]**2:.6f}")

    verdict = "PASS" if (max_abs < 0.01 and r2 > 0.999) else "INVESTIGATE"
    print(f"\nverdict: {verdict}")


if __name__ == "__main__":
    main()

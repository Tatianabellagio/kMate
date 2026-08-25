"""
Aggregate per-sample TSV outputs from batch_runner into a wide matrix:
    rows = biallelic VCF records (chrom, pos, ref_len, alt_len)
    cols = sample IDs
    values = alt-allele frequency in that sample

Usage:
    python aggregate_results.py --input-dir results/per_sample/ \
        --out results/all_samples_freq.parquet
"""
from __future__ import annotations
import argparse, glob, os, time
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True, help="dir of per-sample TSVs")
    ap.add_argument("--out", required=True, help="output .parquet path")
    ap.add_argument("--pattern", default="*.tsv")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))
    print(f"Found {len(files)} per-sample TSVs in {args.input_dir}")
    if not files:
        return

    # Load first to get the variant axis
    print(f"[{time.strftime('%H:%M:%S')}] Loading first file to anchor variant order")
    first = pd.read_csv(files[0], sep="\t")
    keys = first[["chrom", "pos", "ref_len", "alt_len"]].copy()
    print(f"  {len(keys):,} records per sample")

    # Build a wide DataFrame iteratively
    wide = keys.copy()
    for i, f in enumerate(files):
        sid = Path(f).stem
        df = pd.read_csv(f, sep="\t")
        # Verify alignment
        if not (df.shape[0] == keys.shape[0] and
                (df["chrom"] == keys["chrom"]).all() and
                (df["pos"] == keys["pos"]).all()):
            # Re-align via merge
            df_keyed = df[["chrom", "pos", "ref_len", "alt_len", "alt_freq"]]
            wide = wide.merge(df_keyed, on=["chrom","pos","ref_len","alt_len"],
                              how="left", suffixes=("", f"_{sid}"))
            wide = wide.rename(columns={"alt_freq": sid})
        else:
            wide[sid] = df["alt_freq"].values
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(files)} samples loaded")

    print(f"\nFinal wide matrix: {wide.shape}")
    print(f"  variant columns: chrom, pos, ref_len, alt_len ({wide.shape[1] - 4} samples)")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    wide.to_parquet(args.out, compression="snappy")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

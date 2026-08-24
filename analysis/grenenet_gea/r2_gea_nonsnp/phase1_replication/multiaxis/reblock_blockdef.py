#!/usr/bin/env python
"""Reblock the multi-axis per-record model outputs onto an ARBITRARY block definition.

Generalises `reblock_multiaxis.py` (which is hardwired to lib.assign_clq_blocks) so the
coarser/older block definitions can be tested:

  --blockdef hapfire    phase-1 hapFIRE LD haploblocks (lib.assign_ld_blocks). TILES the
                        genome via a nearest-SNP map -> 100% of records keep a block.
  --blockdef clq0.5     BigLD r2>=0.5 interval blocks (coarser than production clq0.9)
  --blockdef clq0.9     production (identical to reblock_multiaxis.py)

Motivation (2026-07-27): A. thaliana is a heavy selfer, so LD extends far and the
production clq0.9 blocks may be far too fine -- measured on the current non-SNP data
they hold a median of 3 records/block and, because clq blocks are intervals with gaps,
they DISCARD 40.7% of records that fall between blocks. hapFIRE keeps 100% (median 8
records/block, 14,043 blocks vs 42,749) -> more evidence per test and a ~3x lighter
multiple-testing burden.

Writes the WZA-input schema (chrom,pos,ref_len,alt_len,MAF,block,pval) into a
blockdef-specific dir so production clq0.9 inputs are never clobbered.

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY reblock_blockdef.py --axis bio1 --cls nonsnp --blockdef hapfire
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import lib

MA = f"{lib.GEA}/phase1_replication/results/multiaxis"

SRC = {
    "kendall":  ("kendall/kendall_{cls}_gen9_{axis}.csv", "pval"),
    "lfmm":     ("lfmm/lfmm_{cls}_gen9_{axis}.csv", "pval"),
    "binomial": ("quasibinom/quasibinom_lf16_{cls}_gen9_{axis}.csv", "pval_quasi"),
}


def assign(blockdef, chrom, pos, how="strict"):
    """how='strict'  -> lib.assign_clq_blocks, interval containment (DROPS gap variants)
       how='tiling'  -> blocks_tiling.assign_tiling, HapFM's gap-free rule (drops nothing)"""
    if blockdef == "hapfire":
        return lib.assign_ld_blocks(chrom, pos)          # already tiles (nearest-SNP map)
    if blockdef.startswith("clq"):
        r2 = float(blockdef[3:])
        if how == "tiling":
            import blocks_tiling as bt
            return bt.merge_small_blocks(bt.assign_tiling(chrom, pos, r2=r2))
        return lib.assign_clq_blocks(chrom, pos, r2=r2)
    raise ValueError(f"unknown blockdef {blockdef}")


def reblock_one(model, cls, axis, blockdef, outdir, how="strict"):
    rel, pcol = SRC[model]
    src = f"{MA}/{rel.format(cls=cls, axis=axis)}"
    if not os.path.exists(src):
        raise FileNotFoundError(src)
    df = pd.read_csv(src).rename(columns={pcol: "pval"})
    df = df.drop(columns=[c for c in ("block",) if c in df.columns])
    n0 = len(df)
    df["block"] = assign(blockdef, df["chrom"].to_numpy(str), df["pos"].to_numpy(np.int64), how)
    df = df[df["block"].to_numpy() != ""].copy()
    keep = [c for c in ("chrom", "pos", "ref_len", "alt_len", "MAF", "block", "pval") if c in df.columns]
    os.makedirs(outdir, exist_ok=True)
    out = f"{outdir}/{model}_{cls}_gen9_{axis}.csv"
    df[keep].to_csv(out, index=False)
    print(f"  {model:9s} {cls:11s} {axis:5s} [{blockdef}/{how}]: {n0:>9,} -> {len(df):>9,} in-block "
          f"({len(df)/n0*100:.1f}%) | {df['block'].nunique():,} blocks -> {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--axis", required=True)
    ap.add_argument("--cls", required=True, choices=["snp", "sv", "smallindel", "nonsnp"])
    ap.add_argument("--blockdef", required=True, choices=["hapfire", "clq0.5", "clq0.7", "clq0.9"])
    ap.add_argument("--how", default="strict", choices=["strict", "tiling"],
                    help="strict = interval containment (drops 40-49%% of records); "
                         "tiling = HapFM's gap-free conversion (drops none)")
    ap.add_argument("--models", nargs="+", default=["kendall", "lfmm", "binomial"])
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    tag = args.blockdef.replace('.', '') + ("_tile" if args.how == "tiling" else "")
    outdir = args.outdir or f"{MA}/wza_in_{tag}"
    for m in args.models:
        reblock_one(m, args.cls, args.axis, args.blockdef, outdir, args.how)


if __name__ == "__main__":
    main()

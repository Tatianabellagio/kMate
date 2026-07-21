#!/usr/bin/env python
"""Reblock the multi-axis per-record model outputs onto clq0.9 blocks (one axis+class).

Uniform over the 3 production models — kendall, lfmm, quasi-binomial (the production
binomial) — reading each model's multiaxis per-record CSV, reassigning the clq0.9
block via lib.assign_clq_blocks, dropping inter-block gaps, and writing the
WZA-input schema (chrom,pos,ref_len,alt_len,MAF,block,pval) to multiaxis/wza_in/.

The quasi-binomial file carries pval_binom/pval_quasi/pval_effN; we take pval_quasi
(the production choice) and write it as `pval`, and name the output `binomial_...`
so run_wza.py --model binomial picks it up (matching the production clq90 naming).

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY reblock_multiaxis.py --axis bio5 --cls snp --r2 0.9
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

MA = f"{lib.GEA}/phase1_replication/results/multiaxis"

# (source relpath template, pval column in source) per output model name
SRC = {
    "kendall":  ("kendall/kendall_{cls}_gen9_{axis}.csv", "pval"),
    "lfmm":     ("lfmm/lfmm_{cls}_gen9_{axis}.csv", "pval"),
    "binomial": ("quasibinom/quasibinom_lf16_{cls}_gen9_{axis}.csv", "pval_quasi"),
}


def reblock_one(model, cls, axis, r2):
    rel, pcol = SRC[model]
    src = f"{MA}/{rel.format(cls=cls, axis=axis)}"
    if not os.path.exists(src):
        print(f"  SKIP {model} {cls} {axis}: missing {src}", flush=True); return
    df = pd.read_csv(src)
    df = df.rename(columns={pcol: "pval"})
    df = df.drop(columns=[c for c in ("block",) if c in df.columns])
    n0 = len(df)
    df["block"] = lib.assign_clq_blocks(df["chrom"].to_numpy(str),
                                        df["pos"].to_numpy(np.int64), r2=r2)
    df = df[df["block"].to_numpy() != ""].copy()
    keep = [c for c in ("chrom", "pos", "ref_len", "alt_len", "MAF", "block", "pval") if c in df.columns]
    outdir = f"{MA}/wza_in"; os.makedirs(outdir, exist_ok=True)
    out = f"{outdir}/{model}_{cls}_gen9_{axis}.csv"
    df[keep].to_csv(out, index=False)
    print(f"  {model:9s} {cls:7s} {axis:5s}: {n0:>9,} -> {len(df):>9,} in-block | "
          f"{df['block'].nunique():,} blocks -> {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--axis", required=True)
    ap.add_argument("--cls", required=True, choices=["snp", "sv", "smallindel"])
    ap.add_argument("--models", nargs="+", default=["kendall", "lfmm", "binomial"])
    ap.add_argument("--r2", type=float, default=0.9)
    args = ap.parse_args()
    for m in args.models:
        reblock_one(m, args.cls, args.axis, args.r2)


if __name__ == "__main__":
    main()

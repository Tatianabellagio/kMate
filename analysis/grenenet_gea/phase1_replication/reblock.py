#!/usr/bin/env python
"""Re-block per-record GEA outputs onto the clq0.9 (r2>=0.9) BigLD blocks.

The clq0.9 replication changes ONLY the window definition (phase-1 hapFIRE blocks
-> finer clq0.9 LD islands) and the WZA aggregation. The per-record GEA statistics
(Kendall tau, LFMM calibrated p) are block-independent, so we reuse the existing
per-record outputs and simply reassign the `block` column here, rather than re-run
the models. (The binomial DOES re-run, because it switches to the site-effective-N
deflation — see run_clq90_binomial.sbatch — but its per-record CSV is still just
re-blocked here.)

For each (model, class) we:
  * load the source per-record CSV (chrom,pos,ref_len,alt_len,MAF,block,stat,pval),
  * drop the old (hapFIRE) block column,
  * assign the clq0.9 block via lib.assign_clq_blocks (interval map),
  * DROP records in inter-block gaps (~40%: clq0.9 blocks are LD islands, they do
    not tile the genome; the natural WZA unit is the defined block),
  * write clq90/wza_in/<model>_<class>_gen9_bio1.csv (a dedicated dir, NOT the raw
    source dir, so the raw binomial is never clobbered), ready for run_wza.py.

Class `nonsnp` (= SV + small-indel together) for kendall has no single source file,
so we concat the existing kendall_sv + kendall_smallindel per-record outputs (their
per-record tau/p are identical whether labelled sv/smallindel or nonsnp; verified
nonsnp records == sv + smallindel exactly).

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY reblock.py --gen 9 --climate bio1 --r2 0.9
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

BASE = f"{lib.GEA}/phase1_replication/results"
CLQ = f"{BASE}/clq90"

# where each (model, class) per-record source lives. binomial is the freshly-built
# site-effective-N run under clq90/; kendall/lfmm are reused from the phase-1-block run.
def source_paths(model: str, cls: str, gen: int, climate: str):
    if model == "binomial":
        return [f"{CLQ}/binomial/binomial_{cls}_gen{gen}_{climate}.csv"]
    if model == "lfmm":
        return [f"{BASE}/lfmm/lfmm_{cls}_gen{gen}_{climate}.csv"]
    if model == "kendall":
        if cls == "nonsnp":                       # concat sv + smallindel per-record
            return [f"{BASE}/kendall/kendall_sv_gen{gen}_{climate}.csv",
                    f"{BASE}/kendall/kendall_smallindel_gen{gen}_{climate}.csv"]
        return [f"{BASE}/kendall/kendall_{cls}_gen{gen}_{climate}.csv"]
    raise ValueError(model)


def reblock(model: str, cls: str, gen: int, climate: str, r2: float):
    paths = source_paths(model, cls, gen, climate)
    for p in paths:
        if not os.path.exists(p):
            print(f"  SKIP {model} {cls}: missing {p}", flush=True)
            return
    df = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)
    n0 = len(df)
    df = df.drop(columns=[c for c in ("block",) if c in df.columns])
    df["block"] = lib.assign_clq_blocks(df["chrom"].to_numpy(str),
                                        df["pos"].to_numpy(np.int64), r2=r2)
    ingap = df["block"].to_numpy() == ""
    df = df[~ingap].copy()
    # Reblocked outputs go to a dedicated dir (never the raw-source dir) so the raw
    # binomial (full, un-gap-dropped) is preserved for re-blocking at another r2.
    outdir = f"{CLQ}/wza_in"
    os.makedirs(outdir, exist_ok=True)
    out = f"{outdir}/{model}_{cls}_gen{gen}_{climate}.csv"
    df.to_csv(out, index=False)
    nblk = df["block"].nunique()
    print(f"  {model:9s} {cls:7s}: {n0:>9,} records -> {len(df):>9,} in-block "
          f"({100*ingap.mean():.0f}% dropped as gaps) | {nblk:,} clq{r2} blocks -> {out}",
          flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=["kendall", "lfmm", "binomial"])
    ap.add_argument("--classes", nargs="+", default=["snp", "nonsnp"])
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--r2", type=float, default=0.9)
    args = ap.parse_args()
    print(f"Re-blocking onto clq{args.r2} BigLD blocks (gen{args.gen}, {args.climate})", flush=True)
    for model in args.models:
        for cls in args.classes:
            reblock(model, cls, args.gen, args.climate, args.r2)


if __name__ == "__main__":
    main()

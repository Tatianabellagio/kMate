#!/usr/bin/env python
"""One-off: reblock the gea_newpanel effective-N (ACER) binomial per-record results
onto OUR clq0.9 partition (lib.assign_clq_blocks), so they're comparable to the
existing clq90 pipeline, and emit per-record diagnostics.

Source: analysis/grenenet_selection/archive/gea_newpanel_snp_nonsnp_fork_retired/results/quasibinom/quasibinom_lf16_{cls}_gen9_bio1.csv
  columns: chrom,pos,ref_len,alt_len,MAF,block,slope,pval_binom,pval_quasi,phi,pval_effN
  (the `block` column here is from a DIFFERENT clq0.9 build (~82k blocks) - IGNORE/drop it)

Output: analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/clq90/binom_fix_test/effn/binomial_{cls}_gen9_bio1.csv
  columns: chrom,pos,MAF,pval,block   (pval = pval_effN, renamed; block = OUR clq0.9 assignment)
  (named binomial_{cls}_... , not effn_{cls}_..., so run_wza.py --model binomial --indir .../effn works)

This is a SCRATCH/test script per the binom_fix_test investigation - does not touch
or overwrite anything under clq90/wza_in, clq90/binomial, clq90/wza, or compare_clq90.csv.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

SRC = f"{lib.GEA}/archive/gea_newpanel_snp_nonsnp_fork_retired/results/quasibinom"
# lib.GEA points at the analysis dir; results live in the mirrored results tree.
RESULTS_GEA = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection"
OUTDIR = f"{RESULTS_GEA}/phase1_replication/clq90/binom_fix_test/effn"
os.makedirs(OUTDIR, exist_ok=True)

GEN = 9
CLIMATE = "bio1"


def gif(p):
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p)]
    if len(p) == 0:
        return np.nan
    return np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)


def process(cls: str):
    src = f"{SRC}/quasibinom_lf16_{cls}_gen{GEN}_{CLIMATE}.csv"
    print(f"== {cls} == loading {src}", flush=True)
    df = pd.read_csv(src)
    n0 = len(df)

    # --- per-record diagnostics BEFORE any filtering/reblocking ---
    p_raw = df["pval_effN"].to_numpy(float)
    n_fin = np.isfinite(p_raw).sum()
    frac_p05 = np.mean(p_raw[np.isfinite(p_raw)] < 0.05)
    min_p = np.nanmin(p_raw)
    g = gif(p_raw)
    print(f"  raw per-record (n={n0:,}, finite={n_fin:,}): "
          f"frac p<0.05={frac_p05:.4f}  min p={min_p:.3e}  GIF-proxy={g:.3f}", flush=True)

    # --- MAF floor (should mostly already hold; enforce it) ---
    maf_ok = df["MAF"] >= 0.05
    n_maf_drop = (~maf_ok).sum()
    df = df[maf_ok].copy()
    print(f"  MAF>=0.05 filter: dropped {n_maf_drop:,}/{n0:,} ({100*n_maf_drop/n0:.2f}%)", flush=True)

    # --- rename pval_effN -> pval, drop the gea_newpanel block col ---
    df = df.drop(columns=[c for c in ("block",) if c in df.columns])
    df = df.rename(columns={"pval_effN": "pval"})

    # --- reassign block via OUR clq0.9 partition ---
    df["block"] = lib.assign_clq_blocks(df["chrom"].to_numpy(str),
                                         df["pos"].to_numpy(np.int64), r2=0.9)
    ingap = df["block"].to_numpy() == ""
    n_gap = ingap.sum()
    df = df[~ingap].copy()

    keep = df["pval"].notna()
    n_nanp = (~keep).sum()
    df = df[keep].copy()

    nblk = df["block"].nunique()
    print(f"  reblock (clq0.9, OUR partition): dropped {n_gap:,} gap records "
          f"({100*n_gap/(len(df)+n_gap):.1f}%), {n_nanp:,} NaN-pval records | "
          f"{len(df):,} in-block records -> {nblk:,} blocks", flush=True)

    out_df = df[["chrom", "pos", "MAF", "pval", "block"]].copy()
    out_path = f"{OUTDIR}/binomial_{cls}_gen{GEN}_{CLIMATE}.csv"
    out_df.to_csv(out_path, index=False)
    print(f"  -> {out_path}  ({len(out_df):,} rows)", flush=True)

    return dict(cls=cls, n0=n0, n_fin=n_fin, frac_p05=frac_p05, min_p=min_p, gif=g,
                n_maf_drop=n_maf_drop, n_gap=n_gap, n_blocks=nblk, n_final=len(out_df))


def main():
    stats_rows = [process(c) for c in ("snp", "nonsnp")]
    summ = pd.DataFrame(stats_rows)
    summ_path = f"{OUTDIR}/reblock_diagnostics.csv"
    summ.to_csv(summ_path, index=False)
    print(f"\nSummary -> {summ_path}")
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()

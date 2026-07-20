#!/usr/bin/env python
"""Run the WZA (../wza_script.py, canonical Booker) on a per-record GEA result.

Consumes a model output CSV (kendall / lfmm / binomial) that carries the columns
`block` (window), `pval` (per-record p, small = significant) and `MAF` (weight),
collapses records to one weighted-Z per LD block, applies the SNP-number
correction and emits empirical block-level Z_pVal.

We run TWO corrections:
  * deg-2  (primary)     = canonical Booker WZA  -> wza_{model}_{class}_gen{g}_{climate}.csv
  * deg-7  (sensitivity) = phase-1 variant       -> ..._deg7.csv
so we can confirm any hit (e.g. CAM5) is not an artifact of the correction order.

The per-record p-value is rank-transformed to an empirical p inside wza_script.py
(small input p -> small empirical p; we do NOT pass --large_i_small_p). MAF<0.05
records are already removed upstream by build_class_matrices, but we still pass
--maf_filter 0.05 to match phase-1 exactly. Rows with no/None block or NaN pval
are dropped before the call so the window grouping stays clean.

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY run_wza.py --model kendall --class sv --gen 1 --climate bio1
"""
from __future__ import annotations
import argparse, os, subprocess, sys, tempfile
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

WZA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "wza_script.py")
PY = sys.executable


def run(in_csv: str, out_csv: str, min_snps: int, verbose: bool,
        poly_deg: int | None = None, sample_snps: int | None = None,
        min_entries: int | None = None):
    # wza_script.py exposes the SNP-number-correction knobs (--poly_deg/--roller/
    # --min_entries) and a SNP cap (--sample_snps). Defaults = canonical Booker deg-2.
    # The phase-1 / last-gen runs (and the wza_investigation finding that deg-7 fits the
    # support curve better, RMSE 0.94 vs deg-2 2.06) use poly_deg=7; pass it through here.
    cmd = [PY, WZA,
           "--correlations", in_csv,
           "--summary_stat", "pval",
           "--window", "block",
           "--MAF", "MAF",
           "--maf_filter", "0.05",
           "--min_snps", str(min_snps),
           "--sep", ",",
           "--retain", "chrom", "pos",
           "--output", out_csv]
    if poly_deg is not None:
        cmd += ["--poly_deg", str(poly_deg)]
    if min_entries is not None:
        cmd += ["--min_entries", str(min_entries)]
    if sample_snps is not None:
        cmd += ["--sample_snps", str(sample_snps)]
    if verbose:
        cmd.append("-v")
    print("  $ " + " ".join(cmd), flush=True)
    # wza_script.py writes before_filtering_wza_df.csv / problematic_windows.csv to CWD;
    # run it inside the output dir so those debris land next to the result, not in repo root.
    subprocess.run(cmd, check=True, cwd=os.path.dirname(out_csv) or ".")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="kendall", choices=["kendall", "lfmm", "binomial", "betabinom"])
    ap.add_argument("--class", dest="cls", required=True,
                    choices=["snp", "sv", "smallindel", "nonsnp"])
    ap.add_argument("--gen", type=int, required=True)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--indir", default=None, help="dir of the model result CSV (default: ../<model>)")
    ap.add_argument("--out", default=f"{lib.GEA}/phase1_replication/results/wza")
    ap.add_argument("--min-snps", type=int, default=2, help="min records per block (WZA --min_snps)")
    ap.add_argument("--poly-deg", type=int, default=None,
                    help="SNP-number-correction degree (None=canonical deg-2; 7=phase-1/last-gen)")
    ap.add_argument("--sample-snps", type=int, default=None,
                    help="cap SNPs/window (0=no cap; e.g. 2000). None=script default")
    ap.add_argument("--min-entries", type=int, default=None,
                    help="min entries per rolling window (e.g. 10 for deg-7, 40 for deg-2)")
    ap.add_argument("--regime", default=None,
                    help="filename suffix tag, e.g. deg7nocap / deg7cap2000 (matches kendall/lfmm)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    indir = args.indir or f"{lib.GEA}/phase1_replication/results/{args.model}"
    in_csv = os.path.abspath(f"{indir}/{args.model}_{args.cls}_gen{args.gen}_{args.climate}.csv")
    if not os.path.exists(in_csv):
        sys.exit(f"missing model result: {in_csv}")
    # run() sets the WZA subprocess cwd to the output dir, so the output path MUST be
    # absolute (a relative --out would resolve against that new cwd and land nowhere).
    args.out = os.path.abspath(args.out)
    os.makedirs(args.out, exist_ok=True)

    df = pd.read_csv(in_csv)
    n0 = len(df)
    df = df[df["block"].notna() & (df["block"].astype(str) != "") &
            (df["block"].astype(str) != "None") & df["pval"].notna()].copy()
    # block must be a clean grouping key (avoid float '123.0' vs '123' mismatches)
    df["block"] = df["block"].astype(str).str.replace(r"\.0$", "", regex=True)
    nblk = df["block"].nunique()
    big = (df.groupby("block").size() >= args.min_snps).sum()
    print(f"{args.model} {args.cls} gen{args.gen} {args.climate}: {len(df):,}/{n0:,} records "
          f"with block+pval | {nblk:,} blocks ({big:,} with >={args.min_snps} records)", flush=True)

    with tempfile.NamedTemporaryFile("w", suffix=".csv", dir=args.out, delete=False) as tf:
        df.to_csv(tf.name, index=False)
        clean = tf.name
    try:
        stem = f"wza_{args.model}_{args.cls}_gen{args.gen}_{args.climate}"
        if args.regime:
            stem += f"_{args.regime}"
        out2 = f"{args.out}/{stem}.csv"
        run(clean, out2, min_snps=args.min_snps, verbose=args.verbose,
            poly_deg=args.poly_deg, sample_snps=args.sample_snps,
            min_entries=args.min_entries)
        _report(out2, args.regime or ("deg%d" % args.poly_deg if args.poly_deg else "deg2"))
    finally:
        os.unlink(clean)


def _report(out_csv: str, tag: str):
    w = pd.read_csv(out_csv)
    pcol = "Z_pVal" if "Z_pVal" in w.columns else None
    if pcol is None:
        print(f"  [{tag}] -> {out_csv} ({len(w):,} blocks; no Z_pVal col — raw scores)", flush=True)
        return
    w = w[w[pcol].notna()]
    nbon = (w[pcol] < 0.05 / len(w)).sum() if len(w) else 0
    top = w.nsmallest(5, pcol)
    gcol = "gene" if "gene" in w.columns else w.columns[0]
    print(f"  [{tag}] -> {out_csv} | {len(w):,} blocks | Bonferroni hits (p<0.05/{len(w)}): {nbon}",
          flush=True)
    for _, r in top.iterrows():
        ch = f" {r['chrom']}:{int(r['pos']):,}" if "chrom" in w.columns and pd.notna(r.get("pos")) else ""
        print(f"      block {r[gcol]}{ch}  Z_pVal={r[pcol]:.2e}", flush=True)


if __name__ == "__main__":
    main()

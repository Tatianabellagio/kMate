#!/usr/bin/env python
"""Per-axis WZA (clq0.9 blocks, deg7-cap2000) for all 3 models x 3 classes, reading the
per-record GEA produced by the multiaxis worker. One --axis at a time.

Reads  analysis/grenenet_gea/gea_newpanel/results/multiaxis_fresh/{binomial,kendall,lfmm}/
       {model}_{cls}_gen9_{axis}.csv   (chrom,pos,...,MAF,pval[,block])
Writes analysis/grenenet_gea/gea_newpanel/results/multiaxis_fresh/wza_{model}_{cls}_{axis}_clq09.csv
"""
from __future__ import annotations
import os, sys, subprocess, tempfile, argparse
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blocks_clq09

PROJ = "/global/scratch/users/tbellg/kmate"
MA = f"{PROJ}/analysis/grenenet_gea/gea_newpanel/results/multiaxis_fresh"
WZA = f"{PROJ}/analysis/grenenet_gea/wza_script.py"
PY = sys.executable
MODELS = ["binomial", "kendall", "lfmm"]
CLASSES = ["snp", "nonsnp", "sv"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--axis", required=True)
    a = ap.parse_args()
    for model in MODELS:
        for cls in CLASSES:
            src = f"{MA}/{model}/{model}_{cls}_gen9_{a.axis}.csv"
            if not os.path.exists(src):
                print(f"  SKIP {model} {cls} {a.axis}: missing {src}", flush=True)
                continue
            d = pd.read_csv(src)
            if "pval" not in d.columns or "MAF" not in d.columns:
                print(f"  SKIP {model} {cls} {a.axis}: cols {list(d.columns)}", flush=True)
                continue
            d["block"] = blocks_clq09.assign_clq09_blocks(d.chrom.to_numpy(), d.pos.to_numpy())
            rec = f"{MA}/{model}_{cls}_gen9_{a.axis}_clq09.csv"
            d.to_csv(rec, index=False)
            wout = f"{MA}/wza_{model}_{cls}_{a.axis}_clq09.csv"
            with tempfile.TemporaryDirectory(prefix=f"wza_{model}_{cls}_{a.axis}_") as td:
                subprocess.run([PY, WZA, "--correlations", rec, "--summary_stat", "pval",
                                "--window", "block", "--MAF", "MAF", "--maf_filter", "0.05",
                                "--sep", ",", "--retain", "chrom", "pos",
                                "--poly_deg", "7", "--min_entries", "10", "--sample_snps", "2000",
                                "--output", wout], cwd=td, check=True)
            w = pd.read_csv(wout); ok = w[w["Z_pVal"].notna()]
            print(f"  [{model} {cls} {a.axis}] WZA {len(w):,} blk, min p={ok['Z_pVal'].min():.1e}", flush=True)


if __name__ == "__main__":
    main()

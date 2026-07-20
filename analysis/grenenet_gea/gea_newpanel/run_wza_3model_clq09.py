#!/usr/bin/env python
"""Reassign clq0.9 haploblocks to the FRESH (rerun_kfw_hb cohort) per-record
climate-GEA p-values for ALL THREE models — binomial, kendall, lfmm — and run the
canonical WZA (deg7-cap2000) per variant class. Generalizes run_wza_clq09.py
(which did LFMM only) so the unified SNP-vs-nonSNP 3x2 raw+WZA grid can be built on
one consistent block partition (blocks_clq09 -> blocks_recompute clq0.9, nearest-gap
assignment), matching the already-fresh LFMM-WZA byte-for-byte on the LFMM row.

Raw per-record sources (all regenerated 2026-07-08 from the rerun_kfw_hb cohort;
each carries chrom,pos,ref_len,alt_len,MAF,block,<stat>,pval):
  binomial : gea_newpanel/binomial/binomial_{cls}_gen9_bio1.csv   (plain full-N GLM, pval)
  kendall  : gea_newpanel/kendall/kendall_{cls}_gen9_bio1.csv     (raw tau, pval)
  lfmm     : phase1_replication/lfmm/lfmm_{cls}_gen9_bio1.csv      (calibrated LFMM p)

Outputs under analysis/grenenet_gea/gea_newpanel/results/:
  {model}_{cls}_gen9_bio1_clq09.csv   per-record, clq0.9 block reassigned (raw view)
  wza_{model}_{cls}_clq09.csv         per-block WZA (block, Z_pVal, chrom, pos)
"""
from __future__ import annotations
import os, sys, subprocess, tempfile
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blocks_clq09

PROJ = "/global/scratch/users/tbellg/kmate"
GNP = f"{PROJ}/analysis/grenenet_gea/gea_newpanel/results"
LFMM_DIR = f"{PROJ}/analysis/grenenet_gea/phase1_replication/results/lfmm"
WZA = f"{PROJ}/analysis/grenenet_gea/wza_script.py"
PY = sys.executable

SOURCES = {
    "binomial": lambda cls: f"{GNP}/binomial/binomial_{cls}_gen9_bio1.csv",
    "kendall":  lambda cls: f"{GNP}/kendall/kendall_{cls}_gen9_bio1.csv",
    "lfmm":     lambda cls: f"{LFMM_DIR}/lfmm_{cls}_gen9_bio1.csv",
}


def main():
    os.makedirs(GNP, exist_ok=True)
    for model, src in SOURCES.items():
        for cls in ["snp", "nonsnp"]:
            in_csv = src(cls)
            if not os.path.exists(in_csv):
                sys.exit(f"MISSING raw source: {in_csv}")
            d = pd.read_csv(in_csv)
            if "pval" not in d.columns or "MAF" not in d.columns:
                sys.exit(f"{in_csv}: need pval+MAF cols, got {list(d.columns)}")
            d["block"] = blocks_clq09.assign_clq09_blocks(d.chrom.to_numpy(), d.pos.to_numpy())
            rec_out = f"{GNP}/{model}_{cls}_gen9_bio1_clq09.csv"
            d.to_csv(rec_out, index=False)
            n_blk = d.loc[d.block != "", "block"].nunique()
            print(f"[{model} {cls}] {len(d):,} records -> {n_blk:,} clq0.9 blocks "
                  f"| raw min p={d.pval.min():.1e} | {rec_out}", flush=True)

            # canonical WZA deg7-cap2000 (identical call to run_wza_clq09.py)
            wout = f"{GNP}/wza_{model}_{cls}_clq09.csv"
            with tempfile.TemporaryDirectory(prefix=f"wza_{model}_{cls}_") as td:
                cmd = [PY, WZA, "--correlations", rec_out, "--summary_stat", "pval",
                       "--window", "block", "--MAF", "MAF", "--maf_filter", "0.05",
                       "--sep", ",", "--retain", "chrom", "pos",
                       "--poly_deg", "7", "--min_entries", "10", "--sample_snps", "2000",
                       "--output", wout]
                subprocess.run(cmd, cwd=td, check=True)
            w = pd.read_csv(wout)
            nan = int(w["Z_pVal"].isna().sum())
            ok = w[w["Z_pVal"].notna()]
            bonf = int((ok["Z_pVal"] < 0.05 / max(len(ok), 1)).sum())
            print(f"[{model} {cls}] WZA -> {wout} | {len(w):,} blocks | NaN={nan} "
                  f"| Bonf<0.05: {bonf} | min Z_pVal={ok['Z_pVal'].min():.2e}\n", flush=True)


if __name__ == "__main__":
    main()

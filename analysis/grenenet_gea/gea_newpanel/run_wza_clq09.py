#!/usr/bin/env python
"""Reassign clq0.9 haploblocks to the raw-LFMM per-record p-values and run the
canonical WZA (deg7-cap2000) per variant class, for the new-panel SNP-vs-nonSNP GEA.

For cls in {snp, nonsnp}:
  1. load results/.../phase1_replication/lfmm/lfmm_{cls}_gen9_bio1.csv
     (chrom,pos,ref_len,alt_len,MAF,block[phase1],pval  -- pval = GIF-calibrated LFMM p)
  2. replace `block` with the clq0.9 block (blocks_clq09.assign_clq09_blocks)
  3. write the reassigned per-record table (for the notebook's raw Manhattan) and a
     WZA `correlations` input, then run wza_script.py --window block deg7-cap2000.

Outputs under results/grenenet_gea/gea_newpanel/:
  lfmm_{cls}_gen9_bio1_clq09.csv   per-record, clq0.9 block  (raw LFMM view)
  wza_{cls}_clq09.csv              per-block WZA (gene=block, Z_pVal, chrom, pos)
"""
from __future__ import annotations
import os, sys, subprocess, tempfile
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blocks_clq09

PROJ = "/global/scratch/users/tbellg/kmate"
LDIR = f"{PROJ}/results/grenenet_gea/phase1_replication/lfmm"      # raw LFMM p (reused)
OUT = f"{PROJ}/results/grenenet_gea/gea_newpanel"
WZA = f"{PROJ}/analysis/grenenet_gea/wza_script.py"
PY = sys.executable


def main():
    os.makedirs(OUT, exist_ok=True)
    for cls in ["snp", "nonsnp"]:
        d = pd.read_csv(f"{LDIR}/lfmm_{cls}_gen9_bio1.csv")
        d["block"] = blocks_clq09.assign_clq09_blocks(d.chrom.to_numpy(), d.pos.to_numpy())
        d["sv_size"] = (d.alt_len - d.ref_len).abs()
        rec_out = f"{OUT}/lfmm_{cls}_gen9_bio1_clq09.csv"
        d.to_csv(rec_out, index=False)
        n_blk = d.loc[d.block != "", "block"].nunique()
        print(f"[{cls}] {len(d):,} records -> {n_blk:,} clq0.9 blocks | {rec_out}", flush=True)

        # WZA (canonical deg7-cap2000), run in a temp dir (writes side CSVs to CWD)
        wout = f"{OUT}/wza_{cls}_clq09.csv"
        with tempfile.TemporaryDirectory(prefix=f"wza_{cls}_") as td:
            cmd = [PY, WZA, "--correlations", rec_out, "--summary_stat", "pval",
                   "--window", "block", "--MAF", "MAF", "--maf_filter", "0.05",
                   "--sep", ",", "--retain", "chrom", "pos",
                   "--poly_deg", "7", "--min_entries", "10", "--sample_snps", "2000",
                   "--output", wout]
            subprocess.run(cmd, cwd=td, check=True)
        w = pd.read_csv(wout)
        print(f"[{cls}] WZA -> {wout} | {len(w):,} blocks | "
              f"min Z_pVal={w['Z_pVal'].min():.2e}\n", flush=True)


if __name__ == "__main__":
    main()

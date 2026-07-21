#!/usr/bin/env python
"""WZA (clq0.9 blocks, deg7-cap2000) on the FRESH quasi-binomial+K16 p-values
(pval_quasi) — the settled binomial-inflation fix — so the 'fixed' grid_3model can
show the corrected binomial row. Mirrors run_wza_3model_clq09.py exactly, but reads
quasibinom_lf16[pval_quasi] as the per-record p.

Outputs under analysis/grenenet_gea/gea_newpanel/results/:
  quasibinom_{cls}_gen9_bio1_clq09.csv   per-record (pval = pval_quasi), clq0.9 block
  wza_quasibinom_{cls}_clq09.csv         per-block WZA
"""
from __future__ import annotations
import os, sys, subprocess, tempfile
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blocks_clq09

PROJ = "/global/scratch/users/tbellg/kmate"
GNP = f"{PROJ}/analysis/grenenet_gea/gea_newpanel/results"
WZA = f"{PROJ}/analysis/grenenet_gea/wza_script.py"
PY = sys.executable


def main():
    for cls in ["snp", "nonsnp"]:
        src = f"{GNP}/quasibinom/quasibinom_lf16_{cls}_gen9_bio1.csv"
        d = pd.read_csv(src, usecols=["chrom", "pos", "ref_len", "alt_len", "MAF", "pval_quasi"])
        d = d.rename(columns={"pval_quasi": "pval"})
        d["block"] = blocks_clq09.assign_clq09_blocks(d.chrom.to_numpy(), d.pos.to_numpy())
        rec_out = f"{GNP}/quasibinom_{cls}_gen9_bio1_clq09.csv"
        d.to_csv(rec_out, index=False)
        n_blk = d.loc[d.block != "", "block"].nunique()
        print(f"[quasibinom {cls}] {len(d):,} records -> {n_blk:,} clq0.9 blocks "
              f"| raw min p={d.pval.min():.1e} | {rec_out}", flush=True)
        wout = f"{GNP}/wza_quasibinom_{cls}_clq09.csv"
        with tempfile.TemporaryDirectory(prefix=f"wza_quasi_{cls}_") as td:
            cmd = [PY, WZA, "--correlations", rec_out, "--summary_stat", "pval",
                   "--window", "block", "--MAF", "MAF", "--maf_filter", "0.05",
                   "--sep", ",", "--retain", "chrom", "pos",
                   "--poly_deg", "7", "--min_entries", "10", "--sample_snps", "2000",
                   "--output", wout]
            subprocess.run(cmd, cwd=td, check=True)
        w = pd.read_csv(wout)
        nan = int(w["Z_pVal"].isna().sum()); ok = w[w["Z_pVal"].notna()]
        bonf = int((ok["Z_pVal"] < 0.05 / max(len(ok), 1)).sum())
        print(f"[quasibinom {cls}] WZA -> {wout} | {len(w):,} blocks | NaN={nan} "
              f"| Bonf<0.05: {bonf} | min Z_pVal={ok['Z_pVal'].min():.2e}\n", flush=True)


if __name__ == "__main__":
    main()

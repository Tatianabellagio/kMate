#!/usr/bin/env python
"""WZA (clq0.9 blocks, deg7-cap2000) on the SV-only class for all 3 GEA models, so
the class-gain tables can compare SNP vs non-SNP vs SV. Mirrors run_wza_3model_clq09.py.

Raw per-record SV sources (bio1, gen9; all fresh on rerun_kfw_hb):
  binomial : phase1_replication/binomial/binomial_sv_gen9_bio1.csv
  kendall  : phase1_replication/kendall/kendall_sv_gen9_bio1.csv
  lfmm     : phase1_replication/lfmm/lfmm_sv_gen9_bio1.csv
Outputs: gea_newpanel/{model}_sv_gen9_bio1_clq09.csv + wza_{model}_sv_clq09.csv
"""
from __future__ import annotations
import os, sys, subprocess, tempfile
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blocks_clq09

PROJ = "/global/scratch/users/tbellg/kmate"
GNP = f"{PROJ}/analysis/grenenet_gea/gea_newpanel/results"
P1 = f"{PROJ}/analysis/grenenet_gea/phase1_replication/results"
WZA = f"{PROJ}/analysis/grenenet_gea/wza_script.py"
PY = sys.executable
SRC = {"binomial": f"{P1}/binomial/binomial_sv_gen9_bio1.csv",
       "kendall":  f"{P1}/kendall/kendall_sv_gen9_bio1.csv",
       "lfmm":     f"{P1}/lfmm/lfmm_sv_gen9_bio1.csv"}


def main():
    for model, src in SRC.items():
        if not os.path.exists(src):
            sys.exit(f"MISSING sv source: {src}")
        d = pd.read_csv(src)
        if "pval" not in d.columns:
            sys.exit(f"{src}: no pval col ({list(d.columns)})")
        d["block"] = blocks_clq09.assign_clq09_blocks(d.chrom.to_numpy(), d.pos.to_numpy())
        rec_out = f"{GNP}/{model}_sv_gen9_bio1_clq09.csv"
        d.to_csv(rec_out, index=False)
        print(f"[{model} sv] {len(d):,} records -> {d.loc[d.block!='','block'].nunique():,} blocks "
              f"| raw min p={d.pval.min():.1e}", flush=True)
        wout = f"{GNP}/wza_{model}_sv_clq09.csv"
        with tempfile.TemporaryDirectory(prefix=f"wza_{model}_sv_") as td:
            subprocess.run([PY, WZA, "--correlations", rec_out, "--summary_stat", "pval",
                            "--window", "block", "--MAF", "MAF", "--maf_filter", "0.05",
                            "--sep", ",", "--retain", "chrom", "pos",
                            "--poly_deg", "7", "--min_entries", "10", "--sample_snps", "2000",
                            "--output", wout], cwd=td, check=True)
        w = pd.read_csv(wout); ok = w[w["Z_pVal"].notna()]
        print(f"[{model} sv] WZA -> {wout} | {len(w):,} blocks | NaN={int(w['Z_pVal'].isna().sum())} "
              f"| Bonf<0.05: {int((ok['Z_pVal']<0.05/max(len(ok),1)).sum())} | min Z_pVal={ok['Z_pVal'].min():.2e}\n", flush=True)


if __name__ == "__main__":
    main()

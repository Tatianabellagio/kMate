#!/usr/bin/env python3
"""DEPRECATED / KNOWN-WRONG — do not use. See benchmarks/SCORING_RULES.md.

This scorer joins truth<->competitor on (chrom,pos) ONLY, which pairs the WRONG allele
at multiallelic positions (truth/kMate carry no REF/ALT, only ref_len/alt_len) and
corrupts R² (e.g. kMate n1_g0 AF 0.937 here vs the correct 1.0000). It also scores the
MAR truth, a closed loop that flatters kMate (RULE 2).

USE INSTEAD: accuracy_vs_competitors/scripts/score_snp_fair.py (shared-panel,
allele-safe, info>=0.9 / fully-called basis). Left in place only for provenance.
"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.stderr.write(
    "\n*** score_competitor_snp.py is DEPRECATED/KNOWN-WRONG (chrom,pos-only join + MAR "
    "closed loop). Use score_snp_fair.py. See benchmarks/SCORING_RULES.md. ***\n"
    "Re-run with KMATE_ALLOW_BAD_SCORER=1 only if you really mean it.\n\n")
import os as _os
if _os.environ.get("KMATE_ALLOW_BAD_SCORER") != "1":
    sys.exit("refusing to run a known-wrong scorer (set KMATE_ALLOW_BAD_SCORER=1 to override)")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_benchmark_table import parse_pool, metrics, COLS  # noqa: E402


def norm_chrom(s):
    return s.astype(str).str.replace("Chr", "", regex=False).str.replace("chr", "", regex=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True)          # hapFIRE | vg-giraffe
    ap.add_argument("--pool", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--est", required=True, help="3-col TSV: chrom pos freq (no header)")
    ap.add_argument("--truth-col", default="truth_af")
    ap.add_argument("--info-min", type=float, default=0.9,
                    help="score only well-informed SNPs (truth info>=this) so MAR==physical "
                         "-- the fair, convention-free basis (BENCHMARK_DESIGN.md 2.1). 0=all.")
    ap.add_argument("--table", required=True)
    a = ap.parse_args()

    meta = parse_pool(a.pool)
    tr = pd.read_csv(a.truth, sep="\t").dropna(subset=[a.truth_col])
    tr = tr[(tr.ref_len == 1) & (tr.alt_len == 1)].copy()      # SNPs only
    if a.info_min > 0 and "info" in tr.columns:
        tr = tr[tr["info"] >= a.info_min]
    tr["k"] = norm_chrom(tr["chrom"]) + ":" + tr["pos"].astype(str)
    est = pd.read_csv(a.est, sep="\t", header=None, names=["chrom", "pos", "est"], na_values=["nan"])
    est["est"] = pd.to_numeric(est["est"], errors="coerce")
    est["k"] = norm_chrom(est["chrom"]) + ":" + est["pos"].astype(str)
    est = est.groupby("k", as_index=False).agg(est=("est", "first"))

    j = tr.merge(est, on="k", how="inner").dropna(subset=["est", a.truth_col])
    mt = metrics(j["est"].values.astype(float), j[a.truth_col].values.astype(float))
    row = {**meta, "tool": a.tool, "mode": "aligned", "var_class": "SNP", **mt,
           "truth_col": a.truth_col, "est_file": Path(a.est).name}
    print(f"[{a.tool}] {a.pool}: n={mt['n']:,} MAE={mt['MAE']:.4f} RMSE={mt['RMSE']:.4f} "
          f"R2={mt['R2']:.4f} r={mt['pearson_r']:.4f}", file=sys.stderr)

    out = pd.DataFrame([row])[COLS]
    tp = Path(a.table)
    out.to_csv(tp, sep="\t", mode="a", header=not tp.exists(), index=False)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Append rows to the master kMate benchmark table (benchmarks/benchmark_table.tsv).

One row per (run, variant-class). Metadata (panel, coverage, founders, generation,
mating, selection, seed) is parsed from the pool name; tool/mode are passed in.
Metrics: MAE / RMSE / R2 / Pearson r / bias / outlier-rate, per class.

Join est<->truth on (chrom,pos,ref_len,alt_len) PLUS a within-key occurrence index,
so duplicate-key SVs (two ALT sequences of equal length at one position) are paired
k-th-to-k-th in panel-meta order -- NOT cross-joined (the bug that crushed SV R2) and
NOT silently deduped (which under-counts SVs). Both files are emitted in the same
panel-meta order, so occurrence-index pairing is exact.

Usage:
  build_benchmark_table.py --tool kMate --mode block \
      --pool cov10_n50_g3_s42_self97_hotspots_dom500_p80_chr1 \
      --truth benchmarks/p80/sims/<pool>/recomb_truth.tsv.gz \
      --est   <kmate_af.tsv> [--truth-col truth_af] [--table benchmarks/benchmark_table.tsv]
"""
import argparse, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

KEYS = ["chrom", "pos", "ref_len", "alt_len"]
COLS = ["tool", "panel", "mode", "coverage", "n_founders", "generation", "mating",
        "selection", "seed", "var_class", "n", "MAE", "RMSE", "R2", "pearson_r",
        "bias", "outlier", "truth_col", "est_file"]


def parse_pool(pool):
    def g(pat, default=None, cast=str):
        m = re.search(pat, pool)
        return cast(m.group(1)) if m else default
    panel = "p80" if "p80" in pool else ("p231" if "p231" in pool else "NA")
    nf = g(r"_n(\d+)_", cast=int)
    if nf is None:  # cactheavy / pgheavy special pools
        nf = g(r"_n(\d+)", cast=int)
    sel = "dom500nr" if "dom500nr" in pool else ("dom500" if "dom500" in pool else "balanced")
    return dict(panel=panel,
                coverage=g(r"cov(\d+)", cast=int),
                n_founders=nf,
                generation=g(r"_g(\d+)", cast=int),
                mating="self97" if "self97" in pool else "outcross",
                selection=sel,
                seed=g(r"_s(\d+)", cast=int))


def var_class(ref_len, alt_len):
    return np.where(np.maximum(ref_len, alt_len) >= 50, "SV",
                    np.where((ref_len == 1) & (alt_len == 1), "SNP", "indel"))


def occ_key(df):
    """Add an occurrence index within each (chrom,pos,ref_len,alt_len) group, in row
    order, so duplicate-key rows get distinct join keys (0,1,2,...)."""
    return df.groupby(KEYS).cumcount()


def metrics(est, truth):
    d = est - truth
    n = len(d)
    if n == 0:
        return dict(n=0, MAE=np.nan, RMSE=np.nan, R2=np.nan, pearson_r=np.nan, bias=np.nan, outlier=np.nan)
    ss_tot = np.sum((truth - truth.mean()) ** 2)
    r2 = 1 - np.sum(d ** 2) / ss_tot if ss_tot > 0 else np.nan
    r = np.corrcoef(est, truth)[0, 1] if n > 1 and est.std() > 0 and truth.std() > 0 else np.nan
    return dict(n=n, MAE=float(np.abs(d).mean()), RMSE=float(np.sqrt((d ** 2).mean())),
                R2=float(r2), pearson_r=float(r), bias=float(d.mean()),
                outlier=float((np.abs(d) > 0.10).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True)
    ap.add_argument("--mode", required=True,
                     choices=["global", "block", "chrom", "ld"],
                     help="global/block are the pre-2026-07-07 --block-mode labels "
                          "(kept for old rows' provenance); chrom/ld are the current "
                          "--unit labels (chrom = selfing production, ld = per-block control)")
    ap.add_argument("--pool", required=True, help="sim pool name (metadata parsed from it)")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--est", required=True)
    ap.add_argument("--truth-col", default="truth_af")
    ap.add_argument("--table", default=str(Path(__file__).resolve().parent / "benchmark_table.tsv"))
    a = ap.parse_args()

    meta = parse_pool(a.pool)
    truth = pd.read_csv(a.truth, sep="\t").dropna(subset=[a.truth_col])
    truth["occ"] = occ_key(truth)
    est = pd.read_csv(a.est, sep="\t")
    est["occ"] = occ_key(est)
    m = truth.merge(est[KEYS + ["occ", "alt_freq"]], on=KEYS + ["occ"], how="inner")
    m["vclass"] = var_class(m["ref_len"].values, m["alt_len"].values)
    print(f"[{a.tool}/{a.mode}] {a.pool}: joined {len(m):,}/{len(truth):,} truth records",
          file=sys.stderr)

    rows = []
    for cls in ["ALL", "SNP", "indel", "SV"]:
        sub = m if cls == "ALL" else m[m.vclass == cls]
        if len(sub) == 0:
            continue
        mt = metrics(sub["alt_freq"].values.astype(float), sub[a.truth_col].values.astype(float))
        rows.append({**meta, "tool": a.tool, "mode": a.mode, "var_class": cls,
                     **mt, "truth_col": a.truth_col, "est_file": Path(a.est).name})
        print(f"  {cls:6} n={mt['n']:>9,} MAE={mt['MAE']:.4f} RMSE={mt['RMSE']:.4f} "
              f"R2={mt['R2']:.4f} r={mt['pearson_r']:.4f}", file=sys.stderr)

    out = pd.DataFrame(rows)[COLS]
    tp = Path(a.table)
    if tp.exists():
        out.to_csv(tp, sep="\t", mode="a", header=False, index=False)
    else:
        out.to_csv(tp, sep="\t", index=False)
    print(f"  -> appended {len(out)} rows to {tp}", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Score one kMate AF run against sim truth, NaN-aware, for the local-only benchmark.

Unlike build_benchmark_table.py, local-only (global-free window) mode emits NaN
alt_freq wherever a window is below the k-mer floor / empty / a record is
unassigned. That is the WHOLE point of the mode: it abstains instead of falling
back to the chrom-wide h. So accuracy and coverage are two different axes and
must be reported separately:

  * call_rate  = fraction of joined records with a FINITE est alt_freq
  * MAE/RMSE/R2/r/bias/outlier = computed ONLY on records where BOTH the truth
    and the estimate are finite ("where it speaks, is it right?")

Join est<->truth on (chrom,pos,ref_len,alt_len) + within-key occurrence index,
exactly as build_benchmark_table.py (duplicate-key SVs paired k-th-to-k-th, not
cross-joined). One row per (run, var_class) appended to a rows TSV.

Usage:
  score_localonly.py --pool <pool> --mode {global,anchored,localonly} \
      --unit {none,dynldK500,w10kb} --truth <recomb_truth_raw.tsv.gz> \
      --est <kmate_af.tsv> --table <rows.tsv>
"""
import argparse, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

KEYS = ["chrom", "pos", "ref_len", "alt_len"]
COLS = ["pool", "panel", "n_founders", "generation", "mating", "selection", "seed",
        "mode", "unit", "var_class", "n_total", "n_called", "call_rate",
        "MAE", "RMSE", "R2", "pearson_r", "bias", "outlier", "truth_col", "est_file"]


def parse_pool(pool):
    def g(pat, cast=str, default=None):
        m = re.search(pat, pool)
        return cast(m.group(1)) if m else default
    nf = g(r"_n(\d+)", cast=int)
    sel = "dom500nr" if "dom500nr" in pool else ("dom500" if "dom500" in pool else "balanced")
    return dict(panel=("p231" if "p231" in pool else ("p80" if "p80" in pool else "NA")),
                n_founders=nf,
                generation=g(r"_g(\d+)", cast=int),
                mating="self97" if "self97" in pool else "outcross",
                selection=sel,
                seed=g(r"_s(\d+)", cast=int))


def var_class(ref_len, alt_len):
    return np.where(np.maximum(ref_len, alt_len) >= 50, "SV",
                    np.where((ref_len == 1) & (alt_len == 1), "SNP", "indel"))


def metrics(est, truth):
    """est/truth: arrays already restricted to finite-both records."""
    d = est - truth
    n = len(d)
    if n == 0:
        return dict(MAE=np.nan, RMSE=np.nan, R2=np.nan, pearson_r=np.nan,
                    bias=np.nan, outlier=np.nan)
    ss_tot = np.sum((truth - truth.mean()) ** 2)
    r2 = 1 - np.sum(d ** 2) / ss_tot if ss_tot > 0 else np.nan
    r = np.corrcoef(est, truth)[0, 1] if n > 1 and est.std() > 0 and truth.std() > 0 else np.nan
    return dict(MAE=float(np.abs(d).mean()), RMSE=float(np.sqrt((d ** 2).mean())),
                R2=float(r2), pearson_r=float(r), bias=float(d.mean()),
                outlier=float((np.abs(d) > 0.10).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--mode", required=True, choices=["global", "anchored", "localonly"])
    ap.add_argument("--unit", required=True, choices=["none", "dynldK500", "w10kb"])
    ap.add_argument("--truth", required=True)
    ap.add_argument("--est", required=True)
    ap.add_argument("--truth-col", default="truth_af")
    ap.add_argument("--table", required=True)
    a = ap.parse_args()

    meta = parse_pool(a.pool)
    truth = pd.read_csv(a.truth, sep="\t").dropna(subset=[a.truth_col])
    truth["occ"] = truth.groupby(KEYS).cumcount()
    est = pd.read_csv(a.est, sep="\t")
    est["occ"] = est.groupby(KEYS).cumcount()
    m = truth.merge(est[KEYS + ["occ", "alt_freq"]], on=KEYS + ["occ"], how="inner")
    m["vclass"] = var_class(m["ref_len"].values, m["alt_len"].values)
    est_finite = np.isfinite(m["alt_freq"].values.astype(float))
    print(f"[{a.mode}/{a.unit}] {a.pool}: joined {len(m):,}/{len(truth):,} truth; "
          f"finite est {est_finite.sum():,} ({100*est_finite.mean():.1f}%)", file=sys.stderr)

    rows = []
    for cls in ["ALL", "SNP", "indel", "SV"]:
        sub = m if cls == "ALL" else m[m.vclass == cls]
        n_total = len(sub)
        if n_total == 0:
            continue
        e = sub["alt_freq"].values.astype(float)
        t = sub[a.truth_col].values.astype(float)
        fin = np.isfinite(e)
        n_called = int(fin.sum())
        mt = metrics(e[fin], t[fin])
        rows.append({**meta, "pool": a.pool, "mode": a.mode, "unit": a.unit,
                     "var_class": cls, "n_total": n_total, "n_called": n_called,
                     "call_rate": n_called / n_total if n_total else np.nan,
                     **mt, "truth_col": a.truth_col, "est_file": Path(a.est).name})
        print(f"  {cls:6} n={n_total:>9,} called={n_called:>9,} "
              f"({100*n_called/n_total:5.1f}%) MAE={mt['MAE']:.4f} R2={mt['R2']:.4f} "
              f"r={mt['pearson_r']:.4f}", file=sys.stderr)

    out = pd.DataFrame(rows)[COLS]
    tp = Path(a.table)
    if tp.exists():
        out.to_csv(tp, sep="\t", mode="a", header=False, index=False)
    else:
        out.to_csv(tp, sep="\t", index=False)
    print(f"  -> appended {len(out)} rows to {tp}", file=sys.stderr)


if __name__ == "__main__":
    main()

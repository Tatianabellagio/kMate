#!/usr/bin/env python3
"""DEPRECATED / KNOWN-WRONG — do not use. See benchmarks/SCORING_RULES.md.

Joins all tools on (chrom,pos) ONLY (load_kmate/load_3col group by norm_chrom:pos),
which mispairs alleles at multiallelic positions and corrupts R²; also defaults to
truth_af_phys which isn't even present in recomb_truth. USE INSTEAD:
accuracy_vs_competitors/scripts/score_snp_fair.py. Left only for provenance.
"""
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.stderr.write(
    "\n*** score_competitors.py is DEPRECATED/KNOWN-WRONG (chrom,pos-only join). "
    "Use score_snp_fair.py. See benchmarks/SCORING_RULES.md. ***\n"
    "Re-run with KMATE_ALLOW_BAD_SCORER=1 only if you really mean it.\n\n")
import os as _os
if _os.environ.get("KMATE_ALLOW_BAD_SCORER") != "1":
    sys.exit("refusing to run a known-wrong scorer (set KMATE_ALLOW_BAD_SCORER=1 to override)")

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # benchmarks/
from _accuracy_panel import grid_panel  # noqa: E402


def norm_chrom(s):
    return s.astype(str).str.replace("Chr", "", regex=False).str.replace("chr", "", regex=False)


def load_truth(p):
    t = pd.read_csv(p, sep="\t")
    t["k"] = norm_chrom(t["chrom"]) + ":" + t["pos"].astype(str)
    return t


def load_kmate(p):
    d = pd.read_csv(p, sep="\t")
    d = d[(d.ref_len == 1) & (d.alt_len == 1)].copy()   # SNPs only
    d["k"] = norm_chrom(d["chrom"]) + ":" + d["pos"].astype(str)
    return d.groupby("k", as_index=False).agg(est=("alt_freq", "first"))


def load_3col(p):
    """chrom, pos, freq (hapFIRE _snp_frequency.txt and vg AF output)."""
    d = pd.read_csv(p, sep="\t", header=None, names=["chrom", "pos", "est"],
                    na_values=["nan"])
    d["est"] = pd.to_numeric(d["est"], errors="coerce")
    d["k"] = norm_chrom(d["chrom"]) + ":" + d["pos"].astype(str)
    return d.groupby("k", as_index=False).agg(est=("est", "first"))


def metrics(truth, est):
    m = np.isfinite(truth) & np.isfinite(est)
    truth, est = truth[m], est[m]
    err = est - truth
    mae = float(np.abs(err).mean())
    rmse = float(np.sqrt((err ** 2).mean()))
    ss_res = float((err ** 2).sum())
    ss_tot = float(((truth - truth.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    pear = float(np.corrcoef(truth, est)[0, 1]) if len(truth) > 1 else np.nan
    return dict(n=int(m.sum()), MAE=mae, RMSE=rmse, R2=r2, pearson_r=pear)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--kmate")
    ap.add_argument("--hapfire")
    ap.add_argument("--vg")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--F", type=int, default=80, help="#founders (for miss filter)")
    ap.add_argument("--miss-thr", type=float, default=None,
                    help="keep truth sites with n_called >= (1-thr)*F")
    ap.add_argument("--truth-col", default="truth_af_phys",
                    choices=["truth_af_phys", "truth_af"],
                    help="physical (missing->REF, fair to read-based tools) vs MAR")
    a = ap.parse_args()

    tr = load_truth(a.truth)
    tr["truth_af"] = tr[a.truth_col]   # select the estimand to score against
    print(f"scoring against truth column: {a.truth_col}", file=sys.stderr)
    if a.miss_thr is not None:
        keep = tr["n_called"] >= (1 - a.miss_thr) * a.F
        print(f"miss filter thr={a.miss_thr}: {keep.sum()}/{len(tr)} sites kept", file=sys.stderr)
        tr = tr[keep]

    loaders = {"kMate": (a.kmate, load_kmate),
               "hapFIRE": (a.hapfire, load_3col),
               "vg-giraffe": (a.vg, load_3col)}
    rows, cells = [], []
    for name, (path, fn) in loaders.items():
        if not path:
            continue
        est = fn(path)
        j = tr.merge(est, on="k", how="inner")
        mt = metrics(j["truth_af"].values, j["est"].values)
        mt["tool"] = name
        rows.append(mt)
        cells.append((f"{name}", j["truth_af"].values, j["est"].values))
        print(f"{name}: n={mt['n']:,} MAE={mt['MAE']:.4f} RMSE={mt['RMSE']:.4f} "
              f"R2={mt['R2']:.4f} r={mt['pearson_r']:.4f}", file=sys.stderr)

    res = pd.DataFrame(rows)[["tool", "n", "MAE", "RMSE", "R2", "pearson_r"]]
    res.to_csv(a.out_prefix + "_metrics.tsv", sep="\t", index=False)
    grid_panel(cells, "", a.out_prefix + "_scatter.png", ncols=len(cells), cell=4.0)
    print(res.to_string(index=False))

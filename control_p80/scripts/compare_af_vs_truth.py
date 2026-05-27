#!/usr/bin/env python3
"""
Compare one or more cactus_em AF estimate TSVs against a recomb_truth TSV.

Join on (chrom, pos, ref_len, alt_len). Reports MAE / RMSE / R^2 / Pearson r /
signed bias / outlier-rate (|est-truth|>0.10), overall and stratified by variant
class (SNP / indel / SV>=50bp). Truth rows with NaN truth_af (info->0) are dropped.

Usage:
  compare_af_vs_truth.py --truth recomb_truth.tsv.gz \
      --est consensus=path/to/consensus.tsv raw=path/to/raw.tsv
"""
import argparse
import numpy as np
import pandas as pd

KEYS = ["chrom", "pos", "ref_len", "alt_len"]


def var_class(ref_len, alt_len):
    out = np.where(np.maximum(ref_len, alt_len) >= 50, "SV",
                   np.where((ref_len == 1) & (alt_len == 1), "SNP", "indel"))
    return out


def metrics(est, truth):
    d = est - truth
    n = len(d)
    mae = np.mean(np.abs(d))
    rmse = np.sqrt(np.mean(d ** 2))
    bias = np.mean(d)
    out = np.mean(np.abs(d) > 0.10)
    ss_res = np.sum(d ** 2)
    ss_tot = np.sum((truth - truth.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    r = np.corrcoef(est, truth)[0, 1] if n > 1 and est.std() > 0 and truth.std() > 0 else float("nan")
    return dict(n=n, MAE=mae, RMSE=rmse, R2=r2, r=r, bias=bias, outlier=out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--est", nargs="+", required=True, help="label=path ...")
    args = ap.parse_args()

    truth = pd.read_csv(args.truth, sep="\t")
    truth = truth.dropna(subset=["truth_af"]).drop_duplicates(KEYS)
    truth["vclass"] = var_class(truth["ref_len"].values, truth["alt_len"].values)
    print(f"[truth] {len(truth):,} records with non-NaN truth_af "
          f"(SNP={ (truth.vclass=='SNP').sum():,}, indel={ (truth.vclass=='indel').sum():,}, "
          f"SV={ (truth.vclass=='SV').sum():,})")

    for spec in args.est:
        label, path = spec.split("=", 1)
        est = pd.read_csv(path, sep="\t").drop_duplicates(KEYS)
        m = truth.merge(est[KEYS + ["alt_freq"]], on=KEYS, how="inner")
        print(f"\n=== {label}  ({path}) ===")
        print(f"  joined {len(m):,} / {len(truth):,} truth records")
        header = f"  {'class':6} {'n':>10} {'MAE':>8} {'RMSE':>8} {'R2':>8} {'r':>7} {'bias':>9} {'|d|>.1':>8}"
        print(header)
        for cls in ["ALL", "SNP", "indel", "SV"]:
            sub = m if cls == "ALL" else m[m.vclass == cls]
            if len(sub) == 0:
                continue
            mt = metrics(sub["alt_freq"].values.astype(float),
                         sub["truth_af"].values.astype(float))
            print(f"  {cls:6} {mt['n']:>10,} {mt['MAE']:>8.4f} {mt['RMSE']:>8.4f} "
                  f"{mt['R2']:>8.4f} {mt['r']:>7.4f} {mt['bias']:>9.5f} {mt['outlier']:>8.4f}")


if __name__ == "__main__":
    main()

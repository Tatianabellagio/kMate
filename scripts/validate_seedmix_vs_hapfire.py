#!/usr/bin/env python
"""Validate kMate SEEDMIX AF estimates against xwu's hapFIRE SNP-frequency truth.

For each SEEDMIX sample S{N}, join the kMate per-record TSV to the hapFIRE
per-sample truth and report agreement.

kMate TSV  (header): chrom pos ref_len alt_len alt_freq info n_called se
  - chrom is "Chr1".."Chr5"; we restrict to SNPs (ref_len==1 & alt_len==1).
hapFIRE truth (no header): chrom_int  pos  freq   (3.24M SNPs; no ref/alt)
  - so the join key is (chrom, pos). kMate positions carrying >1 SNP ALT
    (multi-allelic) are ambiguous against a single truth freq and are dropped.

The truth's allele convention is detected from the sign of the correlation:
if Pearson r < 0 we report the polarity-flipped (1 - truth) agreement instead,
which would mean the truth is REF- rather than ALT-frequency.

Usage:
  python scripts/validate_seedmix_vs_hapfire.py \
      [--kmate-dir analysis/grenenet_gea/common/rerun_kfw_hb/seedmix] \
      [--truth-dir <hapFIRE seed_mix dir>] \
      [--samples 1 2 3 4 5 6 7 8] \
      [--min-called 0] [--out analysis/panel_qc/seedmix_validation_vs_hapfire.tsv]
"""
import argparse, os, sys
import numpy as np
import pandas as pd

DEF_TRUTH = ("/global/scratch/projects/fc_moilab/projects/grenenet-phase1/"
             "frequency/hapFIRE_frequencies/seed_mix")


def load_kmate_snps(path, min_called=0):
    df = pd.read_csv(path, sep="\t",
                     usecols=["chrom", "pos", "ref_len", "alt_len",
                              "alt_freq", "n_called"],
                     dtype={"chrom": str, "pos": np.int64,
                            "ref_len": np.int64, "alt_len": np.int64,
                            "alt_freq": np.float64, "n_called": np.int64})
    snp = df[(df.ref_len == 1) & (df.alt_len == 1) & df.alt_freq.notna()].copy()
    if min_called > 0:
        snp = snp[snp.n_called >= min_called]
    # Chr{N} -> int N
    snp["chrom_i"] = snp.chrom.str.replace("Chr", "", regex=False).astype(np.int64)
    # drop multi-allelic positions (ambiguous vs a single truth freq)
    dup = snp.duplicated(subset=["chrom_i", "pos"], keep=False)
    n_multi = int(dup.sum())
    snp = snp[~dup]
    return snp[["chrom_i", "pos", "alt_freq", "n_called"]], n_multi


def load_truth(path):
    t = pd.read_csv(path, sep="\t", header=None,
                    names=["chrom_i", "pos", "truth"],
                    dtype={"chrom_i": np.int64, "pos": np.int64,
                           "truth": np.float64})
    return t


def metrics(a, b):
    """a=kmate, b=truth. Returns dict of agreement metrics."""
    n = len(a)
    if n < 2:
        return dict(n=n, pearson=np.nan, rmse=np.nan, mae=np.nan, bias=np.nan)
    r = float(np.corrcoef(a, b)[0, 1])
    d = a - b
    return dict(n=n, pearson=r,
                rmse=float(np.sqrt(np.mean(d**2))),
                mae=float(np.mean(np.abs(d))),
                bias=float(np.mean(d)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmate-dir", default="analysis/grenenet_gea/common/rerun_kfw_hb/seedmix")
    ap.add_argument("--truth-dir", default=DEF_TRUTH)
    ap.add_argument("--samples", nargs="+", type=int,
                    default=[1, 2, 3, 4, 5, 6, 7, 8])
    ap.add_argument("--min-called", type=int, default=0,
                    help="only keep SNPs with n_called >= this (panel QC filter)")
    ap.add_argument("--out", default="analysis/panel_qc/seedmix_validation_vs_hapfire.tsv")
    args = ap.parse_args()

    rows = []
    for n in args.samples:
        kpath = os.path.join(args.kmate_dir, f"SEEDMIX_S{n}.tsv")
        tpath = os.path.join(args.truth_dir, f"s{n}_snp_frequency.txt")
        if not os.path.exists(kpath):
            print(f"[S{n}] SKIP — kMate TSV missing: {kpath}", file=sys.stderr)
            continue
        if not os.path.exists(tpath):
            print(f"[S{n}] SKIP — truth missing: {tpath}", file=sys.stderr)
            continue

        kmate, n_multi = load_kmate_snps(kpath, args.min_called)
        truth = load_truth(tpath)
        m = kmate.merge(truth, on=["chrom_i", "pos"], how="inner")

        a = m.alt_freq.to_numpy()
        b = m.truth.to_numpy()
        met = metrics(a, b)
        polarity = "ALT"
        if met["pearson"] < 0:               # truth likely REF-frequency
            met_flip = metrics(a, 1.0 - b)
            polarity = "REF(1-truth)"
            met = met_flip

        row = dict(sample=f"S{n}",
                   n_kmate_snp=len(kmate), n_multi_dropped=n_multi,
                   n_joined=met["n"], polarity=polarity,
                   pearson=met["pearson"], rmse=met["rmse"],
                   mae=met["mae"], bias=met["bias"])
        rows.append(row)
        print(f"[S{n}] joined={met['n']:>8,}  r={met['pearson']:.4f}  "
              f"RMSE={met['rmse']:.4f}  MAE={met['mae']:.4f}  "
              f"bias={met['bias']:+.4f}  polarity={polarity}  "
              f"(multi-allelic dropped={n_multi:,})", flush=True)

    if rows:
        out = pd.DataFrame(rows)
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        out.to_csv(args.out, sep="\t", index=False, float_format="%.5f")
        print(f"\nwrote {len(rows)} rows -> {args.out}")
        # summary line
        print(f"MEAN over samples: r={out.pearson.mean():.4f}  "
              f"RMSE={out.rmse.mean():.4f}  MAE={out.mae.mean():.4f}  "
              f"bias={out.bias.mean():+.4f}")


if __name__ == "__main__":
    main()

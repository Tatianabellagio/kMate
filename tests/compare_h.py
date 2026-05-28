#!/usr/bin/env python3
"""
Compare EM h estimates (from cn_full built by different k-mer indexes) against
the g0 ground-truth founder frequencies.

For each (index variant, sim) the h-only sweep wrote <prefix>.h_sweep.npz with
h_per_alpha[0] = the plain-EM h (alpha=0). We align it to h_truth.tsv by founder
and report:
  MAE, RMSE, Pearson r vs truth
  off-target mass  : sum of estimated h on founders whose truth is 0 (subset sims)
  top-K recall     : fraction of the K truth-nonzero founders that are in the
                     estimator's top-K by h (K = n truth-nonzero)

Usage:
  compare_h.py --truth <h_truth.tsv> --label:PG <prefix.h_sweep.npz> --label:HAP <...> ...
  (any number of  NAME=PATH  positional args)
"""
import argparse
import numpy as np
import sys


def load_truth(path):
    f, v = [], []
    with open(path) as fh:
        next(fh)
        for line in fh:
            a, b = line.rstrip("\n").split("\t")
            f.append(a); v.append(float(b))
    return np.array(f), np.array(v)


def load_h(path):
    d = np.load(path, allow_pickle=True)
    founders = np.asarray(d["founders"]).astype(str)
    h = np.asarray(d["h_per_alpha"])[0].astype(float)  # alpha=0 row
    return founders, h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("pairs", nargs="+", help="NAME=path.h_sweep.npz")
    args = ap.parse_args()

    tf, tv = load_truth(args.truth)
    tmap = dict(zip(tf, tv))
    nz = tv > 0
    K = int(nz.sum())
    print(f"truth: {len(tf)} founders, {K} non-zero (each {tv[nz][0]:.5f} if uniform)\n")

    hdr = f"{'index':<10} {'MAE':>9} {'RMSE':>9} {'pearson':>8} {'offtgt_mass':>12} {'topK_recall':>12}"
    print(hdr); print("-" * len(hdr))
    rows = []
    for pair in args.pairs:
        name, path = pair.split("=", 1)
        hf, h = load_h(path)
        # align to truth founder order
        t = np.array([tmap.get(x, np.nan) for x in hf])
        if np.isnan(t).any():
            miss = np.isnan(t).sum()
            print(f"  WARN {name}: {miss} founders not in truth (skipped from metrics)", file=sys.stderr)
        ok = ~np.isnan(t)
        h_ = h[ok]; t_ = t[ok]; hf_ = hf[ok]
        mae = np.mean(np.abs(h_ - t_))
        rmse = np.sqrt(np.mean((h_ - t_) ** 2))
        r = np.corrcoef(h_, t_)[0, 1] if h_.std() > 0 else float("nan")
        nz_ = t_ > 0
        offt = float(h_[~nz_].sum())
        kk = int(nz_.sum())
        topk = set(np.argsort(h_)[::-1][:kk].tolist())
        truth_set = set(np.where(nz_)[0].tolist())
        recall = len(topk & truth_set) / kk if kk else float("nan")
        print(f"{name:<10} {mae:>9.5f} {rmse:>9.5f} {r:>8.4f} {offt:>12.5f} {recall:>12.3f}")
        rows.append((name, mae, rmse, r, offt, recall))

    if len(rows) > 1:
        base = rows[0]
        print(f"\nrelative to {base[0]} (baseline):")
        for nm, mae, rmse, r, offt, rec in rows[1:]:
            print(f"  {nm:<10} dMAE={mae-base[1]:+.5f} ({(mae/base[1]-1)*100:+.1f}%)  "
                  f"dRMSE={rmse-base[2]:+.5f}  d_offtgt={offt-base[4]:+.5f}")


if __name__ == "__main__":
    main()

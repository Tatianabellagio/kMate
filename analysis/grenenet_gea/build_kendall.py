#!/usr/bin/env python
"""Naive climate GEA: per-SV Kendall-tau of Delta-p vs site temperature.

Delta-p = (generation-3 pool AF) - p0 (founding SEEDMIX). For each non-SNP record
we compute Kendall's tau-b between Delta-p across the 193 gen-3 `site_gen_plot`
pools and each pool's site temperature (bio1, the figure's x-axis). Positive tau =
allele rises in warm gardens / falls in cold = candidate climate-adaptive SV.

NAIVE: pools within a site share one climate value (ties; tau-b handles them) and
are non-independent (pseudoreplication NOT corrected here) — a first pass to rank
candidates for the poster, not a structure-corrected test (no LFMM/WZA yet).

Output (--out dir):
  kendall_gen3_bio1.npz   tau, pval (+ record meta chrom/pos/ref_len/alt_len,
                          sv_size, p0, dp_mean) for all non-SNP records
  kendall_gen3_bio1.top.csv   top candidates (true SV >50bp), sorted by tau
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy.stats import kendalltau
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

_X = None      # climate vector (per pool), shared via fork
_DPT = None    # Delta-p transposed [records x pools], shared via fork


def _chunk(rng):
    a, b = rng
    tau = np.full(b - a, np.nan); pv = np.full(b - a, np.nan)
    for i in range(a, b):
        y = _DPT[i]
        if np.ptp(y) > 0:                    # skip invariant records (tau undefined)
            t, p = kendalltau(_X, y)
            tau[i - a] = t; pv[i - a] = p
    return a, tau, pv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", type=int, default=3)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--kind", default="nonsnp")
    ap.add_argument("--pooldir", default=f"{lib.GEA}/pool_matrices")
    ap.add_argument("--out", default=f"{lib.GEA}/gea")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--sv-min-bp", type=int, default=50)
    ap.add_argument("--sv-only", action="store_true",
                    help="compute on true SVs (>sv-min-bp) ONLY (clean SV npz)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    mt = pd.read_csv(f"{args.pooldir}/pool_gen{args.gen}_{args.kind}.meta.csv")
    x = mt[args.climate].to_numpy(float)
    af = np.load(f"{args.pooldir}/pool_gen{args.gen}_{args.kind}_af.npy")  # [pools x rec]
    p0 = np.load(f"{lib.AF_STORE}/p0_{args.kind}.npy")
    idx = dict(np.load(f"{lib.AF_STORE}/index_{args.kind}.npz", allow_pickle=True))
    # --sv-only: restrict the WHOLE computation to true SVs (>sv_min_bp), so the
    # saved arrays + GIF/QQ/Bonferroni are unambiguously the SV set (no small indels).
    if args.sv_only:
        keep = np.abs(idx["alt_len"] - idx["ref_len"]) > args.sv_min_bp
        af = af[:, keep]; p0 = p0[keep]
        for k in ("chrom", "pos", "ref_len", "alt_len"):
            idx[k] = idx[k][keep]
        print(f"--sv-only: restricted to {int(keep.sum()):,} SVs >{args.sv_min_bp}bp", flush=True)
    n_rec = af.shape[1]
    print(f"gen{args.gen} {args.kind}: {af.shape[0]} pools x {n_rec:,} records vs "
          f"{args.climate} (range {x.min():.1f}-{x.max():.1f})", flush=True)

    global _X, _DPT
    _X = x
    _DPT = np.ascontiguousarray((af - p0[None, :]).T)   # [records x pools] float32
    del af

    step = 20000
    ranges = [(a, min(a + step, n_rec)) for a in range(0, n_rec, step)]
    tau = np.empty(n_rec, np.float32); pv = np.empty(n_rec, np.float32)
    with Pool(args.threads) as pool:
        for a, t, p in pool.imap_unordered(_chunk, ranges):
            tau[a:a + len(t)] = t; pv[a:a + len(p)] = p
    print("kendall done", flush=True)

    ref_len = idx["ref_len"]; alt_len = idx["alt_len"]
    sv_size = np.abs(alt_len - ref_len)
    dp_mean = np.nanmean(_DPT, axis=1).astype(np.float32)
    tag = f"_svonly{args.sv_min_bp}" if args.sv_only else ""
    base = f"{args.out}/kendall_gen{args.gen}_{args.climate}{tag}"
    np.savez(f"{base}.npz",
             tau=tau, pval=pv, chrom=idx["chrom"], pos=idx["pos"],
             ref_len=ref_len, alt_len=alt_len, sv_size=sv_size,
             p0=p0, dp_mean=dp_mean)
    # top candidates: true SVs, finite tau, ranked by tau (most up-in-warm first)
    sv = sv_size > args.sv_min_bp
    df = pd.DataFrame(dict(
        chrom=idx["chrom"][sv], pos=idx["pos"][sv], ref_len=ref_len[sv],
        alt_len=alt_len[sv], sv_size=sv_size[sv], tau=tau[sv], pval=pv[sv],
        p0=p0[sv], dp_mean=dp_mean[sv]))
    df["rec_index"] = np.where(sv)[0]
    df = df[np.isfinite(df.tau)].sort_values("tau", ascending=False)
    df.head(500).to_csv(f"{base}.top.csv", index=False)
    print(f"SVs(>{args.sv_min_bp}bp) tested: {int(sv.sum()):,} | "
          f"tau>0 & p<0.05: {int(((df.tau>0)&(df.pval<0.05)).sum()):,} | "
          f"top tau={df.tau.iloc[0]:.3f} (p={df.pval.iloc[0]:.1e})")
    print(f"-> {base}.npz + .top.csv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Is the SV > SNP per-variant |s| excess REAL, or a k-mer-support (noise) artifact?

SVs are estimated from fewer k-mers -> higher AF estimation noise -> a noisier AF
trajectory makes a bigger spurious logit-slope |s|, even with no real selection. This
script asks whether the ~1.07x SV/SNP median-|s| excess survives once we also control
for per-variant SUPPORT (kMate n_called; higher = more data = lower AF SE), not just p0.

Three readouts (site 4 by default):
  1. do SVs have lower support than SNPs (confirm the premise)?
  2. p0- AND support-matched median-|s| SV/SNP ratio: if it collapses to ~1 -> the excess
     was noise; if it persists -> real (SVs move more even at matched precision).
  3. noise-floor model: predicted pure-noise |s| ~ sqrt(p0(1-p0)/n_eff) * c; compare the
     SIGNAL |s|-minus-floor across classes. n_eff from support.

INPUT: af_store (p0_{snp,nonsnp}, nc_{snp,nonsnp}/<sample>.npy), pool_matrices (via site_scoef).
OUTPUT: site_temporal/site{S}_noise_check.{json,png}
  PY=<plotting env for numbers; render fig in basic>; $PY noise_check_sv_s.py --site 4
"""
from __future__ import annotations
import argparse, glob, json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from site_variant_temporal_scoef import site_scoef, STORE, SV_MIN_BP


def mean_support(kind, n_samples):
    """Per-variant mean n_called over n_samples cohort files (support ~ AF precision)."""
    fs = sorted(glob.glob(f"{STORE}/nc_{kind}/*.npy"))[:n_samples]
    acc = None
    for f in fs:
        a = np.load(f).astype(np.float32)
        acc = a if acc is None else acc + a
    return acc / len(fs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=int, default=4)
    ap.add_argument("--out", default=f"{lib.GEA}/site_temporal")
    ap.add_argument("--n-support", type=int, default=30, help="#samples to average n_called")
    ap.add_argument("--min-p0", type=float, default=0.02)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    idx_non = np.load(f"{STORE}/index_nonsnp.npz")
    dlen = np.abs(idx_non["alt_len"].astype(np.int64) - idx_non["ref_len"].astype(np.int64))
    Rsnp = site_scoef(args.site, "snp"); Rnon = site_scoef(args.site, "nonsnp")
    sup_snp = mean_support("snp", args.n_support)
    sup_non = mean_support("nonsnp", args.n_support)

    def frame(name, R, sup, sel):
        p0 = R["p0"][sel]
        keep = (p0 >= args.min_p0) & (p0 <= 1 - args.min_p0) & np.isfinite(R["s"][sel])
        return pd.DataFrame(dict(cls=name, p0=p0[keep], abs_s=np.abs(R["s"][sel])[keep],
                                 support=sup[sel][keep]))
    V = pd.concat([
        frame("snp", Rsnp, sup_snp, np.ones(Rsnp["s"].shape[0], bool)),
        frame("smallindel", Rnon, sup_non, dlen <= SV_MIN_BP),
        frame("sv", Rnon, sup_non, dlen > SV_MIN_BP)], ignore_index=True)

    # 1. support by class
    sup_med = V.groupby("cls").support.median()
    print("[1] median n_called:  " + "  ".join(f"{c}={sup_med[c]:.0f}" for c in
          ("snp", "smallindel", "sv")))

    # 2. p0- AND support-matched SV/SNP |s| ratio (2D deciles)
    V["p0b"] = pd.qcut(V.p0, 10, labels=False, duplicates="drop")
    V["supb"] = pd.qcut(V.support, 10, labels=False, duplicates="drop")
    med = V.groupby(["p0b", "supb", "cls"]).abs_s.median().unstack()
    cnt = V[V.cls == "sv"].groupby(["p0b", "supb"]).size()
    ratio_cells = (med["sv"] / med["snp"]).dropna()
    w = cnt.reindex(ratio_cells.index).fillna(0)
    matched2d = float(np.nansum(ratio_cells * w) / w.sum())
    p0only = float(_p0_matched(V, "sv", "snp"))
    raw = float(V[V.cls == "sv"].abs_s.median() / V[V.cls == "snp"].abs_s.median())
    print(f"[2] |s| SV/SNP ratio:  raw={raw:.3f}  p0-matched={p0only:.3f}  "
          f"p0+support-matched={matched2d:.3f}")

    # 3. noise-floor model: pure-noise |s| ~ sqrt(p0(1-p0)/support); compare signal
    V["floor"] = np.sqrt(V.p0 * (1 - V.p0) / np.maximum(V.support, 1.0))
    # scale floor to a slope: half-normal mean over 4 points, const cancels in class ratio
    fl = V.groupby("cls").floor.median()
    print(f"[3] noise floor proxy (sqrt(p0(1-p0)/support)) median: " +
          "  ".join(f"{c}={fl[c]:.4f}" for c in ("snp", "smallindel", "sv")))
    print(f"    -> SV/SNP floor ratio={fl['sv']/fl['snp']:.3f}  vs observed |s| ratio "
          f"raw={raw:.3f} (if similar, excess is noise)")

    summ = dict(site=args.site, median_support={c: float(sup_med[c]) for c in
                ("snp", "smallindel", "sv")},
                ratio_raw=raw, ratio_p0matched=p0only, ratio_p0_support_matched=matched2d,
                noise_floor_ratio_sv_snp=float(fl["sv"] / fl["snp"]))
    json.dump(summ, open(f"{args.out}/site{args.site}_noise_check.json", "w"), indent=2)

    # plot data (render in basic env)
    sb = V.groupby(["supb", "cls"]).abs_s.median().unstack()
    supx = V.groupby("supb").support.median()
    np.savez(f"{args.out}/site{args.site}_noise_check_plotdata.npz",
             supx=supx.values, s_snp=sb["snp"].values, s_sv=sb["sv"].values,
             s_smallindel=sb["smallindel"].values,
             ratio_raw=raw, ratio_p0=p0only, ratio_2d=matched2d,
             floor_ratio=float(fl["sv"] / fl["snp"]))
    print(f"-> {args.out}/site{args.site}_noise_check.json + _plotdata.npz")


def _p0_matched(V, a, b):
    med = V.groupby(["p0b", "cls"]).abs_s.median().unstack()
    cnt = V[V.cls == a].groupby("p0b").size()
    r = (med[a] / med[b]).reindex(cnt.index)
    return np.nansum(r * cnt) / cnt[np.isfinite(r)].sum()


if __name__ == "__main__":
    main()

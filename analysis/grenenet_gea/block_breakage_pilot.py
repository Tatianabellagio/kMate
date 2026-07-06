#!/usr/bin/env python
"""Per-sample block 'breakage' vs climate (PILOT, Chr1).

kMate models each pool as a mixture of the 231 founder haplotypes (freqs h per chrom)
and the .tsv alt_freq is the OBSERVED k-mer frequency. An intact founder mixture would
reproduce the h-projection exactly, so |observed - (h·V_pa)| is the local mosaic/breakage
signal. Per block we take the mean residual; a block is 'broken' if it exceeds THR.
Per sample we report % blocks broken + mean residual, joined to bio1 + coverage.

Output: results/grenenet_gea/blocks_mcf90/breakage_pilot.csv  (+ caller plots it)
Usage: block_breakage_pilot.py --n 30 [--thr 0.08] [--min-cov 3] [--samples a,b,c]
"""
import os, sys, glob, argparse
import numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

OUT = "results/grenenet_kmate_arch3"
BLOCKS = "results/grenenet_gea/blocks_mcf90/chr1_clq0.9_blocks_clq0.9.tsv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--thr", type=float, default=0.08)
    ap.add_argument("--min-cov", type=float, default=3.0)
    ap.add_argument("--samples", default=None, help="comma list to override selection")
    a = ap.parse_args()

    print("loading panel V_pa (Chr1) + blocks ...", flush=True)
    vp = sp.load_npz("panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz").astype(np.float64)  # 231 x R
    pos = np.load("panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz", allow_pickle=True)["pos"].astype(np.int64)
    blocks = pd.read_csv(BLOCKS, sep="\t")
    bstart = blocks.start_pos.values; bend = blocks.end_pos.values

    avail = sorted({os.path.basename(p).replace("_Chr1.tsv", "")
                    for p in glob.glob(f"{OUT}/*_Chr1.tsv")})
    meta = lib.cohort_meta(avail).dropna(subset=["bio1"])
    meta = meta[meta.coverage >= a.min_cov]
    # latest generation per plot (most evolved), then even spread across bio1
    meta = meta.sort_values("generation").groupby(["site", "plot"]).tail(1).sort_values("bio1")
    if a.samples:
        pick = meta.loc[meta.index.intersection(a.samples.split(","))]
    else:
        idx = np.linspace(0, len(meta) - 1, min(a.n, len(meta))).astype(int)
        pick = meta.iloc[np.unique(idx)]
    print(f"{len(avail)} samples available; pilot picks {len(pick)} across bio1 "
          f"[{pick.bio1.min():.1f}, {pick.bio1.max():.1f}] C", flush=True)

    rows = []
    for sid, mr in pick.iterrows():
        try:
            h = np.load(f"{OUT}/{sid}_Chr1.h_per_chrom.npz", allow_pickle=True)["Chr1"]
            tsv = pd.read_csv(f"{OUT}/{sid}_Chr1.tsv", sep="\t", usecols=["pos", "alt_freq"])
        except Exception as e:
            print(f"  skip {sid}: {e}"); continue
        proj = np.asarray(h @ vp).ravel()                       # expected alt freq per record
        pj = pd.DataFrame({"pos": pos, "proj": proj}).groupby("pos", as_index=False).proj.mean()
        ob = tsv.groupby("pos", as_index=False).alt_freq.mean()
        d = pj.merge(ob, on="pos", how="inner")
        d["resid"] = (d.alt_freq - d.proj).abs()
        P = d.pos.values; R = d.resid.values; o = np.argsort(P); P = P[o]; R = R[o]
        # mean residual per block
        broke, tot, resids = 0, 0, []
        for s, e in zip(bstart, bend):
            lo = np.searchsorted(P, s); hi = np.searchsorted(P, e, side="right")
            if hi - lo >= 2:
                mr_blk = R[lo:hi].mean(); resids.append(mr_blk)
                tot += 1; broke += mr_blk > a.thr
        rows.append(dict(sample=sid, site=int(mr.site), generation=int(mr.generation),
                         bio1=float(mr.bio1), coverage=float(mr.coverage),
                         n_blocks=tot, pct_broken=100*broke/max(tot, 1),
                         mean_resid=float(np.mean(resids))))
        print(f"  {sid} site{int(mr.site)} bio1={mr.bio1:.1f} cov={mr.coverage:.1f} "
              f"-> {100*broke/max(tot,1):.1f}% broken, meanresid {np.mean(resids):.4f}", flush=True)

    df = pd.DataFrame(rows)
    out = "results/grenenet_gea/blocks_mcf90/breakage_pilot.csv"
    df.to_csv(out, index=False)
    print(f"\nwrote {len(df)} samples -> {out}")
    if len(df) > 3:
        print(f"corr(bio1, %broken)     = {df.bio1.corr(df.pct_broken):.3f}")
        print(f"corr(coverage, %broken) = {df.coverage.corr(df.pct_broken):.3f}  (confound check)")
        print(f"corr(bio1, mean_resid)  = {df.bio1.corr(df.mean_resid):.3f}")


if __name__ == "__main__":
    main()

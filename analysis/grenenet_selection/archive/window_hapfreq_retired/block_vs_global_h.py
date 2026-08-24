#!/usr/bin/env python
"""Recombination diagnostic: is block-local founder h the SAME as global (chrom-wide) h?

For each replicate sample at a site, every window npz stores both the locally-estimated
per-block founder frequencies (Chr*_h_blocks, U x 231) and the chromosome-wide global
founder frequencies (Chr*_global_h, 231) for the SAME sample. If recombination is
negligible (whole ecotype genomes segregate as units), each block's h == global h ->
r ~ 1 everywhere. Recombination decouples a block's local ancestry from the genome-wide
ecotype composition -> lower r.

Per status-0 (genuinely local) block we compute r(h_block, global_h) and L1 distance
across the 231 founders, for each gen-3 site-4 replicate. Reports the distribution and a
genome map of the median-over-replicates similarity, plus whether low-r blocks are
CONSISTENT across replicates (shared recombination/structure) vs idiosyncratic (noise).

Env: kmate. SITE / GEN via env. Writes site<ID>_block_vs_global_h.png + csv.
"""
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

WIN = "results/grenenet_kmate_window"
GEA = "results/grenenet_gea/hapfreq"
SITE = int(os.environ.get("SITE", 4))
GEN = int(os.environ.get("GEN", 3))
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]


def corr_rows(A, g):
    """row-wise Pearson r of each row of A (m x k) against vector g (k,)."""
    A = A - A.mean(1, keepdims=True); g = g - g.mean()
    num = (A * g).sum(1)
    den = np.sqrt((A**2).sum(1) * (g**2).sum()) + 1e-12
    return num / den


def load_sample(samp):
    """Stack all chroms: per-block (chrom,start,end,status,r,L1) for one sample."""
    rows = []
    for ch in CHROMS:
        f = f"{WIN}/{samp}_{ch}.h_blocks_per_chrom.npz"
        if not os.path.exists(f):
            return None
        z = np.load(f, allow_pickle=True)
        hb = z[f"{ch}_h_blocks"].astype(np.float64)
        gh = z[f"{ch}_global_h"].astype(np.float64)
        st = z[f"{ch}_status"]
        r = corr_rows(hb, gh)
        l1 = np.abs(hb - gh[None, :]).sum(1)
        rows.append(pd.DataFrame({
            "chrom": z[f"{ch}_block_chrom"].astype(str),
            "start": z[f"{ch}_block_start"], "end": z[f"{ch}_block_end"],
            "status": st, "r": r, "l1": l1}))
    return pd.concat(rows, ignore_index=True)


def main():
    pt = lib.pool_table()
    samp = pt[(pt.site == SITE) & (pt.generation == GEN)].sampleid.astype(str).tolist()
    samp = [s for s in samp if os.path.exists(f"{WIN}/{s}_Chr1.h_blocks_per_chrom.npz")]
    print(f"site {SITE} gen {GEN}: {len(samp)} replicate samples")

    pieces = []
    for s in samp:
        d = load_sample(s)
        d["samp"] = s
        pieces.append(d)
    A = pd.concat(pieces, ignore_index=True)
    loc = A[A.status == 0].copy()                 # genuinely local blocks only
    key = loc.chrom + ":" + loc.start.astype(str)

    print(f"\nstatus-0 (locally estimated) block-replicate pairs: {len(loc):,}")
    print(f"  r(block_h, global_h):  median {loc.r.median():.3f}  "
          f"mean {loc.r.mean():.3f}  | r>0.95: {100*(loc.r>0.95).mean():.0f}%  "
          f"r>0.90: {100*(loc.r>0.90).mean():.0f}%  r<0.80: {100*(loc.r<0.80).mean():.0f}%")
    print(f"  L1 |block_h - global_h|: median {loc.l1.median():.3f} "
          f"(0=identical, 2=disjoint)")

    # per-block median across replicates + consistency (how many reps decoupled)
    loc["key"] = key
    g = loc.groupby("key").agg(chrom=("chrom", "first"), start=("start", "first"),
                               end=("end", "first"), r_med=("r", "median"),
                               r_iqr=("r", lambda x: x.quantile(.75) - x.quantile(.25)),
                               n_rep=("r", "size"),
                               frac_dec=("r", lambda x: (x < 0.8).mean())).reset_index()
    print(f"\nper-block (median over {len(samp)} reps): "
          f"r_med>0.95 {100*(g.r_med>0.95).mean():.0f}%  r_med<0.80 {100*(g.r_med<0.80).mean():.0f}%")
    print(f"  of blocks decoupled (r_med<0.8): consistently (frac_dec>0.7 of reps) "
          f"{int(((g.r_med<0.8)&(g.frac_dec>0.7)).sum())} vs idiosyncratic "
          f"{int(((g.r_med<0.8)&(g.frac_dec<=0.7)).sum())}")
    g.sort_values("r_med").to_csv(f"{GEA}/site{SITE}_block_vs_global_h.csv", index=False)

    # genome-cumulative x for the map
    off, centers, x = 0.0, [], np.zeros(len(g))
    g = g.sort_values(["chrom", "start"]).reset_index(drop=True)
    mid = (g.start + g.end) / 2
    for ch in CHROMS:
        m = (g.chrom == ch).to_numpy()
        if not m.any():
            continue
        x[m] = mid[m] + off; centers.append(off + mid[m].max() / 2); off += mid[m].max() * 1.02
    g["x"] = x

    fig, (a0, a1) = plt.subplots(2, 1, figsize=(11, 7))
    a0.hist(loc.r, bins=120, color="#34495e")
    a0.axvline(loc.r.median(), color="firebrick", ls="--",
               label=f"median r = {loc.r.median():.3f}")
    a0.set_xlabel("r( block-local h , global h )  per block-replicate")
    a0.set_ylabel("count"); a0.set_title(
        f"Site {SITE} gen {GEN}: is block-local founder h the same as genome-wide h? "
        f"({len(samp)} replicates, {len(loc):,} local blocks)", loc="left", fontsize=11)
    a0.legend(frameon=False); a0.spines[["top", "right"]].set_visible(False)

    for i, ch in enumerate(CHROMS):
        m = (g.chrom == ch).to_numpy()
        a1.scatter(g.x[m], g.r_med[m], s=6, c=["#3b4cc0", "#7aa0c4"][i % 2],
                   alpha=.5, edgecolors="none", rasterized=True)
    a1.axhline(0.8, color="firebrick", lw=.8, ls=":", label="decoupled (r<0.8)")
    a1.set_xticks(centers); a1.set_xticklabels(CHROMS)
    a1.set_ylabel("median r over replicates"); a1.set_xlabel("genome position")
    a1.set_title("Where does block-local ancestry decouple from the whole-genome ecotype "
                 "composition? (low = recombination-broken)", loc="left", fontsize=10)
    a1.legend(frameon=False, fontsize=8); a1.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = f"{GEA}/site{SITE}_block_vs_global_h.png"
    fig.savefig(out, dpi=150)
    print(f"\n[done] {out}")


if __name__ == "__main__":
    main()

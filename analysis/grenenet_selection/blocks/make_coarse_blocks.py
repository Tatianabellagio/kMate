#!/usr/bin/env python
"""Coarsen the LD-block map by greedily merging ADJACENT blocks (within a chrom) until
each merged window clears a minimum-variant floor. Produces the candidate window maps for
the 'how coarse must we go for kMate to estimate h' benchmark.

NOTE: floor is on n_variants (the proxy we have); the benchmark agent should ideally also
sweep a min-OBSERVED-KMER floor (they have the per-block k-mer counts). Variant-floor maps
are the starting sweep.

Input  : analysis/grenenet_selection/blocks/results/blocks_mcf90/chr{N}_clq0.9_blocks_clq0.9.tsv
Output : analysis/grenenet_selection/blocks/results/blocks_mcf90/coarse/chr{N}_floor{F}.tsv  (chrom start end n_variants n_merged)
Usage  : make_coarse_blocks.py --floors 8,15,25,40
"""
import os, argparse
import numpy as np, pandas as pd

BR = "analysis/grenenet_selection/blocks/results/blocks_mcf90"
OUTD = f"{BR}/coarse"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]


def coarsen(b, floor):
    """b sorted by start_pos -> merged intervals each with >=floor variants (last may be short)."""
    out = []
    s = e = None; acc = nm = 0
    for _, r in b.iterrows():
        if s is None:
            s = r.start_pos
        e = r.end_pos; acc += r.n_variants; nm += 1
        if acc >= floor:
            out.append((s, e, acc, nm)); s = None; acc = nm = 0
    if s is not None:                      # trailing remainder -> fold into previous window
        if out:
            ps, pe, pa, pn = out[-1]; out[-1] = (ps, e, pa + acc, pn + nm)
        else:
            out.append((s, e, acc, nm))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--floors", default="8,15,25,40")
    a = ap.parse_args()
    floors = [int(x) for x in a.floors.split(",")]
    os.makedirs(OUTD, exist_ok=True)

    print(f"{'floor':>6}{'#windows':>10}{'medVar':>8}{'q25':>6}{'%>=floor':>10}{'med span(bp)':>13}")
    base_total = 0
    for F in floors:
        rows = []
        for ch in CHROMS:
            cl = ch.lower()
            b = pd.read_csv(f"{BR}/{cl}_clq0.9_blocks_clq0.9.tsv", sep="\t").sort_values("start_pos")
            for s, e, nv, nm in coarsen(b, F):
                rows.append((ch, int(s), int(e), int(nv), int(nm)))
        df = pd.DataFrame(rows, columns=["chrom", "start_pos", "end_pos", "n_variants", "n_merged"])
        df.to_csv(f"{OUTD}/floor{F}.tsv", sep="\t", index=False)
        # also per-chrom for array jobs
        for ch in CHROMS:
            df[df.chrom == ch].to_csv(f"{OUTD}/{ch.lower()}_floor{F}.tsv", sep="\t", index=False)
        nv = df.n_variants.values; sp = (df.end_pos - df.start_pos).values
        print(f"{F:>6}{len(df):>10}{int(np.median(nv)):>8}{int(np.percentile(nv,25)):>6}"
              f"{100*(nv>=F).mean():>9.0f}%{int(np.median(sp)):>13}")
    # reference: base map
    nbase = sum(len(pd.read_csv(f"{BR}/{ch.lower()}_clq0.9_blocks_clq0.9.tsv", sep="\t")) for ch in CHROMS)
    print(f"\n(base CLQ0.9 map = {nbase:,} blocks, median 7 var)")
    print(f"coarse maps -> {OUTD}/floor{{F}}.tsv  and per-chrom {OUTD}/chr{{N}}_floor{{F}}.tsv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Gap-free ("tiling") block assignment, following hapFIRE's own rule.

WHY THIS EXISTS
---------------
BigLD defines LD blocks on a FILTERED common variant set (`recompute_blocks.py`:
MAF>0.05, call-rate>=90%, one variant per position). On Chr1 that is 291,133
variants, whereas kMate has 699,626 records. The resulting block TSVs are
therefore intervals ("LD islands") that do not tile the genome.

`lib.assign_clq_blocks` assigns by STRICT INTERVAL CONTAINMENT, so every record
falling between two islands is dropped. Measured on the current data that is:

    snp 39.7%   smallindel 40.3%   nonsnp 40.7%   sv 49.3%   of all records

...silently discarded before WZA ever sees them. Those records are NOT in
distant, unlinked territory -- the median dropped record is 731 bp from the
nearest block edge and 57% are within 1 kb (median inter-block gap: 267 bp).
They are dropped only because they were not in the set used to DEFINE the
blocks. Since the per-variant GEA models (kendall / lfmm / quasi-binomial) run
on every record and WZA is only an aggregation step to control inflation, there
is no methodological reason to drop them.

HOW HapFM DOES IT (the rule implemented here)
---------------------------------------------
HapFM -- the software that actually builds the haploblocks -- runs BigLD on the
common set and then converts the breakpoints to GENOME-WIDE breakpoints that
cover every variant. See
`/global/scratch/users/tbellg/HapFM/bin/utility_functions.py
 ::convert_fine_genomewide_breakpoints`:

    if i == 0:                 left = 0                                  # absorb the head
    elif i == last:            left = common_index[prev_right] + 1; right = r-1   # absorb the tail
    else:                      left = common_index[prev_right] + 1       # absorb the gap

i.e. each block is EXTENDED to begin immediately after the previous block ends,
so inter-island variants are absorbed into the following block rather than
discarded. HapFM then merges any block spanning < 2 variants into a neighbour.
(hapFIRE carries a copy of this same function -- shared code lineage, same
author -- differing only in using a merge floor of 6.)

Note HapFM works in VARIANT-INDEX space (`left = prev_right + 1`); this module
implements the identical partition in POSITION space, where block i owns the
half-open interval (end_{i-1}, end_i].

This module reproduces that rule in position space: block i owns the half-open
interval (end_{i-1}, end_i], the first block owns everything up to end_0, and
the last block owns everything after end_{n-2}. Result: 100% of records keep a
block, while the r2-derived LD boundaries are preserved exactly.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
HAPFM_MIN_VARIANTS = 2        # HapFM's own small-block merge floor (span < 2 variants)


def load_blocks(r2: float, chrom: str) -> pd.DataFrame:
    ci = int(chrom.replace("Chr", ""))
    tag = f"clq{r2}"
    f = f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv"
    return pd.read_csv(f, sep="\t").sort_values("start_pos").reset_index(drop=True)


def assign_tiling(chrom: np.ndarray, pos: np.ndarray, r2: float = 0.9) -> np.ndarray:
    """Assign EVERY (chrom,pos) to a block id, hapFIRE-style (no gaps, no drops).

    Block i owns (end_{i-1}, end_i]; first block absorbs the head, last the tail.
    Block ids keep the same `Chr{n}_{idx}` convention as lib.assign_clq_blocks, so
    ids are directly comparable to the strict-containment assignment.
    """
    chrom = np.asarray(chrom, dtype=str)
    pos = np.asarray(pos, dtype=np.int64)
    out = np.empty(len(pos), dtype=object)
    out[:] = ""
    for ch in np.unique(chrom):
        if ch not in CHROMS:
            continue
        b = load_blocks(r2, ch)
        ends = b["end_pos"].to_numpy(np.int64)
        m = chrom == ch
        # searchsorted on block END positions: the first block whose end >= pos owns it.
        idx = np.searchsorted(ends, pos[m], side="left")
        idx = np.clip(idx, 0, len(ends) - 1)          # tail variants -> last block
        out[m] = [f"{ch}_{i}" for i in idx]
    return out


def merge_small_blocks(block_ids: np.ndarray, min_variants: int = HAPFM_MIN_VARIANTS) -> np.ndarray:
    """Merge blocks holding < min_variants records into the previous block on the
    same chromosome (HapFM merges sub-2-variant windows into a neighbour this way)."""
    block_ids = np.asarray(block_ids, dtype=object)
    counts = pd.Series(block_ids).value_counts()
    remap = {}
    for ch in CHROMS:
        ids = sorted([b for b in counts.index if isinstance(b, str) and b.startswith(ch + "_")],
                     key=lambda s: int(s.split("_")[1]))
        prev_keep = None
        for bid in ids:
            n = counts[bid] + sum(counts[k] for k, v in remap.items() if v == bid)
            if counts[bid] < min_variants and prev_keep is not None:
                remap[bid] = prev_keep
            else:
                prev_keep = remap.get(bid, bid)
    if not remap:
        return block_ids
    return np.array([remap.get(b, b) for b in block_ids], dtype=object)


def _report(r2):
    MA = f"{lib.GEA}/phase1_replication/results/multiaxis"
    print(f"\n===== clq{r2} =====")
    print(f"{'class':11s} {'records':>10s} | {'STRICT kept':>12s} {'%drop':>7s} | "
          f"{'TILING kept':>12s} {'%drop':>7s} | {'blocks':>8s} {'med/blk':>8s}")
    print("-" * 92)
    for cls in ["snp", "sv", "smallindel", "nonsnp"]:
        f = f"{MA}/kendall/kendall_{cls}_gen9_bio1.csv"
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f, usecols=["chrom", "pos"])
        c = d["chrom"].to_numpy(str); p = d["pos"].to_numpy(np.int64)
        strict = lib.assign_clq_blocks(c, p, r2=r2)
        tiled = assign_tiling(c, p, r2=r2)
        tiled = merge_small_blocks(tiled)
        n = len(d); ks = (strict != "").sum(); kt = (tiled != "").sum()
        per = pd.Series(tiled[tiled != ""]).value_counts()
        print(f"{cls:11s} {n:10,} | {ks:12,} {100*(n-ks)/n:6.1f}% | "
              f"{kt:12,} {100*(n-kt)/n:6.1f}% | {len(per):8,} {per.median():8.0f}")


if __name__ == "__main__":
    for r2 in (0.9, 0.5):
        _report(r2)

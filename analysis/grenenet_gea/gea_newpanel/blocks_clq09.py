#!/usr/bin/env python
"""clq0.9 haploblock assignment for the new-panel GEA.

Replaces the phase-1 hapFIRE LD blocks (lib.assign_ld_blocks) with the recomputed
**BigLD CLQ-cut 0.9** haploblocks (`blocks_recompute/chr{N}_clq0.9_blocks_clq0.9.tsv`),
which are much finer (~82k blocks genome-wide vs 16.7k phase-1) and are the project's
trusted haplotype-scale unit.

Block file format: chrom, start_pos, end_pos, n_variants  (one row per block).
Block id = f"{chrom}_{rank}" with rank = 1-based order along the chromosome.

`assign_clq09_blocks(chrom, pos)` -> string array of block ids. A record inside a
block interval gets that block; a record in a between-block gap inherits the NEAREST
block by boundary distance (same spirit as the phase-1 nearest-SNP rule). This keeps
every record assignable so the WZA windows lose no loci.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

BLOCK_DIR = "/global/scratch/users/tbellg/kmate/results/grenenet_gea/blocks_recompute"
_CHROMS = [f"Chr{i}" for i in range(1, 6)]
_CACHE: dict | None = None


def load_blocks() -> pd.DataFrame:
    """All clq0.9 blocks (chrom, start_pos, end_pos, block, mid) sorted per chrom."""
    frames = []
    for c in _CHROMS:
        n = c.replace("Chr", "")
        b = pd.read_csv(f"{BLOCK_DIR}/chr{n}_clq0.9_blocks_clq0.9.tsv", sep="\t")
        b = b.sort_values("start_pos").reset_index(drop=True)
        b["block"] = c + "_" + (b.index + 1).astype(str)
        b["mid"] = (b.start_pos + b.end_pos) / 2.0
        frames.append(b)
    return pd.concat(frames, ignore_index=True)


def assign_clq09_blocks(chrom, pos) -> np.ndarray:
    """Map each (chrom, pos) to its clq0.9 block id (nearest block if in a gap)."""
    global _CACHE
    if _CACHE is None:
        _CACHE = {c: g for c, g in load_blocks().groupby("chrom")}
    chrom = np.asarray(chrom, dtype=object)
    pos = np.asarray(pos, dtype=np.int64)
    out = np.full(len(pos), "", dtype=object)
    for c, g in _CACHE.items():
        m = chrom == c
        if not m.any():
            continue
        st = g.start_pos.to_numpy(); en = g.end_pos.to_numpy(); bid = g.block.to_numpy()
        p = pos[m]
        j = np.searchsorted(st, p, side="right") - 1        # last block with start<=pos
        j = np.clip(j, 0, len(st) - 1)
        contained = (p >= st[j]) & (p <= en[j])
        # for gaps: nearer of block j (left) and j+1 (right) by boundary distance
        jn = np.clip(j + 1, 0, len(st) - 1)
        d_left = np.abs(p - en[j])
        d_right = np.abs(st[jn] - p)
        pick = np.where(contained, j, np.where(d_right < d_left, jn, j))
        out[np.where(m)[0]] = bid[pick]
    return out


if __name__ == "__main__":
    # smoke test / stats
    b = load_blocks()
    print(f"clq0.9 blocks: {len(b):,} total")
    print(b.groupby("chrom").size().to_string())
    span = (b.end_pos - b.start_pos)
    print(f"block span bp: median={span.median():.0f} mean={span.mean():.0f} "
          f"p90={span.quantile(0.9):.0f}")
    print(f"variants/block: median={b.n_variants.median():.0f} max={b.n_variants.max()}")

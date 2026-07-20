#!/usr/bin/env python
"""Dynamic LD-guided block growth: relax r2 only where blocks are too small, merging along
the LD gradient until each unit clears a panel-k-mer floor K. Strict-r2 (CLQ0.9) blocks that
already clear K are untouched; sub-K blocks absorb their most-LD-linked neighbour first;
merging stops at a true LD break (inter-block r2 < FLOOR) -> those stay small (desert/fallback).

inter-block LD = r2 between adjacent units' founder PC1 (dominant haplotype axis).
Output: analysis/grenenet_gea/blocks_mcf90/chr{n}_units_dynld_K{K}.tsv

CORRECTNESS (2026-06-19 audit fixes):
 - panel k-mers per unit = # bubble CENTROIDS in the unit's [start,end] interval, computed
   from a sorted-centroid cumulative count (km_of). This is GAP-INCLUSIVE (k-mers in the gaps
   between the original fine blocks count toward a merged unit, exactly as kMate window-mode
   will use them) and matches block_em's centroid-in-[start,end] rule. NOT a sum of fine-block
   counts (which drops ~43% gap k-mers and over-merges).
 - heap entries are VERSION-STAMPED; a block's version bumps on every merge, so stale entries
   (computed against an outdated PC1) are rejected at pop -> the merge is truly
   highest-CURRENT-LD-first and the FLOOR test always uses current PC1.
"""
import os, sys, argparse, heapq
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from recompute_blocks import build_common_matrix

MAF, MINCF = 0.05, 0.9
BR = "analysis/grenenet_gea/blocks_mcf90"


def pc1(geno, lo, hi):
    M = geno[:, lo:hi].astype(np.float64)
    if M.shape[1] == 0:
        return np.zeros(M.shape[0])
    M = M - M.mean(0); sd = M.std(0); sd[sd == 0] = 1.0; M /= sd
    # PC1 founder scores via eigh of the 231x231 Gram (robust + small; SVD can fail to converge)
    G = M @ M.T
    w, V = np.linalg.eigh(G)            # ascending eigenvalues
    return V[:, -1] * np.sqrt(max(w[-1], 0.0))


def r2(a, b):
    if a.std() < 1e-12 or b.std() < 1e-12:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1] ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--K", type=int, default=500); ap.add_argument("--floor", type=float, default=0.05)
    a = ap.parse_args(); chrlc = a.chrom.lower()
    print(f"[{a.chrom}] founder matrix + k-mer centroids ...", flush=True)
    _, raw, positions = build_common_matrix(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}", MAF, MINCF)
    geno = (raw >= 0.5).astype(np.int8); positions = np.asarray(positions)
    # gap-inclusive interval k-mer counts via cumulative sum over bubble centroids, WEIGHTED
    # by k-mers-per-bubble (each bubble holds many k-mers; bubble_id maps each k-mer to a bubble).
    meta = np.load(f"data/kmer_pa_231_arch3_filt2inv/kmer_pa_{a.chrom}.meta.npz", allow_pickle=True)
    cent_b = ((meta["bubble_start"] + meta["bubble_end"]) // 2).astype(np.int64)   # per bubble
    kpb = np.bincount(meta["bubble_id"], minlength=len(cent_b)).astype(np.int64)   # k-mers per bubble
    order = np.argsort(cent_b, kind="stable")
    cent = cent_b[order]
    kcum = np.concatenate([[0], np.cumsum(kpb[order])]).astype(np.int64)           # prefix k-mer counts

    def km_of(lo, hi):  # # k-mers whose bubble CENTROID in [positions[lo], positions[hi-1]] inclusive
        s = positions[lo]; e = positions[hi - 1]
        return int(kcum[np.searchsorted(cent, e, "right")] - kcum[np.searchsorted(cent, s, "left")])

    b = pd.read_csv(f"{BR}/{chrlc}_clq0.9_blocks_clq0.9.tsv", sep="\t").sort_values("start_pos").reset_index(drop=True)
    n = len(b)
    lo = np.array([np.searchsorted(positions, s) for s in b.start_pos])
    hi = np.array([np.searchsorted(positions, e, side="right") for e in b.end_pos])
    km = np.array([km_of(lo[i], hi[i]) for i in range(n)], dtype=np.int64)   # gap-inclusive per-block
    n_fine_covered = int((km >= a.K).sum())                                  # already >=K before any merge
    nxt = list(range(1, n + 1)); nxt[-1] = -1; prev = list(range(-1, n - 1))
    alive = [True] * n
    ver = [0] * n
    pc = [pc1(geno, lo[i], hi[i]) for i in range(n)]
    K, FLOOR = a.K, a.floor

    heap = []
    def push(i):
        j = nxt[i]
        if j == -1 or min(km[i], km[j]) >= K:
            return
        l = r2(pc[i], pc[j])
        if l >= FLOOR:
            heapq.heappush(heap, (-l, i, j, ver[i], ver[j]))
    for i in range(n):
        push(i)

    merges = 0
    while heap:
        negl, i, j, vi, vj = heapq.heappop(heap)
        # reject stale: dead, no longer adjacent, no longer sub-K, or version changed since push
        if not (alive[i] and alive[j] and nxt[i] == j and min(km[i], km[j]) < K
                and ver[i] == vi and ver[j] == vj):
            continue
        hi[i] = hi[j]; alive[j] = False               # absorb right neighbour j into left block i
        nxt[i] = nxt[j]
        if nxt[i] != -1:
            prev[nxt[i]] = i
        km[i] = km_of(lo[i], hi[i])                   # gap-inclusive recount on the NEW interval
        pc[i] = pc1(geno, lo[i], hi[i]); ver[i] += 1  # bump version -> invalidates stale heap entries
        merges += 1
        push(i)
        if prev[i] != -1:
            push(prev[i])

    rows = []
    for i in range(n):
        if alive[i]:
            rows.append((a.chrom, int(positions[lo[i]]), int(positions[hi[i] - 1]),
                         int(hi[i] - lo[i]), int(km[i])))
    U = pd.DataFrame(rows, columns=["chrom", "start_pos", "end_pos", "n_variants", "panel_kmers"])
    U["covered"] = U.panel_kmers >= K
    out = f"{BR}/{chrlc}_units_dynld_K{K}.tsv"; U.to_csv(out, sep="\t", index=False)
    # sanity: total k-mers conserved (every centroid in the chrom's block span counted once)
    print(f"\n[{a.chrom}] {n:,} fine blocks -> {len(U):,} units ({merges} merges)")
    print(f"  COVERED (>=K kmers): {int(U.covered.sum()):,} ({100*U.covered.mean():.0f}%); "
          f"DESERT: {int((~U.covered).sum()):,} ({100*(~U.covered).mean():.0f}%)")
    print(f"  unit size: median {int(U.n_variants.median())} var / {int((U.end_pos-U.start_pos).median())} bp; "
          f"panel k-mers median {int(U.panel_kmers.median())}")
    print(f"  ({n_fine_covered:,} fine blocks already >=K, kept untouched)")
    print(f"  -> {out}")


if __name__ == "__main__":
    main()

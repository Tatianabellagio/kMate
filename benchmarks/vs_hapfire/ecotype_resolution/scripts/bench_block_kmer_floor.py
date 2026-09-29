#!/usr/bin/env python3
"""How many k-mers does a BLOCK need to resolve the ecotype mixture?

Block-mode kMate with the global fallback DISABLED (min_kmers_per_block=1) and
the global anchor OFF (anchor=0) -> each block estimates its founder mixture h
from ONLY its own k-mers. For a g0 (non-recombinant) pool the true per-block
mixture is the SAME for every block = the global pool_weights. So per-block
accuracy vs (#nonzero k-mers in the block) is a clean "k-mer floor" curve, and
comparing pools of different mixture size shows how the floor scales with the
number of ecotypes you are trying to differentiate.

Writes one row per local-fit block: (pool, true_n, cov, n_nz_kmers, n_panel_kmers,
cosine, pearson, jaccard_topN, mass_on_true, eff_n, h_rmse, would_fallback).
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/global/scratch/users/tbellg/kmate/src")
from kmate.per_sample_per_chrom import _count_and_load_kmer_pa_dense  # noqa: E402
from kmate.block_em import (define_windows, assign_kmers_to_blocks,  # noqa: E402
                            solve_em_per_block)

ROOT = Path("/global/scratch/users/tbellg/kmate")
PREFIX = str(ROOT / "benchmarks/p80/data/kmer_pa_p80_filt2/kmer_pa")


def load_truth_h(pool, panel, founders):
    """Map pool_weights.tsv -> truth h vector in kmer_pa founder order."""
    w = pd.read_csv(ROOT / f"benchmarks/{panel}/sims/{pool}/pool_weights.tsv", sep="\t")
    wmap = {str(f): float(x) for f, x in zip(w["founder"], w["weight"])}
    h = np.array([wmap.get(str(f), 0.0) for f in founders], dtype=float)
    s = h.sum()
    return (h / s) if s > 0 else h


def block_metrics(h_b, truth, true_set, true_n):
    hs = h_b.sum()
    if hs <= 0:
        return None
    h_b = h_b / hs
    # cosine + pearson to the truth mixture
    cos = float(h_b @ truth / (np.linalg.norm(h_b) * np.linalg.norm(truth) + 1e-12))
    pr = float(np.corrcoef(h_b, truth)[0, 1]) if h_b.std() > 0 else np.nan
    eff_n = float(1.0 / np.sum(h_b ** 2))
    top = set(np.argsort(h_b)[::-1][:true_n].tolist())
    jac = len(top & true_set) / len(top | true_set)
    mass = float(h_b[list(true_set)].sum())
    rmse = float(np.sqrt(np.mean((h_b - truth) ** 2)))
    return dict(cosine=cos, pearson=pr, jaccard_topN=jac,
                mass_on_true=mass, eff_n=eff_n, h_rmse=rmse)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--panel", default="p80")
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--window-bp", type=int, default=5000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--kmer-weight", default="inv_mb", choices=["uniform", "inv_mb"])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # true_n from pool name (..._nN_g0_...) or count from weights
    true_n = None
    for tok in a.pool.split("_"):
        if tok.startswith("n") and tok[1:].isdigit():
            true_n = int(tok[1:]); break
    cov = None
    for tok in a.pool.split("_"):
        if tok.startswith("cov") and tok[3:].isdigit():
            cov = int(tok[3:]); break

    reads = [str(ROOT / f"benchmarks/{a.panel}/sims/{a.pool}/reads/r1.fq"),
             str(ROOT / f"benchmarks/{a.panel}/sims/{a.pool}/reads/r2.fq")]

    t0 = time.time()
    print(f"[{a.pool}] counting k-mers ...", flush=True)
    kmer_pa_dense, counts, meta, coverage, F, K = _count_and_load_kmer_pa_dense(
        a.chrom, PREFIX, reads, a.threads)
    if kmer_pa_dense is None:
        sys.exit(f"chrom {a.chrom} missing for {a.pool}")
    founders = meta["founders"] if "founders" in meta else np.arange(F)
    truth = load_truth_h(a.pool, a.panel, founders)
    true_set = set(np.where(truth > 0)[0].tolist())
    if true_n is None:
        true_n = len(true_set)
    print(f"[{a.pool}] F={F} K={K:,} cov~{coverage:.1f} true_n={true_n} "
          f"(weights nonzero={len(true_set)}) count {time.time()-t0:.0f}s", flush=True)

    blocks = define_windows(meta["bubble_chrom"], meta["bubble_start"],
                            meta["bubble_end"], window_bp=a.window_bp)
    n_blocks = len(blocks)
    kmer_block = assign_kmers_to_blocks(meta["bubble_id"], meta["bubble_chrom"],
                                        meta["bubble_start"], meta["bubble_end"], blocks)
    omega = None
    if a.kmer_weight == "inv_mb":
        m_b = np.bincount(meta["bubble_id"])[meta["bubble_id"]].astype(np.float32)
        omega = (1.0 / m_b).astype(np.float32)

    nz = counts > 0
    n_nz_per_block = np.array([int(nz[kmer_block == b].sum()) for b in range(n_blocks)])
    n_panel_per_block = np.bincount(kmer_block[kmer_block >= 0], minlength=n_blocks)
    print(f"[{a.pool}] {n_blocks} windows of {a.window_bp}bp; "
          f"nonzero-kmer/block: median={np.median(n_nz_per_block):.0f} "
          f"max={n_nz_per_block.max()}; "
          f"blocks with >=1 nz kmer={int((n_nz_per_block>=1).sum())}", flush=True)

    # PURE LOCAL: no fallback (min=1), no anchor (0.0). status 0 = real local fit.
    t = time.time()
    h_blocks, status, global_h = solve_em_per_block(
        counts.astype(np.float32), kmer_pa_dense, kmer_block, n_blocks, coverage,
        min_kmers_per_block=1, global_anchor_weight=0.0, omega=omega,
        n_workers=max(1, a.threads // 2), verbose=True)
    print(f"[{a.pool}] pure-local block-EM {time.time()-t:.0f}s "
          f"({(status==0).sum()} local, {(status==2).sum()} empty)", flush=True)

    rows = []
    for b in range(n_blocks):
        if status[b] != 0:
            continue
        m = block_metrics(h_blocks[b], truth, true_set, true_n)
        if m is None:
            continue
        rows.append(dict(pool=a.pool, true_n=true_n, cov=cov,
                         block=b, n_nz_kmers=int(n_nz_per_block[b]),
                         n_panel_kmers=int(n_panel_per_block[b]),
                         would_fallback=int(n_nz_per_block[b] < 200), **m))
    df = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, sep="\t", index=False)
    print(f"[{a.pool}] wrote {len(df):,} local-fit blocks -> {a.out} "
          f"(total {time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

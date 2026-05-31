#!/usr/bin/env python3
"""Level-B index equivalence check: K_pa(ours index) vs K_pa(PanGenie index).

Both are built from the SAME arch3 merged_231 VCF + REF, filtered identically
(filt2inv); the only difference is the k-mer index (in-house `ours_Chr{N}` vs
PanGenie `pang_135_pangenie_index_Chr{N}`). If our indexer is downstream-equivalent
to PanGenie's, the two K_pa matrices should be near-identical: same founders, ~full
k-mer overlap, ~1.0 per-founder count correlation, ~100% carrier agreement on
shared k-mers. This is the downstream (Level-B) confirmation that Level-A parity
(index-level) carries through to the matrix the EM actually consumes.

Usage:
  python scripts/compare_index_pg_vs_ours.py            # all 5 chroms
  python scripts/compare_index_pg_vs_ours.py Chr1 Chr2  # subset
"""
import sys
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz

ROOT = Path(__file__).resolve().parents[1]
OURS = ROOT / "data/kmer_pa_231_arch3_filt2inv"
PG   = ROOT / "data/kmer_pa_231_arch3_pgidx_filt2inv"


def load(d, chrom):
    M = load_npz(d / f"kmer_pa_{chrom}.kmer_pa.npz").tocsc()
    m = np.load(d / f"kmer_pa_{chrom}.meta.npz", allow_pickle=True)
    return M, np.asarray(m["founders"]).astype(str), np.asarray(m["kmer_index"]).astype(str)


def compare_chrom(chrom):
    print(f"\n{'='*70}\n{chrom}\n{'='*70}")
    try:
        cnA, foA, kmA = load(OURS, chrom)   # ours
        cnB, foB, kmB = load(PG, chrom)     # pangenie
    except FileNotFoundError as e:
        print(f"  SKIP — missing matrix: {e}")
        return None

    print(f"  ours : F={len(foA)} K={cnA.shape[1]:,} nnz={cnA.nnz:,}")
    print(f"  pg   : F={len(foB)} K={cnB.shape[1]:,} nnz={cnB.nnz:,}")

    founders_ok = np.array_equal(foA, foB)
    print(f"  founders identical order: {founders_ok}")

    sA, sB = set(kmA), set(kmB)
    inter = sA & sB
    overlap_pct = 100 * len(inter) / max(len(sA), len(sB))
    print(f"  k-mer sets: ours={len(sA):,} pg={len(sB):,} shared={len(inter):,} "
          f"({overlap_pct:.2f}% of larger)  ours-only={len(sA-sB):,} pg-only={len(sB-sA):,}")

    count_corr = carrier_agree = identical_cols = None
    if founders_ok:
        kcA = np.asarray(cnA.sum(axis=1)).flatten()
        kcB = np.asarray(cnB.sum(axis=1)).flatten()
        count_corr = float(np.corrcoef(kcA, kcB)[0, 1])
        print(f"  per-founder k-mer count corr(ours,pg)={count_corr:.6f}  "
              f"median|Δ|={np.median(np.abs(kcA-kcB)):.0f}")

        if inter:
            iA = {k: i for i, k in enumerate(kmA)}
            iB = {k: i for i, k in enumerate(kmB)}
            shared = np.array(sorted(inter))
            samp = (shared if len(shared) <= 200_000
                    else np.random.default_rng(0).choice(shared, 200_000, replace=False))
            colA = cnA[:, [iA[k] for k in samp]].toarray().astype(bool)
            colB = cnB[:, [iB[k] for k in samp]].toarray().astype(bool)
            carrier_agree = float((colA == colB).mean())
            identical_cols = float((colA == colB).all(axis=0).mean())
            print(f"  carrier agreement on {len(samp):,} shared k-mers: "
                  f"per-cell {100*carrier_agree:.4f}%  fully-identical cols {100*identical_cols:.4f}%")

    # pass criterion: equivalent indexes -> high overlap + ~1.0 corr + ~100% carrier
    ok = (founders_ok and overlap_pct >= 98.0
          and (count_corr or 0) >= 0.999
          and (carrier_agree or 0) >= 0.999)
    print(f"  => {chrom}: {'PASS (equivalent)' if ok else 'REVIEW (see numbers above)'}")
    return ok


if __name__ == "__main__":
    chroms = sys.argv[1:] or ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
    results = {c: compare_chrom(c) for c in chroms}
    print(f"\n{'='*70}\nSUMMARY\n{'='*70}")
    for c, r in results.items():
        print(f"  {c}: {'PASS' if r else ('SKIP' if r is None else 'REVIEW')}")
    print("\nPASS on all chroms => PanGenie and our index produce the same K_pa "
          "(downstream-equivalent); arch3 K_pa can use the in-house index.")

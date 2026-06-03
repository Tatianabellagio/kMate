#!/usr/bin/env python3
"""Generic K_pa-vs-K_pa comparator for two kmer_pa dirs.

Same metrics as compare_index_pg_vs_ours.py (k-mer set overlap, per-founder
count correlation, carrier agreement on shared k-mers) but with arbitrary dirs
and labels — so it works for nocap-vs-capped (this session's no-caps check),
not just ours-vs-PG.

For the no-caps check specifically: the capped selection is a sorted-order
truncation of the no-cap selection, so the capped k-mer set should be a near
subset of the no-cap set. The headline number is therefore `b-only` (k-mers the
caps DROPPED) and whether carriers on shared k-mers stay identical.

Usage:
  python scripts/compare_kmer_pa_dirs.py \
      --a data/kmer_pa_231_arch3_nocap_filt2inv --a-label nocap \
      --b data/kmer_pa_231_arch3_filt2inv       --b-label capped \
      --chroms Chr1
"""
import argparse
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz


def load(d: Path, chrom: str):
    M = load_npz(d / f"kmer_pa_{chrom}.kmer_pa.npz").tocsc()
    m = np.load(d / f"kmer_pa_{chrom}.meta.npz", allow_pickle=True)
    return M, np.asarray(m["founders"]).astype(str), np.asarray(m["kmer_index"]).astype(str)


def compare_chrom(chrom, dirA, dirB, labA, labB):
    print(f"\n{'='*70}\n{chrom}\n{'='*70}")
    try:
        cnA, foA, kmA = load(dirA, chrom)
        cnB, foB, kmB = load(dirB, chrom)
    except FileNotFoundError as e:
        print(f"  SKIP — missing matrix: {e}")
        return None

    print(f"  {labA:8s}: F={len(foA)} K={cnA.shape[1]:,} nnz={cnA.nnz:,}")
    print(f"  {labB:8s}: F={len(foB)} K={cnB.shape[1]:,} nnz={cnB.nnz:,}")
    dK = cnA.shape[1] - cnB.shape[1]
    print(f"  ΔK ({labA}-{labB}) = {dK:+,}  "
          f"({100*dK/cnB.shape[1]:+.2f}% vs {labB})   "
          f"Δnnz = {cnA.nnz - cnB.nnz:+,}")

    founders_ok = np.array_equal(foA, foB)
    print(f"  founders identical order: {founders_ok}")

    sA, sB = set(kmA), set(kmB)
    inter = sA & sB
    aonly, bonly = sA - sB, sB - sA
    print(f"  k-mer sets: {labA}={len(sA):,} {labB}={len(sB):,} shared={len(inter):,}")
    print(f"    {labA}-only={len(aonly):,}   {labB}-only={len(bonly):,}")
    if sB:
        print(f"    fraction of {labB} retained in {labA} (subset check): "
              f"{100*len(inter)/len(sB):.4f}%")

    if founders_ok:
        kcA = np.asarray(cnA.sum(axis=1)).flatten()
        kcB = np.asarray(cnB.sum(axis=1)).flatten()
        corr = float(np.corrcoef(kcA, kcB)[0, 1])
        print(f"  per-founder k-mer count corr({labA},{labB})={corr:.6f}  "
              f"median|Δ|={np.median(np.abs(kcA-kcB)):.0f}  "
              f"mean {labA}/{labB} ratio={float(kcA.sum())/max(1,kcB.sum()):.4f}")

        if inter:
            iA = {k: i for i, k in enumerate(kmA)}
            iB = {k: i for i, k in enumerate(kmB)}
            shared = np.array(sorted(inter))
            samp = (shared if len(shared) <= 200_000
                    else np.random.default_rng(0).choice(shared, 200_000, replace=False))
            colA = cnA[:, [iA[k] for k in samp]].toarray().astype(bool)
            colB = cnB[:, [iB[k] for k in samp]].toarray().astype(bool)
            print(f"  carrier agreement on {len(samp):,} shared k-mers: "
                  f"per-cell {100*float((colA==colB).mean()):.4f}%  "
                  f"fully-identical cols {100*float((colA==colB).all(axis=0).mean()):.4f}%")
    return True


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True, help="dir A (e.g. nocap)")
    ap.add_argument("--b", required=True, help="dir B (e.g. capped baseline)")
    ap.add_argument("--a-label", default="A")
    ap.add_argument("--b-label", default="B")
    ap.add_argument("--chroms", default="Chr1,Chr2,Chr3,Chr4,Chr5")
    args = ap.parse_args()
    dirA, dirB = Path(args.a), Path(args.b)
    chroms = [c for c in args.chroms.split(",") if c]
    for c in chroms:
        compare_chrom(c, dirA, dirB, args.a_label, args.b_label)
    print(f"\n{'='*70}\nNOTE\n{'='*70}")
    print(f"  '{args.b_label}-only' = k-mers the caps DROPPED that no-caps keeps "
          f"(expected ~0 if {args.b_label} ⊆ {args.a_label}).")
    print(f"  '{args.a_label}-only' = extra k-mers no-caps adds over {args.b_label}.")
    print("  100% per-cell carrier agreement on shared k-mers => caps only change "
          "k-mer COUNT, not which founder carries which shared k-mer.")

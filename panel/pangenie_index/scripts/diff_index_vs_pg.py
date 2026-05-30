#!/usr/bin/env python -u
"""
Level A index-parity diff: our standalone k-mer index vs PanGenie's production
index, bubble by bubble.

Both indexes are built from the SAME inputs (pang_135 `pang_1001gplus_all.vcf.gz`
+ `TAIR10.chr.fa`, k=31), so this is the true apples-to-apples indexer test:
any difference is attributable to the two indexers, not to differing inputs.

For each chromosome we stream both `kmers.tsv.gz` (position-sorted) in a
two-pointer merge keyed on (start, end). For each matched bubble we compare the
`unique_kmers` field and the `unique_kmers_overhang` field as ORDER-INSENSITIVE
sets (PanGenie's emission order is a std::map artefact; downstream
`build_kmer_pa.py` consumes them as an unordered set anyway).

Per field we bucket each bubble into:
  identical        : same set (incl. both empty/`nan`)
  ours_subset      : ours ⊂ PG     (ours kept fewer; never a wrong k-mer)
  pg_subset        : PG ⊂ ours      (ours kept extra; PG never emitted them)
  overlap          : non-empty intersection, neither a subset
  disjoint         : non-empty sets, empty intersection
  ours_empty_only  : ours nan/empty, PG non-empty
  pg_empty_only    : PG nan/empty, ours non-empty

Usage:
  diff_index_vs_pg.py --pg <pg_prefix> --ours <our_prefix> --chroms Chr1[,Chr2,...]
                      [--examples N] [--out report.tsv]
where a prefix expands to `<prefix>_<chrom>_kmers.tsv.gz`.
"""
from __future__ import annotations
import argparse
import gzip
import sys
from collections import Counter


def parse_kmer_field(field: str) -> frozenset:
    """A kmers / overhang TSV field -> frozenset of k-mers.
    Empty field, 'nan', and '' all map to the empty set (PanGenie emits the
    literal 'nan' for bubbles with no surviving k-mers; so do we)."""
    if not field or field == "nan":
        return frozenset()
    return frozenset(field.split(","))


def iter_bubbles(path: str):
    """Yield ((start, end), kmers_set, overhang_set) per data row, in file order."""
    with gzip.open(path, "rt") as f:
        next(f)  # header
        for line in f:
            p = line.rstrip("\n").split("\t")
            yield (int(p[1]), int(p[2])), parse_kmer_field(p[3]), \
                  parse_kmer_field(p[4]) if len(p) > 4 else frozenset()


def classify(ours: frozenset, pg: frozenset) -> str:
    if ours == pg:
        return "identical"
    if not ours and pg:
        return "ours_empty_only"
    if ours and not pg:
        return "pg_empty_only"
    if ours < pg:
        return "ours_subset"
    if pg < ours:
        return "pg_subset"
    if ours & pg:
        return "overlap"
    return "disjoint"


def diff_chrom(pg_path: str, our_path: str, chrom: str, n_examples: int):
    """Two-pointer merge on (start,end). Returns (stats dict, examples list)."""
    kmer_buckets = Counter()
    oh_buckets = Counter()
    coord_only_pg = 0
    coord_only_ours = 0
    matched = 0
    # k-mer-level tallies on matched bubbles
    n_kmers_pg = n_kmers_ours = inter_kmers = 0
    examples = []

    gp = iter_bubbles(pg_path)
    go = iter_bubbles(our_path)
    pg = next(gp, None)
    ou = next(go, None)
    while pg is not None or ou is not None:
        if ou is None or (pg is not None and pg[0] < ou[0]):
            coord_only_pg += 1
            pg = next(gp, None)
            continue
        if pg is None or ou[0] < pg[0]:
            coord_only_ours += 1
            ou = next(go, None)
            continue
        # coordinates match
        matched += 1
        (_, pg_k, pg_oh) = pg
        (_, ou_k, ou_oh) = ou
        kc = classify(ou_k, pg_k)
        oc = classify(ou_oh, pg_oh)
        kmer_buckets[kc] += 1
        oh_buckets[oc] += 1
        n_kmers_pg += len(pg_k)
        n_kmers_ours += len(ou_k)
        inter_kmers += len(ou_k & pg_k)
        if kc != "identical" and len(examples) < n_examples:
            examples.append((pg[0], kc,
                             sorted(ou_k - pg_k)[:3], sorted(pg_k - ou_k)[:3]))
        pg = next(gp, None)
        ou = next(go, None)

    stats = {
        "chrom": chrom, "matched": matched,
        "coord_only_pg": coord_only_pg, "coord_only_ours": coord_only_ours,
        "kmer_buckets": dict(kmer_buckets), "oh_buckets": dict(oh_buckets),
        "n_kmers_pg": n_kmers_pg, "n_kmers_ours": n_kmers_ours,
        "inter_kmers": inter_kmers,
    }
    return stats, examples


def pct(n, d):
    return f"{100.0 * n / d:.3f}%" if d else "n/a"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--pg", required=True, help="PG index prefix (-> <prefix>_<chrom>_kmers.tsv.gz)")
    ap.add_argument("--ours", required=True, help="our index prefix")
    ap.add_argument("--chroms", required=True, help="comma-separated, e.g. Chr1,Chr2")
    ap.add_argument("--examples", type=int, default=8, help="example differing bubbles to print per chrom")
    args = ap.parse_args()

    for chrom in args.chroms.split(","):
        pg_path = f"{args.pg}_{chrom}_kmers.tsv.gz"
        our_path = f"{args.ours}_{chrom}_kmers.tsv.gz"
        print(f"\n{'='*72}\n{chrom}: PG={pg_path}\n      OURS={our_path}\n{'='*72}", file=sys.stderr)
        stats, examples = diff_chrom(pg_path, our_path, chrom, args.examples)

        m = stats["matched"]
        print(f"\n## {chrom}")
        print(f"matched bubbles (same coords): {m:,}")
        print(f"coord-only-in-PG:   {stats['coord_only_pg']:,}")
        print(f"coord-only-in-ours: {stats['coord_only_ours']:,}")
        print(f"\nunique_kmers field (per matched bubble):")
        for k in ("identical", "ours_subset", "pg_subset", "overlap", "disjoint",
                  "ours_empty_only", "pg_empty_only"):
            v = stats["kmer_buckets"].get(k, 0)
            print(f"  {k:18s} {v:>10,}  ({pct(v, m)})")
        print(f"\noverhang field (per matched bubble):")
        for k in ("identical", "ours_subset", "pg_subset", "overlap", "disjoint",
                  "ours_empty_only", "pg_empty_only"):
            v = stats["oh_buckets"].get(k, 0)
            print(f"  {k:18s} {v:>10,}  ({pct(v, m)})")
        print(f"\nk-mer totals on matched bubbles:")
        print(f"  PG   unique_kmers: {stats['n_kmers_pg']:,}")
        print(f"  OURS unique_kmers: {stats['n_kmers_ours']:,}  "
              f"(ours/PG = {stats['n_kmers_ours']/stats['n_kmers_pg']:.4f})"
              if stats['n_kmers_pg'] else "")
        print(f"  intersection:      {stats['inter_kmers']:,}  "
              f"(recall vs PG = {pct(stats['inter_kmers'], stats['n_kmers_pg'])}, "
              f"precision = {pct(stats['inter_kmers'], stats['n_kmers_ours'])})")
        if examples:
            print(f"\nexample differing bubbles (coord, bucket, ours-only[:3], pg-only[:3]):")
            for coord, bucket, oo, po in examples:
                print(f"  {coord} {bucket}\n     ours-only: {oo}\n     pg-only:   {po}")


if __name__ == "__main__":
    main()

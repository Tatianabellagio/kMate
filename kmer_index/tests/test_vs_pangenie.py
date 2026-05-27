#!/usr/bin/env python3
"""
Regression test: compare build_kmers_tsv.py output against PanGenie-index output.

Runs against the pre-built test fixtures in:
  timing_pg2/pg_Chr1_kmers.tsv.gz          (PanGenie-index v4.2.1)
  timing_ours_v2/ours_prof_Chr1_kmers.tsv.gz  (our implementation)

Expected result (2 Mb Chr1 diploid test panel, 1985 bubbles):
  - All 1985 bubble coordinates match exactly
  - ≥99% of bubbles have identical unique_kmers AND unique_kmers_overhang sets
  - Remaining diffs: PG emits boundary-spanning overhang k-mers that straddle
    bubble_start (off-by-one in PG's left-flank window); our output excludes them,
    which is semantically correct.
  - No cases where our output has a k-mer PG does not (we are a strict subset of PG
    on the main kmer column, never a superset).

Usage:
  python tests/test_vs_pangenie.py
  # or from repo root:
  python -m pytest kmer_index/tests/test_vs_pangenie.py -v
"""
import gzip
import sys
from pathlib import Path

WORK = Path(__file__).resolve().parent.parent

PG_TSV  = WORK / "timing_pg2"   / "pg_Chr1_kmers.tsv.gz"
OUR_TSV = WORK / "timing_ours_v2" / "ours_prof_Chr1_kmers.tsv.gz"


def parse_tsv(path: Path) -> dict:
    rows = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = (parts[0], int(parts[1]), int(parts[2]))
            kmers = frozenset(parts[3].split(",")) if parts[3] != "nan" else frozenset()
            oh    = frozenset(parts[4].split(",")) if parts[4] != "nan" else frozenset()
            rows[key] = (kmers, oh)
    return rows


def run():
    if not PG_TSV.exists():
        print(f"SKIP: PG fixture not found at {PG_TSV}")
        sys.exit(0)
    if not OUR_TSV.exists():
        print(f"SKIP: our fixture not found at {OUR_TSV}")
        sys.exit(0)

    pg  = parse_tsv(PG_TSV)
    our = parse_tsv(OUR_TSV)

    # --- coordinate completeness ---
    only_pg  = set(pg) - set(our)
    only_our = set(our) - set(pg)
    assert len(only_pg)  == 0, f"Bubbles in PG not in ours: {list(only_pg)[:5]}"
    assert len(only_our) == 0, f"Bubbles in ours not in PG: {list(only_our)[:5]}"

    shared = set(pg) & set(our)
    n = len(shared)
    assert n == 1985, f"Expected 1985 bubbles, got {n}"

    exact          = 0
    oh_pg_extra    = 0   # PG has an overhang kmer we don't — boundary-spanning quirk
    kmer_our_extra = []  # we have a main kmer PG doesn't — should be 0

    for key in shared:
        pk, po = pg[key]
        ok, oo = our[key]
        if pk == ok and po == oo:
            exact += 1
            continue
        # Main kmer column: we must never be a superset of PG
        extra_in_ours = ok - pk
        if extra_in_ours:
            kmer_our_extra.append((key, extra_in_ours))
        # Overhang column: PG-only extras are the boundary-spanning quirk
        oh_only_pg  = po - oo
        oh_only_our = oo - po
        if oh_only_pg and not oh_only_our and not (ok - pk) and not (pk - ok):
            oh_pg_extra += 1

    pct_exact = exact / n * 100

    print(f"Bubbles total:          {n}")
    print(f"Exact match:            {exact} ({pct_exact:.1f}%)")
    print(f"PG boundary-overhang:   {oh_pg_extra}  (expected; PG off-by-one at bubble_start)")
    print(f"Our-only main kmers:    {len(kmer_our_extra)}  (must be 0)")
    print()

    assert pct_exact >= 99.0, f"Exact match {pct_exact:.1f}% < 99% threshold"
    assert len(kmer_our_extra) == 0, (
        f"Our output has {len(kmer_our_extra)} bubble(s) with k-mers PG doesn't have: "
        f"{kmer_our_extra[:3]}"
    )
    assert exact + oh_pg_extra == n, (
        f"Unexpected diff type: {n - exact - oh_pg_extra} bubble(s) have unexplained differences"
    )

    print("PASS")


if __name__ == "__main__":
    run()

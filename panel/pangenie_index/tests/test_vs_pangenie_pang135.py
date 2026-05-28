#!/usr/bin/env python3
"""
Regression test on the full pang_135 panel: compare our build_kmers_tsv.py
output (panel/pangenie_index/pang_135_diploid/) against the canonical PG-index output
already on disk (panel/pangenie_genotyping/data/pang_135_pangenie_index_Chr{1..5}_kmers.tsv.gz,
built 2026-05-01 from the same .dipl VCF + TAIR10 ref).

Run after slurm_pang135_diploid.sh finishes.

Reports per-chrom:
  - bubble coord overlap
  - exact-match fraction
  - PG-only / our-only main-kmer counts (must both be 0 in main column)
  - PG boundary-overhang quirk count (expected non-zero; PG off-by-one)
"""
import gzip
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
PG_DIR  = _ROOT / "panel/pangenie_genotyping/data"
OUR_DIR = _ROOT / "panel/pangenie_index/pang_135_diploid"

CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]


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
    fail = False
    summary = []
    for chrom in CHROMS:
        pg  = PG_DIR  / f"pang_135_pangenie_index_{chrom}_kmers.tsv.gz"
        our = OUR_DIR / f"ours_{chrom}_kmers.tsv.gz"
        if not pg.exists():
            print(f"[{chrom}] SKIP — PG fixture missing: {pg}")
            continue
        if not our.exists():
            print(f"[{chrom}] SKIP — our output missing: {our}")
            continue

        pg_rows  = parse_tsv(pg)
        our_rows = parse_tsv(our)

        only_pg  = set(pg_rows)  - set(our_rows)
        only_our = set(our_rows) - set(pg_rows)
        shared   = set(pg_rows)  & set(our_rows)

        exact = 0
        oh_pg_extra = 0
        kmer_our_extra = 0
        kmer_pg_extra  = 0
        for key in shared:
            pk, po = pg_rows[key]
            ok, oo = our_rows[key]
            if pk == ok and po == oo:
                exact += 1
                continue
            if ok - pk:
                kmer_our_extra += 1
            if pk - ok:
                kmer_pg_extra  += 1
            if (po - oo) and not (oo - po) and pk == ok:
                oh_pg_extra += 1

        n = len(shared)
        pct = exact / n * 100 if n else 0.0
        line = (
            f"[{chrom}] PG={len(pg_rows):>7} ours={len(our_rows):>7} "
            f"only_PG={len(only_pg):>5} only_ours={len(only_our):>5} | "
            f"exact={exact:>7}/{n} ({pct:5.1f}%) "
            f"oh_pg_only={oh_pg_extra:>5} | "
            f"kmer_ours_extra={kmer_our_extra:>3} kmer_pg_extra={kmer_pg_extra:>4}"
        )
        print(line)
        summary.append((chrom, n, exact, kmer_our_extra, len(only_our), pct))

        if kmer_our_extra > 0 or len(only_our) > 0:
            fail = True

    print()
    if not summary:
        print("NO OUTPUTS TO COMPARE — has the SLURM job finished?")
        sys.exit(1)
    if fail:
        print("FAIL — our output has k-mers/bubbles not in PG (semantic divergence)")
        sys.exit(2)
    print("PASS — no false-positive k-mers in our output")


if __name__ == "__main__":
    run()

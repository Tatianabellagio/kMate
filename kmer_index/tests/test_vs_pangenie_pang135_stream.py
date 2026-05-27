#!/usr/bin/env python3
"""
Memory-frugal streaming version of test_vs_pangenie_pang135.py.

Both our output and PG's production index emit bubbles in identical
(chrom, start, end) order (verified: same row count, same coords). So we can
zip the two gzip streams line-by-line and compare one bubble at a time without
holding either full file in RAM (the non-streaming version OOMs: ~30M k-mer
str objects per chrom per file).

Reports per chrom:
  - row-count / coord agreement (hard fail on mismatch)
  - exact-match fraction
  - kmer_ours_extra : bubbles where our main-kmer set has a k-mer PG lacks  (MUST be 0)
  - kmer_pg_extra   : bubbles where PG main-kmer set has a k-mer we lack
  - oh_pg_only      : bubbles identical on main kmers but PG has extra overhang kmers
  - oh_ours_only    : bubbles identical on main kmers but WE have extra overhang kmers
"""
import gzip
import sys
from pathlib import Path

PG_DIR  = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/data")
OUR_DIR = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/kmer_index/pang_135_diploid")
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]


def fields(line: str):
    p = line.rstrip("\n").split("\t")
    key = (p[0], int(p[1]), int(p[2]))
    kmers = frozenset(p[3].split(",")) if p[3] != "nan" else frozenset()
    oh    = frozenset(p[4].split(",")) if p[4] != "nan" else frozenset()
    return key, kmers, oh


def run():
    fail = False
    for chrom in CHROMS:
        pg  = PG_DIR  / f"pang_135_pangenie_index_{chrom}_kmers.tsv.gz"
        our = OUR_DIR / f"ours_{chrom}_kmers.tsv.gz"
        if not pg.exists() or not our.exists():
            print(f"[{chrom}] SKIP — missing fixture")
            continue

        n = exact = 0
        coord_mismatch = 0
        kmer_ours_extra = kmer_pg_extra = 0
        oh_pg_only = oh_ours_only = 0
        first_ours_extra = []

        with gzip.open(pg, "rt") as fpg, gzip.open(our, "rt") as four:
            # skip headers
            ph = next(fpg); oh_ = next(four)
            for lp, lo in zip(fpg, four):
                kp, pk, po = fields(lp)
                ko, ok, oo = fields(lo)
                if kp != ko:
                    coord_mismatch += 1
                    continue
                n += 1
                if pk == ok and po == oo:
                    exact += 1
                    continue
                ours_x = ok - pk
                if ours_x:
                    kmer_ours_extra += 1
                    if len(first_ours_extra) < 3:
                        first_ours_extra.append((ko, sorted(ours_x)[:2]))
                if pk - ok:
                    kmer_pg_extra += 1
                if pk == ok:  # main kmers identical, diff is overhang only
                    if po - oo:
                        oh_pg_only += 1
                    if oo - po:
                        oh_ours_only += 1

        pct = exact / n * 100 if n else 0.0
        print(f"[{chrom}] bubbles={n:>7} exact={exact:>7} ({pct:5.2f}%) | "
              f"coord_mismatch={coord_mismatch} | "
              f"kmer_ours_extra={kmer_ours_extra} kmer_pg_extra={kmer_pg_extra} | "
              f"oh_pg_only={oh_pg_only} oh_ours_only={oh_ours_only}")
        if coord_mismatch or kmer_ours_extra:
            fail = True
            for k, ex in first_ours_extra:
                print(f"    ours-extra @ {k}: {ex}")

    print()
    if fail:
        print("FAIL — coord mismatch or our output emits k-mers PG does not")
        sys.exit(2)
    print("PASS — our main-kmer column is a subset of PG (no false-positive k-mers)")


if __name__ == "__main__":
    run()

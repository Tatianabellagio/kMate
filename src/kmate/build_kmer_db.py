"""`kmate build-kmer-db` — build a Jellyfish k-mer DB from a read pool, once.

Scanning the reads is the expensive step. When the same pool is queried against
several disjoint k-mer sets (one per chromosome), build the DB once here and pass
it to `kmate run --kmer-db <db>`, instead of re-scanning the reads for every
chromosome. Counts are byte-identical to the per-chrom count path (same canonical
hash). See the count-once / query-per-chrom note in kmer_count.build_kmer_db.
"""
from __future__ import annotations
import argparse
from .kmer_count import build_kmer_db


def main():
    ap = argparse.ArgumentParser(
        prog="kmate build-kmer-db",
        description="Build a canonical Jellyfish k-mer DB from a read pool (count once, "
                    "then query per chromosome with `kmate run --kmer-db`).")
    ap.add_argument("--reads", required=True, nargs="+",
                    help="FASTQ/FASTA file(s) (plain or .gz), or a single .bam.")
    ap.add_argument("--out", required=True, help="output Jellyfish DB path (e.g. pool.jf)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--hash-size", default="3G",
                    help="initial Jellyfish hash size (auto-grows); lower for small "
                         "pools / memory-capped jobs. Production: 3G.")
    args = ap.parse_args()
    reads = args.reads if len(args.reads) > 1 else args.reads[0]
    db = build_kmer_db(reads, args.out, k=31, threads=args.threads,
                       hash_size=args.hash_size)
    print(f"Wrote {db}")


if __name__ == "__main__":
    main()

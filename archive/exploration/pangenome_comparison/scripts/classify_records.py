"""
classify_records.py — read a cactus pangenome VCF, classify every record × ALT
into size classes and emit per-(chrom, size_class) counts to TSV.

Streams via `bcftools query` (stdin); avoids pysam header-parse cost on 5M+
record files. Writes:
    <out>.records.tsv     per-record summary (chrom, pos, n_alts, max_size, max_class)
    <out>.size_class.tsv  aggregated counts per (chrom, size_class), counted at the ALT level
                          (so a multi-allelic bubble contributes once per ALT)

Size classes (max(ref_len, alt_len)):
    SNP          : ref_len == 1 and alt_len == 1
    small_indel  : max < 50
    small_sv     : 50 <= max < 500
    medium_sv    : 500 <= max < 5000
    large_sv     : max >= 5000
"""
import sys, argparse, gzip
from collections import Counter

CLASSES = ["SNP", "small_indel", "small_sv", "medium_sv", "large_sv"]

def classify(ref_len: int, alt_len: int) -> str:
    if ref_len == 1 and alt_len == 1:
        return "SNP"
    m = max(ref_len, alt_len)
    if m < 50:    return "small_indel"
    if m < 500:   return "small_sv"
    if m < 5000:  return "medium_sv"
    return "large_sv"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in",    dest="inp",  required=True, help="bcftools query stream: chrom\\tpos\\tref\\talt")
    ap.add_argument("--out",   required=True, help="output prefix")
    ap.add_argument("--label", required=True, help="panel label (e.g. pang_82, pang_135)")
    args = ap.parse_args()

    # (chrom, size_class) -> count of ALTs
    alt_counts = Counter()
    # chrom -> [n_records, n_alts]
    chrom_counts = Counter()
    chrom_alts = Counter()
    n_records = 0
    n_alts_total = 0
    rec_path = args.out + ".records.tsv"
    opener = open
    with opener(args.inp) as inp, open(rec_path, "w") as recf:
        recf.write("chrom\tpos\tn_alts\tmax_size\tmax_class\n")
        for line in inp:
            line = line.rstrip("\n")
            if not line or line.startswith("#"): continue
            fields = line.split("\t")
            if len(fields) < 4: continue
            chrom, pos, ref, alts = fields[0], fields[1], fields[2], fields[3]
            n_records += 1
            chrom_counts[chrom] += 1
            ref_len = len(ref)
            max_size = 0
            max_class = "SNP"
            n_alts = 0
            for alt in alts.split(","):
                if alt == "" or alt == ".": continue
                n_alts += 1
                alt_len = len(alt)
                cls = classify(ref_len, alt_len)
                alt_counts[(chrom, cls)] += 1
                m = max(ref_len, alt_len)
                if m > max_size:
                    max_size = m
                    max_class = cls
            chrom_alts[chrom] += n_alts
            n_alts_total += n_alts
            recf.write(f"{chrom}\t{pos}\t{n_alts}\t{max_size}\t{max_class}\n")

    sc_path = args.out + ".size_class.tsv"
    with open(sc_path, "w") as f:
        f.write("panel\tchrom\tsize_class\tn_alts\n")
        chroms = sorted(chrom_counts.keys())
        for chrom in chroms:
            for cls in CLASSES:
                f.write(f"{args.label}\t{chrom}\t{cls}\t{alt_counts.get((chrom, cls), 0)}\n")

    chrom_path = args.out + ".chrom.tsv"
    with open(chrom_path, "w") as f:
        f.write("panel\tchrom\tn_records\tn_alts\n")
        for chrom in sorted(chrom_counts.keys()):
            f.write(f"{args.label}\t{chrom}\t{chrom_counts[chrom]}\t{chrom_alts[chrom]}\n")

    print(f"[done] {args.label}: {n_records:,} records, {n_alts_total:,} ALTs")
    print(f"  wrote {rec_path}")
    print(f"  wrote {sc_path}")
    print(f"  wrote {chrom_path}")

if __name__ == "__main__":
    main()

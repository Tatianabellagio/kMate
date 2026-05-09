"""
loo_concordance.py — per-record concordance between a PanGenie genotype VCF
(short reads → catalog) and the cactus truth VCF (long-read assembly) for one
LOO sample.

Reports genotype concordance (GC) and non-Reference Discordance (nRD), broken
out by variant size class:
    SNP          : ref_len == 1 and alt_len == 1
    small_indel  : max(ref_len, alt_len) < 50
    small_sv     : 50 <= max < 500
    medium_sv    : 500 <= max < 5000
    large_sv     : max >= 5000

Usage:
  python loo_concordance.py \\
    --pangenie-vcf data/loo_genotyped/9947.vcf.gz \\
    --truth-vcf    /home/tbellagio/scratch/pang/.../pang_1001gplus_all.vcf.gz \\
    --sample       9947 \\
    --truth-sample 100612 \\
    --out          data/loo_concordance/9947

Outputs:
  <out>_summary.tsv    one row per size class with GC, nRD, counts
  <out>_records.tsv.gz per-record (chrom, pos, size_class, truth_gt, pg_gt, match)
"""
from __future__ import annotations
import argparse, gzip, sys
from pathlib import Path
import pysam


SIZE_CLASSES = ["SNP", "small_indel", "small_sv", "medium_sv", "large_sv"]


def classify(ref_len: int, alt_len: int) -> str:
    if ref_len == 1 and alt_len == 1:
        return "SNP"
    m = max(ref_len, alt_len)
    if m < 50:    return "small_indel"
    if m < 500:   return "small_sv"
    if m < 5000:  return "medium_sv"
    return "large_sv"


def gt_to_dose(gt) -> int | None:
    """0|0 → 0, 0|1 / 1|0 → 1, 1|1 → 2, ./. → None."""
    if gt is None or len(gt) == 0:
        return None
    vals = [g for g in gt if g is not None]
    if not vals:
        return None
    return sum(1 for g in vals if g > 0)


def index_truth(truth_vcf: str, sample: str) -> dict:
    """Return dict (chrom, pos, ref, alt) → truth dose for one sample."""
    idx = {}
    v = pysam.VariantFile(truth_vcf)
    if sample not in v.header.samples:
        sys.exit(f"ERROR: sample '{sample}' not in truth VCF samples")
    skipped_no_alt = 0
    for rec in v.fetch():
        if not rec.alts:
            skipped_no_alt += 1
            continue
        gt = rec.samples[sample]["GT"]
        dose = gt_to_dose(gt)
        # Index every (alt) — for multi-allelic, dose mapping per alt is fragile;
        # treat any non-ref allele as "carrier" (dose >= 1)
        for alt in rec.alts:
            idx[(rec.chrom, rec.pos, rec.ref, alt)] = dose
    print(f"[truth] {len(idx):,} (chrom,pos,ref,alt) keys, {skipped_no_alt} skipped (no alt)")
    return idx


def compare(pangenie_vcf: str, sample: str, truth_idx: dict, out_prefix: str):
    v = pysam.VariantFile(pangenie_vcf)
    if sample not in v.header.samples:
        sys.exit(f"ERROR: sample '{sample}' not in pangenie VCF samples")

    # tallies per size_class
    n_total       = {c: 0 for c in SIZE_CLASSES}
    n_matched     = {c: 0 for c in SIZE_CLASSES}   # variants where both truth+pg have a call
    n_concordant  = {c: 0 for c in SIZE_CLASSES}   # truth_dose == pg_dose
    n_nonref_either = {c: 0 for c in SIZE_CLASSES} # truth_dose>0 OR pg_dose>0
    n_nonref_disc   = {c: 0 for c in SIZE_CLASSES} # nonref_either AND truth_dose != pg_dose
    n_truth_missing = {c: 0 for c in SIZE_CLASSES}
    n_pg_missing    = {c: 0 for c in SIZE_CLASSES}

    rec_path = out_prefix + "_records.tsv.gz"
    Path(rec_path).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(rec_path, "wt") as recf:
        recf.write("chrom\tpos\tref_len\talt_len\tsize_class\ttruth_dose\tpg_dose\tmatch\n")
        for rec in v.fetch():
            if not rec.alts:
                continue
            for alt in rec.alts:
                cls = classify(len(rec.ref), len(alt))
                n_total[cls] += 1
                truth_dose = truth_idx.get((rec.chrom, rec.pos, rec.ref, alt))
                pg_gt = rec.samples[sample]["GT"]
                pg_dose = gt_to_dose(pg_gt)
                if truth_dose is None:
                    n_truth_missing[cls] += 1
                if pg_dose is None:
                    n_pg_missing[cls] += 1
                if truth_dose is None or pg_dose is None:
                    continue
                n_matched[cls] += 1
                concord = (truth_dose == pg_dose)
                if concord:
                    n_concordant[cls] += 1
                if truth_dose > 0 or pg_dose > 0:
                    n_nonref_either[cls] += 1
                    if not concord:
                        n_nonref_disc[cls] += 1
                recf.write(f"{rec.chrom}\t{rec.pos}\t{len(rec.ref)}\t{len(alt)}\t{cls}\t"
                           f"{truth_dose}\t{pg_dose}\t{int(concord)}\n")

    # Summary
    sum_path = out_prefix + "_summary.tsv"
    with open(sum_path, "w") as f:
        f.write("size_class\tn_total\tn_matched\tn_concordant\tGC\t"
                "n_nonref_either\tn_nonref_disc\tnRD\tn_truth_missing\tn_pg_missing\n")
        for cls in SIZE_CLASSES:
            t = n_total[cls]; m = n_matched[cls]; c = n_concordant[cls]
            ne = n_nonref_either[cls]; nd = n_nonref_disc[cls]
            tm = n_truth_missing[cls]; pm = n_pg_missing[cls]
            gc = c / m if m > 0 else float("nan")
            nrd = nd / ne if ne > 0 else float("nan")
            f.write(f"{cls}\t{t}\t{m}\t{c}\t{gc:.4f}\t{ne}\t{nd}\t{nrd:.4f}\t{tm}\t{pm}\n")

    print(f"\nWrote {sum_path}")
    print(f"Wrote {rec_path}")
    print("\n=== Headline ===")
    print(f"{'class':<14}{'n_total':>10}{'n_matched':>12}{'GC':>9}{'nRD':>9}")
    for cls in SIZE_CLASSES:
        t = n_total[cls]; m = n_matched[cls]; c = n_concordant[cls]
        ne = n_nonref_either[cls]; nd = n_nonref_disc[cls]
        gc = c / m if m > 0 else float("nan")
        nrd = nd / ne if ne > 0 else float("nan")
        print(f"{cls:<14}{t:>10}{m:>12}{gc:>9.4f}{nrd:>9.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pangenie-vcf", required=True)
    ap.add_argument("--truth-vcf",    required=True)
    ap.add_argument("--sample",       required=True, help="sample name in pangenie VCF (1001G ID)")
    ap.add_argument("--truth-sample", required=True, help="sample name in truth VCF (cactus assembly ID)")
    ap.add_argument("--out",          required=True, help="output prefix")
    args = ap.parse_args()

    truth_idx = index_truth(args.truth_vcf, args.truth_sample)
    compare(args.pangenie_vcf, args.sample, truth_idx, args.out)


if __name__ == "__main__":
    main()

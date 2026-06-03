"""
grenenet_snp_concordance.py — per-ecotype SNP concordance: PanGenie genotype
VCF vs GrENE-Net 231-founder SNP VCF (independent ground truth).

Two key wrinkles compared to loo_concordance.py:

  1. GrENE-Net chrom names are '1' .. '5' (no 'Chr' prefix); cactus/PanGenie
     output uses 'Chr1' .. 'Chr5'. Add chrom-rename map.
  2. PanGenie output is diploid (e.g. 1/1 for inbred carriers), GrENE-Net
     is also diploid (e.g. 1/1, 0/1, 0/0). Both sides yield comparable doses
     under gt_to_dose; no diploidization needed.

Concordance is computed only at SNP records (PanGenie has SVs which GrENE-Net
doesn't carry). For each shared (chrom_translated, pos, ref, alt), compare
truth_dose vs pg_dose, then aggregate GC + nRD per ecotype.

Output:
  <out>_summary.tsv   single row: panel, ecotype, n_pg_snp, n_truth_snp,
                      n_overlap, n_concordant, n_nonref_either, n_nonref_disc,
                      GC, nRD, n_truth_missing, n_pg_missing
"""
import argparse, sys
import pysam
from collections import Counter

# PanGenie/cactus -> GrENE-Net chrom translation
CHROM_MAP = {f"Chr{i}": str(i) for i in range(1, 6)}

def gt_to_dose(gt):
    if gt is None or len(gt) == 0: return None
    vals = [g for g in gt if g is not None]
    if not vals: return None
    return sum(1 for g in vals if g > 0)

def is_snp(ref, alt):
    return len(ref) == 1 and len(alt) == 1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pangenie-vcf", required=True)
    ap.add_argument("--truth-vcf",    required=True, help="GrENE-Net VCF (bgzipped+tabix)")
    ap.add_argument("--sample",       required=True)
    ap.add_argument("--panel",        required=True, help="main or loo")
    ap.add_argument("--out",          required=True, help="output prefix")
    args = ap.parse_args()

    # Step 1: build truth dict for this ecotype only
    print(f"[truth] reading GrENE-Net VCF for sample={args.sample}", file=sys.stderr)
    tv = pysam.VariantFile(args.truth_vcf)
    if args.sample not in tv.header.samples:
        sys.exit(f"sample '{args.sample}' not in truth VCF")
    tv.subset_samples([args.sample])
    truth = {}     # (chrom_grenenet, pos, ref, alt) -> dose
    n_truth_snp = 0
    for rec in tv.fetch():
        if not rec.alts: continue
        for alt in rec.alts:
            if not is_snp(rec.ref, alt): continue
            n_truth_snp += 1
            d = gt_to_dose(rec.samples[args.sample]["GT"])
            truth[(rec.chrom, rec.pos, rec.ref, alt)] = d
    print(f"[truth] {n_truth_snp:,} truth SNP entries", file=sys.stderr)

    # Step 2: stream PanGenie VCF, compare at SNP records
    pv = pysam.VariantFile(args.pangenie_vcf)
    if args.sample not in pv.header.samples:
        sys.exit(f"sample '{args.sample}' not in pangenie VCF")
    pv.subset_samples([args.sample])

    n_pg_snp = 0
    n_overlap = 0
    n_concordant = 0
    n_nonref_either = 0
    n_nonref_disc = 0
    n_truth_missing = 0  # in truth dict but truth dose is None
    n_pg_missing = 0     # PanGenie GT is None
    n_no_truth_record = 0  # this PanGenie SNP isn't in GrENE-Net at all
    for rec in pv.fetch():
        if not rec.alts: continue
        chrom_g = CHROM_MAP.get(rec.chrom, rec.chrom)
        for alt in rec.alts:
            if not is_snp(rec.ref, alt): continue
            n_pg_snp += 1
            key = (chrom_g, rec.pos, rec.ref, alt)
            if key not in truth:
                n_no_truth_record += 1
                continue
            n_overlap += 1
            truth_dose = truth[key]
            pg_dose = gt_to_dose(rec.samples[args.sample]["GT"])
            if truth_dose is None: n_truth_missing += 1
            if pg_dose is None:    n_pg_missing += 1
            if truth_dose is None or pg_dose is None: continue
            concord = (truth_dose == pg_dose)
            if concord: n_concordant += 1
            if truth_dose > 0 or pg_dose > 0:
                n_nonref_either += 1
                if not concord: n_nonref_disc += 1

    matched = n_overlap - n_truth_missing - n_pg_missing
    GC  = n_concordant / matched if matched > 0 else float("nan")
    nRD = n_nonref_disc / n_nonref_either if n_nonref_either > 0 else float("nan")

    out = open(args.out + "_summary.tsv", "w")
    out.write("panel\tecotype\tn_pg_snp\tn_truth_snp\tn_overlap\tn_no_truth_record\t"
              "n_truth_missing\tn_pg_missing\tn_matched\tn_concordant\tn_nonref_either\tn_nonref_disc\tGC\tnRD\n")
    out.write(f"{args.panel}\t{args.sample}\t{n_pg_snp}\t{n_truth_snp}\t{n_overlap}\t{n_no_truth_record}\t"
              f"{n_truth_missing}\t{n_pg_missing}\t{matched}\t{n_concordant}\t{n_nonref_either}\t{n_nonref_disc}\t{GC:.6f}\t{nRD:.6f}\n")
    out.close()
    print(f"\n=== {args.sample} ({args.panel}) ===")
    print(f"  PanGenie SNPs:     {n_pg_snp:,}")
    print(f"  GrENE-Net SNPs:    {n_truth_snp:,}")
    print(f"  PG SNPs in truth:  {n_overlap:,} ({n_overlap/max(n_pg_snp,1)*100:.1f}%)")
    print(f"  matched (both call): {matched:,}")
    print(f"  GC:    {GC:.4f}  ({n_concordant:,}/{matched:,})")
    print(f"  nRD:   {nRD:.4f}  ({n_nonref_disc:,}/{n_nonref_either:,})")
    print(f"  wrote {args.out}_summary.tsv")

if __name__ == "__main__":
    main()

"""
merged_vcf_stats.py — comprehensive stats + integrity check on the merged
231-founder VCF (output of merge_vcfs.sh).

Computes:
  - integrity:          n_samples, n_records, gz/tabix sanity
  - per-record:         AC (carriers) and AN (called alleles), missingness
  - per-sample:         n_records_called, n_records_with_alt, n_missing
  - by size class:      SNP / small_indel / small_sv / medium_sv / large_sv
                        breakdown of carriers and missingness
  - per-cohort:         cactus-80 vs PanGenie-151 — separate accounting since
                        they use different ploidy and different evidence sources

For the 151 PanGenie founders specifically:
  - "How many SVs were called for the new ecotypes?"
    = per-ecotype count of records where GT has any non-zero allele,
      restricted to size_class in (small_sv, medium_sv, large_sv)
  - "How many missing that need Beagle imputation?"
    = per-record count of `./. ` cells across the 231 panel
"""
import argparse
from collections import Counter, defaultdict
import pysam

CACTUS_80_PATH = "/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work/sample_rename.txt"

def size_class(ref_len, alt_len):
    if ref_len == 1 and alt_len == 1: return "SNP"
    m = max(ref_len, alt_len)
    if m < 50:    return "small_indel"
    if m < 500:   return "small_sv"
    if m < 5000:  return "medium_sv"
    return "large_sv"

def has_alt(gt):
    if gt is None or len(gt) == 0: return None  # missing
    vals = [g for g in gt if g is not None]
    if not vals: return None
    return any(g > 0 for g in vals)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--max-records", type=int, default=0, help="0 = all (debug)")
    args = ap.parse_args()

    # Cohort assignment: 80 cactus assembly-overlap vs 151 PanGenie
    cactus_set = set()
    with open(CACTUS_80_PATH) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2: cactus_set.add(fields[1])  # 1001G ID

    v = pysam.VariantFile(args.vcf)
    samples = list(v.header.samples)
    F = len(samples)
    cohort = ["cactus" if s in cactus_set else "pangenie" for s in samples]
    n_cactus = sum(1 for c in cohort if c == "cactus")
    n_pang   = sum(1 for c in cohort if c == "pangenie")
    print(f"== panel ==  n_samples={F}  cactus={n_cactus}  pangenie={n_pang}")

    # Per-sample tallies
    sample_called  = [0] * F
    sample_alt     = [0] * F
    sample_missing = [0] * F
    sample_alt_by_class = defaultdict(lambda: [0]*F)  # cls -> [count_per_sample]

    # Per-record tallies (we don't need to store per record; aggregate by class)
    cls_record_count   = Counter()
    cls_total_called   = Counter()  # sum of n_called over records
    cls_total_alt      = Counter()  # sum of n_alt-carrier over records
    cls_total_missing  = Counter()
    cls_records_with_any_missing = Counter()
    cls_records_with_all_called  = Counter()
    # For "imputation candidates": records where missing fraction > threshold
    cls_records_missing_gt5pct  = Counter()
    cls_records_missing_gt10pct = Counter()

    n_records = 0
    for rec in v.fetch():
        if not rec.alts: continue
        # Use longest ALT for size classification (multi-allelic conservative)
        alt_max = max(len(a) for a in rec.alts)
        cls = size_class(len(rec.ref), alt_max)
        cls_record_count[cls] += 1
        n_records += 1

        n_alt = 0; n_called = 0; n_missing = 0
        for f_idx, sname in enumerate(samples):
            ha = has_alt(rec.samples[sname]["GT"])
            if ha is None:
                sample_missing[f_idx] += 1
                n_missing += 1
            else:
                sample_called[f_idx] += 1
                n_called += 1
                if ha:
                    sample_alt[f_idx] += 1
                    sample_alt_by_class[cls][f_idx] += 1
                    n_alt += 1
        cls_total_called[cls]  += n_called
        cls_total_alt[cls]     += n_alt
        cls_total_missing[cls] += n_missing
        if n_missing == 0: cls_records_with_all_called[cls] += 1
        else:              cls_records_with_any_missing[cls] += 1
        if n_missing / F > 0.05:  cls_records_missing_gt5pct[cls]  += 1
        if n_missing / F > 0.10:  cls_records_missing_gt10pct[cls] += 1

        if args.max_records and n_records >= args.max_records: break

    # Write per-class summary
    classes = ["SNP","small_indel","small_sv","medium_sv","large_sv"]
    with open(args.out + "_by_size_class.tsv", "w") as f:
        f.write("size_class\tn_records\tn_called_total\tn_alt_total\tn_missing_total\t"
                "records_all_called\trecords_any_missing\trecords_missing_gt5pct\trecords_missing_gt10pct\t"
                "frac_called\tfrac_alt_of_called\tfrac_missing\n")
        for c in classes:
            nr = cls_record_count.get(c, 0)
            ncalled = cls_total_called.get(c, 0)
            nalt = cls_total_alt.get(c, 0)
            nmiss = cls_total_missing.get(c, 0)
            f.write(f"{c}\t{nr}\t{ncalled}\t{nalt}\t{nmiss}\t"
                    f"{cls_records_with_all_called.get(c,0)}\t{cls_records_with_any_missing.get(c,0)}\t"
                    f"{cls_records_missing_gt5pct.get(c,0)}\t{cls_records_missing_gt10pct.get(c,0)}\t"
                    f"{ncalled/(nr*F) if nr>0 else 0:.6f}\t"
                    f"{nalt/ncalled if ncalled>0 else 0:.6f}\t"
                    f"{nmiss/(nr*F) if nr>0 else 0:.6f}\n")

    # Per-sample TSV
    with open(args.out + "_by_sample.tsv", "w") as f:
        f.write("sample\tcohort\tn_called\tn_alt\tn_missing\tn_SNP_alt\tn_small_indel_alt\t"
                "n_small_sv_alt\tn_medium_sv_alt\tn_large_sv_alt\tn_sv_alt_total\n")
        for f_idx, sname in enumerate(samples):
            sv_total = sum(sample_alt_by_class.get(c, [0]*F)[f_idx] for c in ('small_sv','medium_sv','large_sv'))
            f.write(f"{sname}\t{cohort[f_idx]}\t{sample_called[f_idx]}\t{sample_alt[f_idx]}\t{sample_missing[f_idx]}\t"
                    f"{sample_alt_by_class['SNP'][f_idx]}\t{sample_alt_by_class['small_indel'][f_idx]}\t"
                    f"{sample_alt_by_class['small_sv'][f_idx]}\t{sample_alt_by_class['medium_sv'][f_idx]}\t"
                    f"{sample_alt_by_class['large_sv'][f_idx]}\t{sv_total}\n")

    # Headline summary to stdout
    print(f"\n== records ==  total={n_records:,}")
    for c in classes:
        nr = cls_record_count.get(c, 0)
        nmiss = cls_total_missing.get(c, 0)
        print(f"  {c:<13s} {nr:>10,}  any_missing={cls_records_with_any_missing.get(c,0):>9,}"
              f"  >5%_missing={cls_records_missing_gt5pct.get(c,0):>9,}  >10%_missing={cls_records_missing_gt10pct.get(c,0):>9,}")

    print(f"\n== headlining: how many SV calls per cohort ==")
    pang_samples = [(i, samples[i]) for i in range(F) if cohort[i] == "pangenie"]
    cact_samples = [(i, samples[i]) for i in range(F) if cohort[i] == "cactus"]
    for label, samps in [("PanGenie 151", pang_samples), ("cactus 80", cact_samples)]:
        sv_per_sample = [sum(sample_alt_by_class[c][i] for c in ('small_sv','medium_sv','large_sv')) for i,_ in samps]
        if sv_per_sample:
            print(f"  {label}: SV-alt per sample  median={int(sorted(sv_per_sample)[len(sv_per_sample)//2]):,}  "
                  f"min={min(sv_per_sample):,}  max={max(sv_per_sample):,}  total={sum(sv_per_sample):,}")

    print(f"\n  wrote {args.out}_by_size_class.tsv")
    print(f"  wrote {args.out}_by_sample.tsv")

if __name__ == "__main__":
    main()

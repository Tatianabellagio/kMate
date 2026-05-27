"""
qc_genotyped_one.py — per-sample QC of one PanGenie-genotyped VCF.

Reports:
  - n_records, n_called, n_missing
  - GT category breakdown: hom_ref, het, hom_alt, compound_het, missing
  - GQ summary: median, mean, frac >=30, >=100, ==10000
  - Records by SV size class (using REF/ALT lengths)

Output: one-line TSV row appended to a master file.
"""
import argparse, sys
from collections import Counter
import pysam

def classify_gt(gt):
    if gt is None or len(gt) == 0: return "missing"
    vals = [g for g in gt if g is not None]
    if not vals: return "missing"
    n_alt = sum(1 for g in vals if g > 0)
    if n_alt == 0: return "hom_ref"
    if n_alt == len(vals):
        # all non-zero — could be hom_alt or compound (different alts)
        return "hom_alt" if len(set(vals)) == 1 else "compound_het"
    return "het"

def size_class(ref_len, alt_len):
    if ref_len == 1 and alt_len == 1: return "SNP"
    m = max(ref_len, alt_len)
    if m < 50:    return "small_indel"
    if m < 500:   return "small_sv"
    if m < 5000:  return "medium_sv"
    return "large_sv"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    vf = pysam.VariantFile(args.vcf)
    if args.sample not in vf.header.samples:
        sys.exit(f"sample '{args.sample}' not in VCF")

    gt_counts = Counter()
    sc_counts = Counter()  # (size_class, gt_category) -> count
    gq_vals = []
    n_records = 0
    n_alts_total = 0
    sample_idx = list(vf.header.samples).index(args.sample)

    for rec in vf.fetch():
        if not rec.alts: continue
        n_records += 1
        s = rec.samples[args.sample]
        gt = s.get("GT")
        gq = s.get("GQ")
        cat = classify_gt(gt)
        gt_counts[cat] += 1
        if gq is not None and gq != ".":
            try:
                gq_vals.append(int(gq))
            except (TypeError, ValueError):
                pass
        # By-size-class breakdown — use the single longest ALT
        cls = size_class(len(rec.ref), max(len(a) for a in rec.alts))
        sc_counts[(cls, cat)] += 1

    n_total = sum(gt_counts.values())
    n_missing = gt_counts.get("missing", 0)
    n_called = n_total - n_missing
    gq_sorted = sorted(gq_vals)
    median_gq = gq_sorted[len(gq_sorted)//2] if gq_sorted else 0
    mean_gq = sum(gq_vals)/len(gq_vals) if gq_vals else 0

    out = open(args.out, "w")
    out.write("sample\tn_records\tn_called\tn_missing\thom_ref\thet\thom_alt\tcompound_het\t"
              "median_gq\tmean_gq\tgq_ge30_pct\tgq_ge100_pct\tgq_max_n\t"
              + "\t".join(f"{c}_{g}" for c in ['SNP','small_indel','small_sv','medium_sv','large_sv']
                          for g in ['hom_ref','het','hom_alt','compound_het','missing']) + "\n")
    out.write(f"{args.sample}\t{n_total}\t{n_called}\t{n_missing}\t"
              f"{gt_counts.get('hom_ref',0)}\t{gt_counts.get('het',0)}\t{gt_counts.get('hom_alt',0)}\t{gt_counts.get('compound_het',0)}\t"
              f"{median_gq}\t{mean_gq:.1f}\t"
              f"{sum(1 for g in gq_vals if g>=30)/max(len(gq_vals),1)*100:.2f}\t"
              f"{sum(1 for g in gq_vals if g>=100)/max(len(gq_vals),1)*100:.2f}\t"
              f"{sum(1 for g in gq_vals if g>=10000)}\t")
    cells = []
    for c in ['SNP','small_indel','small_sv','medium_sv','large_sv']:
        for g in ['hom_ref','het','hom_alt','compound_het','missing']:
            cells.append(str(sc_counts.get((c,g), 0)))
    out.write("\t".join(cells) + "\n")
    out.close()
    print(f"wrote {args.out}: {n_total} records, GT={dict(gt_counts)}")

if __name__ == "__main__":
    main()

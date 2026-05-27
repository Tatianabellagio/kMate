"""
merged_vcf_stats_aggregate.py — read bcftools-query stream from stdin, fast.

Input format (one record per line):
  CHROM  POS  REF  ALT[,...]  GT_sample1  GT_sample2 ...
GT format: '0/0', '0/1', '1/1', '0', '1', '.', './.', etc.

Outputs:
  <out>_by_size_class.tsv  per (size_class) × cohort: n_records, AC, AN, missing,
                           records_any_missing, records_miss_gt5pct/gt10pct
  <out>_by_sample.tsv      per sample: cohort, n_called, n_alt, n_missing,
                           n_alt_by_size_class
  <out>_per_record.tsv.gz  per record: chrom, pos, size_class, ac/an/missing
                           for full + each cohort
"""
import sys, argparse, gzip, subprocess

CLASSES = ["SNP","small_indel","small_sv","medium_sv","large_sv"]

BCF = "/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools"
VCF = "/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/data/merged/founders_231_chr.vcf.gz"

def size_class(ref_len, alt_len):
    if ref_len == 1 and alt_len == 1: return "SNP"
    m = max(ref_len, alt_len)
    if m < 50:    return "small_indel"
    if m < 500:   return "small_sv"
    if m < 5000:  return "medium_sv"
    return "large_sv"

def parse_gt(gt_str):
    """Return (called: bool, alt: bool)."""
    if not gt_str or gt_str == '.' or gt_str == './.' or gt_str == '.|.':
        return (False, False)
    has_call = False; has_alt = False
    for a in gt_str.replace('|','/').split('/'):
        if a == '.':
            continue
        elif a == '0':
            has_call = True
        else:
            has_call = True
            has_alt = True
    return (has_call, has_alt)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cactus-samples", required=True)
    ap.add_argument("--pangenie-samples", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--vcf", default=VCF, help="VCF whose query results are read from stdin (used for sample-order header lookup)")
    args = ap.parse_args()

    # Sample order = VCF header order = bcftools query column order
    samples = subprocess.check_output([BCF, "query", "-l", args.vcf]).decode().rstrip("\n").split('\n')
    F = len(samples)
    cactus_set = set(open(args.cactus_samples).read().split())
    pang_set   = set(open(args.pangenie_samples).read().split())
    is_pang = [1 if samples[i] in pang_set else 0 for i in range(F)]
    print(f"[stats] F={F}  cactus={F-sum(is_pang)}  pangenie={sum(is_pang)}", file=sys.stderr)

    # Aggregates
    cls_agg = {c: {k: 0 for k in
        ['n_records','ac_all','an_all','missing_all','records_all_called','records_any_missing',
         'records_miss_gt5pct','records_miss_gt10pct','records_pang_any_missing','records_cact_any_missing',
         'ac_pang','an_pang','missing_pang','ac_cact','an_cact','missing_cact']} for c in CLASSES}
    n_called  = [0] * F
    n_alt     = [0] * F
    n_missing = [0] * F
    n_alt_by_class = {c: [0]*F for c in CLASSES}

    rec_path = args.out_prefix + "_per_record.tsv.gz"
    rec_f = gzip.open(rec_path, "wt")
    rec_f.write("chrom\tpos\tsize_class\tref_len\talt_max_len\tac_all\tan_all\tmissing_all\t"
                "ac_pang\tan_pang\tmissing_pang\tac_cact\tan_cact\tmissing_cact\n")

    n_records = 0
    for line in sys.stdin:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4 + F: continue
        chrom = fields[0]; pos = fields[1]
        ref_len = len(fields[2])
        alts = fields[3].split(",")
        alt_max_len = max(len(a) for a in alts)
        cls = size_class(ref_len, alt_max_len)

        ac_all = an_all = miss_all = 0
        ac_pang = an_pang = miss_pang = 0
        ac_cact = an_cact = miss_cact = 0
        for i, g in enumerate(fields[4:4+F]):
            called, is_alt = parse_gt(g)
            if called:
                an_all += 1
                if is_alt:
                    ac_all += 1
                    n_alt[i] += 1
                    n_alt_by_class[cls][i] += 1
                n_called[i] += 1
                if is_pang[i]:
                    an_pang += 1
                    if is_alt: ac_pang += 1
                else:
                    an_cact += 1
                    if is_alt: ac_cact += 1
            else:
                miss_all += 1
                n_missing[i] += 1
                if is_pang[i]: miss_pang += 1
                else:          miss_cact += 1

        a = cls_agg[cls]
        a['n_records']   += 1
        a['ac_all']      += ac_all
        a['an_all']      += an_all
        a['missing_all'] += miss_all
        a['ac_pang']     += ac_pang
        a['an_pang']     += an_pang
        a['missing_pang']+= miss_pang
        a['ac_cact']     += ac_cact
        a['an_cact']     += an_cact
        a['missing_cact']+= miss_cact
        if miss_all == 0: a['records_all_called'] += 1
        if miss_all > 0:  a['records_any_missing'] += 1
        if miss_all / F > 0.05:  a['records_miss_gt5pct']  += 1
        if miss_all / F > 0.10:  a['records_miss_gt10pct'] += 1
        if miss_pang > 0: a['records_pang_any_missing'] += 1
        if miss_cact > 0: a['records_cact_any_missing'] += 1

        rec_f.write(f"{chrom}\t{pos}\t{cls}\t{ref_len}\t{alt_max_len}\t"
                    f"{ac_all}\t{an_all}\t{miss_all}\t{ac_pang}\t{an_pang}\t{miss_pang}\t"
                    f"{ac_cact}\t{an_cact}\t{miss_cact}\n")

        n_records += 1
        if n_records % 500000 == 0:
            print(f"  ... {n_records:,} records", file=sys.stderr)
    rec_f.close()

    cls_path = args.out_prefix + "_by_size_class.tsv"
    with open(cls_path, "w") as f:
        f.write("size_class\tn_records\tac_all\tan_all\tmissing_all\trecords_all_called\trecords_any_missing\t"
                "records_miss_gt5pct\trecords_miss_gt10pct\trecords_pang_any_missing\trecords_cact_any_missing\t"
                "ac_pang\tan_pang\tmissing_pang\tac_cact\tan_cact\tmissing_cact\n")
        for c in CLASSES:
            a = cls_agg[c]
            f.write("\t".join([c] + [str(a[k]) for k in
                ['n_records','ac_all','an_all','missing_all','records_all_called','records_any_missing',
                 'records_miss_gt5pct','records_miss_gt10pct','records_pang_any_missing','records_cact_any_missing',
                 'ac_pang','an_pang','missing_pang','ac_cact','an_cact','missing_cact']]) + "\n")

    samp_path = args.out_prefix + "_by_sample.tsv"
    with open(samp_path, "w") as f:
        f.write("sample\tcohort\tn_called\tn_alt\tn_missing\tn_SNP_alt\tn_small_indel_alt\t"
                "n_small_sv_alt\tn_medium_sv_alt\tn_large_sv_alt\tn_sv_alt_total\n")
        for i, s in enumerate(samples):
            cohort = "pangenie" if is_pang[i] else "cactus"
            sv_total = n_alt_by_class['small_sv'][i] + n_alt_by_class['medium_sv'][i] + n_alt_by_class['large_sv'][i]
            f.write(f"{s}\t{cohort}\t{n_called[i]}\t{n_alt[i]}\t{n_missing[i]}\t"
                    f"{n_alt_by_class['SNP'][i]}\t{n_alt_by_class['small_indel'][i]}\t"
                    f"{n_alt_by_class['small_sv'][i]}\t{n_alt_by_class['medium_sv'][i]}\t"
                    f"{n_alt_by_class['large_sv'][i]}\t{sv_total}\n")

    print(f"[stats] DONE  n_records={n_records:,}", file=sys.stderr)
    print(f"  wrote {cls_path}", file=sys.stderr)
    print(f"  wrote {samp_path}", file=sys.stderr)
    print(f"  wrote {rec_path}", file=sys.stderr)

if __name__ == "__main__":
    main()

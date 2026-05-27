"""
compute_ld.py — pairwise LD (r²) sampler for one VCF, one chromosome.

Reads VCF via cyvcf2, extracts per-marker dose vectors over all samples,
samples N anchor markers, computes r² with neighbors within max_dist bp,
emits TSV of (pos1, pos2, distance, r2, class1, class2).

Handles mixed ploidy (haploid-from-cactus + diploid-from-PanGenie) by
mapping each haploid allele to {0, 2} for dose-equivalent comparison with
diploid 0/1/2.

Marker class labeled per record:
  SNP          ref_len == 1 and alt_len == 1
  small_indel  max(ref,alt) < 50
  small_sv     50 <= max < 500
  medium_sv    500 <= max < 5000
  large_sv     max >= 5000

Filter: retain only biallelic records with maf >= --min-maf and
non-missing fraction >= 0.7. Multi-allelic records are SKIPPED to keep
dose computation simple (LD on multi-allelic requires per-alt handling).
"""
import argparse, sys, time
import numpy as np
import cyvcf2

def size_class(ref_len, alt_len):
    if ref_len == 1 and alt_len == 1: return "SNP"
    m = max(ref_len, alt_len)
    if m < 50:    return "small_indel"
    if m < 500:   return "small_sv"
    if m < 5000:  return "medium_sv"
    return "large_sv"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf",      required=True)
    ap.add_argument("--chrom",    required=True, help="e.g. Chr1 or 1")
    ap.add_argument("--out",      required=True)
    ap.add_argument("--n-anchors", type=int, default=10000,
        help="number of anchor markers to use (random sample)")
    ap.add_argument("--max-dist",  type=int, default=1_000_000,
        help="neighbor distance window in bp (default 1Mb)")
    ap.add_argument("--min-maf",   type=float, default=0.05)
    ap.add_argument("--min-call-frac", type=float, default=0.7)
    ap.add_argument("--seed",      type=int, default=42)
    ap.add_argument("--anchor-classes", default="SNP,small_indel,small_sv,medium_sv,large_sv",
        help="comma list of size_classes eligible to be sampled as anchors")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    anchor_classes = set(args.anchor_classes.split(","))

    # Pass 1: collect marker metadata + dose vectors
    print(f"[load] reading {args.vcf} chrom={args.chrom}", file=sys.stderr)
    t0 = time.time()
    v = cyvcf2.VCF(args.vcf)
    pos_arr = []
    cls_arr = []
    dose_rows = []
    n_total = 0; n_kept = 0
    for rec in v(args.chrom):
        n_total += 1
        if len(rec.ALT) != 1: continue        # multi-allelic skipped
        ref = rec.REF; alt = rec.ALT[0]
        cls = size_class(len(ref), len(alt))
        # gt_types: 0=hom_ref, 1=het, 2=unknown, 3=hom_alt for cyvcf2
        # gt_alt_freqs gives 0..1; or use rec.genotypes which is list of [a0,a1,phased]
        # Use rec.genotype.array() for fast access — shape (n_samples, 3)
        gts = rec.genotype.array()
        # For each sample, dose = sum of non-zero alleles. Mixed ploidy:
        #   haploid sample has gts[s] = [a0, -1, phased] in cyvcf2 → only first allele real
        # Simpler: use gt_types
        # gt_types codes: 0 hom_ref, 1 het, 2 unknown, 3 hom_alt — but for haploid,
        # cyvcf2 maps `0` -> 0 (hom_ref), `1` -> 3 (hom_alt). UNKNOWN_PLOIDY=2.
        gtt = np.asarray(rec.gt_types)
        dose = np.full(gtt.shape[0], -1, dtype=np.int8)  # -1 = missing
        dose[gtt == 0] = 0  # hom_ref
        dose[gtt == 1] = 1  # het
        dose[gtt == 3] = 2  # hom_alt
        # Filter low call rate
        n_called = (dose >= 0).sum()
        call_frac = n_called / dose.shape[0]
        if call_frac < args.min_call_frac: continue
        # Compute MAF on non-missing samples
        called_dose = dose[dose >= 0]
        af = called_dose.sum() / (2 * called_dose.shape[0])
        maf = min(af, 1 - af)
        if maf < args.min_maf: continue
        pos_arr.append(rec.POS)
        cls_arr.append(cls)
        dose_rows.append(dose)
        n_kept += 1
        if n_kept % 50000 == 0:
            print(f"  ... {n_kept:,} kept (of {n_total:,} seen) in {time.time()-t0:.0f}s", file=sys.stderr)
    print(f"[load] done: {n_kept:,} markers kept (of {n_total:,} total) in {time.time()-t0:.0f}s", file=sys.stderr)
    if n_kept == 0:
        sys.exit("no markers passed filters")

    pos = np.asarray(pos_arr, dtype=np.int64)
    cls = np.asarray(cls_arr)
    dose = np.stack(dose_rows)  # shape (n_kept, n_samples)
    n_samples = dose.shape[1]
    print(f"  dose matrix: {dose.shape}  n_samples={n_samples}", file=sys.stderr)

    # Sample anchor markers from eligible classes
    eligible = np.where(np.isin(cls, list(anchor_classes)))[0]
    if len(eligible) == 0:
        sys.exit(f"no markers with class in {anchor_classes}")
    n_anchors = min(args.n_anchors, len(eligible))
    anchor_idx = rng.choice(eligible, size=n_anchors, replace=False)
    anchor_idx.sort()
    print(f"  sampled {n_anchors} anchors from {len(eligible)} eligible", file=sys.stderr)

    # Pass 2: per anchor, compute r² with all neighbors within max_dist
    out = open(args.out, "w")
    out.write("pos1\tpos2\tdistance\tclass1\tclass2\tn_pairs_called\tr2\n")
    n_written = 0
    t1 = time.time()
    # Pre-compute mean and std (treating -1 as missing → mask per pair)
    for ai_idx, ai in enumerate(anchor_idx):
        a_pos = pos[ai]
        a_dose = dose[ai].astype(np.float32)
        a_mask = a_dose >= 0
        # Find neighbors: positions in [a_pos+1, a_pos+max_dist]
        # Search bounds via np.searchsorted on pos array
        lo = ai + 1
        hi_idx = np.searchsorted(pos, a_pos + args.max_dist, side='right')
        if lo >= hi_idx: continue
        # Loop neighbors
        for ni in range(lo, hi_idx):
            n_pos = pos[ni]
            d = n_pos - a_pos
            n_dose = dose[ni].astype(np.float32)
            n_mask = n_dose >= 0
            joint = a_mask & n_mask
            n_pair = int(joint.sum())
            if n_pair < int(0.5 * n_samples): continue
            x = a_dose[joint]; y = n_dose[joint]
            xs = x.std(); ys = y.std()
            if xs == 0 or ys == 0: continue
            xm = x - x.mean(); ym = y - y.mean()
            r = (xm * ym).sum() / (np.sqrt((xm*xm).sum() * (ym*ym).sum()))
            r2 = r * r
            out.write(f"{a_pos}\t{n_pos}\t{d}\t{cls[ai]}\t{cls[ni]}\t{n_pair}\t{r2:.6f}\n")
            n_written += 1
        if (ai_idx + 1) % 500 == 0:
            print(f"  anchors {ai_idx+1}/{n_anchors}  pairs written {n_written:,}  "
                  f"elapsed {time.time()-t1:.0f}s", file=sys.stderr)
    out.close()
    print(f"[done] {n_written:,} pairs written -> {args.out} (total {time.time()-t0:.0f}s)", file=sys.stderr)

if __name__ == "__main__":
    main()

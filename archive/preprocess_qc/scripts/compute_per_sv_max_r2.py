"""
compute_per_sv_max_r2.py — for each SV in the production 231-founder VCF
(`founders_231_chr.vcf.gz`: 80 cactus + 151 PanGenie, no Beagle imputation),
find the max r² to any SNP within ±--window-bp.

Two modes:
  --mode within  → SNP neighbors come from the same production VCF
  --mode cross   → SNP neighbors come from the xwu GrENE-Net SNP-only catalog
                   (sample-aligned by founder name; chrom name translated
                   between Chr1 ↔ 1 etc.)

Multi-allelic decomposition: each ALT becomes its own biallelic record.
Per-founder dose: 1 if any allele in the call equals this ALT index, 0 if
called and not, -1 if missing. This treats every founder as 'carrier vs not',
which matches the cactus_em / cn_var consumer and avoids the
haploid-cactus-vs-diploid-PanGenie ploidy mismatch in r² computation.

Filtering:
  - drop ploidy-aware singletons (AC ≤ 1 or AC ≥ N-1) at scan time
  - require call rate ≥ --min-call-frac across the founder axis
  - MAF ≥ --min-maf

Output (one row per biallelic SV record passing filters):
  chrom, sv_pos, sv_class, sv_size, sv_ac, max_r2,
  best_tag_snp_pos, best_tag_snp_dist_bp, n_snps_in_window
"""
import argparse, sys, time
import numpy as np
import cyvcf2


SV_CLASSES = {"small_sv", "medium_sv", "large_sv"}


def size_class(ref_len, alt_len):
    if ref_len == 1 and alt_len == 1:
        return "SNP"
    m = max(ref_len, alt_len)
    if m < 50:
        return "small_indel"
    if m < 500:
        return "small_sv"
    if m < 5000:
        return "medium_sv"
    return "large_sv"


def per_alt_carrier_dose(rec, n_samples):
    """Return list of (alt_idx, dose) for every ALT allele in this record.
    `dose` is int8 of shape (n_samples,), values in {-1, 0, 1}, where 1 means
    the founder genotype contains this ALT index, 0 means the founder is
    called and does NOT contain this ALT, -1 means missing.

    Handles mixed ploidy uniformly: cyvcf2.genotype.array() returns one row
    per sample with [a0, a1, phased]. For haploid samples cyvcf2 fills a1 with
    -1 (i.e. alt position not present), so the carrier check (a0==k OR a1==k)
    still resolves correctly without misclassifying the haploid filler.
    """
    g = rec.genotype.array()           # (n_samples, 3) int — alleles + phased
    a0 = g[:, 0]
    a1 = g[:, 1]
    out = []
    for k in range(1, len(rec.ALT) + 1):
        carrier = ((a0 == k) | (a1 == k)).astype(np.int8)
        # Missing if BOTH allele slots are -1 (true ./.) OR if the only-allele
        # slot (haploid) is -1. cyvcf2 codes haploid `0` as a0=0, a1=-1 (so a1
        # is filler, not missing). True missing haploid is a0=-1, a1=-1.
        missing = (a0 < 0) & (a1 < 0)
        # For diploid `./.`, both are -1 too — same condition.
        # For haploid `.` (cactus side), cyvcf2 also yields a0=-1 a1=-1.
        dose = carrier
        dose[missing] = -1
        out.append((k, dose))
    return out


def stream_vcf_decomposed(vcf_path, target_chrom, founder_order,
                          min_maf, min_call_frac, sv_only=False, snp_only=False):
    """Stream a VCF, decompose multi-allelics per ALT, return arrays sorted by
    POS. founder_order may be None to skip reordering.

    Returns a dict with keys: pos, cls, size, ac, dose (n_records × n_founders).
    Records that fail call-rate or MAF filters are dropped.
    """
    v = cyvcf2.VCF(vcf_path)
    vcf_samples = [s for s in v.samples if s.strip()]
    if founder_order is not None:
        name_to_idx = {s: i for i, s in enumerate(vcf_samples)}
        missing = [f for f in founder_order if f not in name_to_idx]
        if missing:
            sys.exit(f"[fatal] {len(missing)} founders not in VCF (e.g. {missing[:3]})")
        src_idx = np.asarray([name_to_idx[f] for f in founder_order], dtype=np.int64)
        n_samples = len(founder_order)
    else:
        src_idx = None
        n_samples = len(vcf_samples)

    pos_l, cls_l, size_l, ac_l, dose_l = [], [], [], [], []
    n_recs = 0; n_alts = 0; t0 = time.time()
    for rec in v(target_chrom):
        n_recs += 1
        ref = rec.REF
        for k, dose in per_alt_carrier_dose(rec, n_samples=len(vcf_samples)):
            n_alts += 1
            alt = rec.ALT[k - 1]
            cls = size_class(len(ref), len(alt))
            if sv_only and cls not in SV_CLASSES:
                continue
            if snp_only and cls != "SNP":
                continue
            if src_idx is not None:
                dose = dose[src_idx]
            n_called = int((dose >= 0).sum())
            if n_called == 0 or n_called / n_samples < min_call_frac:
                continue
            ac = int((dose == 1).sum())
            # haploid-style symmetric MAF (carrier-AF, since each founder is
            # represented as carrier 0/1 after decomposition)
            af = ac / n_called
            maf = min(af, 1.0 - af)
            if maf < min_maf:
                continue
            # drop singletons on either edge
            if ac <= 1 or ac >= n_called - 1:
                continue
            pos_l.append(rec.POS)
            cls_l.append(cls)
            size_l.append(max(len(ref), len(alt)))
            ac_l.append(ac)
            dose_l.append(dose)
        if n_recs % 200_000 == 0:
            print(f"  ... {n_recs:,} recs / {n_alts:,} ALTs / {len(pos_l):,} kept "
                  f"in {time.time()-t0:.0f}s", file=sys.stderr)
    print(f"[stream] {target_chrom}: {n_recs:,} records, {n_alts:,} ALTs, "
          f"{len(pos_l):,} kept after filters in {time.time()-t0:.0f}s",
          file=sys.stderr)
    if not pos_l:
        return {"pos": np.empty(0, dtype=np.int64),
                "cls": np.empty(0, dtype="<U16"),
                "size": np.empty(0, dtype=np.int32),
                "ac":  np.empty(0, dtype=np.int32),
                "dose": np.empty((0, n_samples), dtype=np.int8)}
    pos = np.asarray(pos_l, dtype=np.int64)
    cls = np.asarray(cls_l)
    size = np.asarray(size_l, dtype=np.int32)
    ac = np.asarray(ac_l, dtype=np.int32)
    dose = np.stack(dose_l).astype(np.int8)
    order = np.argsort(pos, kind="stable")
    return {"pos": pos[order], "cls": cls[order], "size": size[order],
            "ac": ac[order], "dose": dose[order]}


def compute_max_r2(sv_pos, sv_dose, snp_pos, snp_dose,
                   window_bp, min_pair_frac):
    """Per-SV max r² scan. Both sv_dose and snp_dose are int8 carrier vectors
    in {-1, 0, 1}; -1 marks missing samples that are pairwise-excluded."""
    n_sv = sv_pos.shape[0]
    n_founders = sv_dose.shape[1]
    min_pair = int(min_pair_frac * n_founders)

    max_r2 = np.full(n_sv, np.nan, dtype=np.float32)
    best_pos = np.full(n_sv, -1, dtype=np.int64)
    best_dist = np.zeros(n_sv, dtype=np.int64)
    n_in_win = np.zeros(n_sv, dtype=np.int32)

    if snp_pos.size == 0:
        return max_r2, best_pos, best_dist, n_in_win

    snp_dose_f = snp_dose.astype(np.float32)
    snp_mask   = snp_dose >= 0
    sv_dose_f  = sv_dose.astype(np.float32)
    sv_mask    = sv_dose >= 0

    t0 = time.time()
    for i in range(n_sv):
        a_pos = sv_pos[i]
        lo = np.searchsorted(snp_pos, a_pos - window_bp, side="left")
        hi = np.searchsorted(snp_pos, a_pos + window_bp, side="right")
        if lo >= hi:
            continue
        a_full = sv_dose_f[i]
        am = sv_mask[i]
        cur_max = -1.0
        cur_pos = -1
        cur_dist = 0
        cnt = 0
        for j in range(lo, hi):
            n_pos = snp_pos[j]
            if n_pos == a_pos:
                continue
            mask_j = am & snp_mask[j]
            n_pair = int(mask_j.sum())
            if n_pair < min_pair:
                continue
            x = a_full[mask_j]
            y = snp_dose_f[j][mask_j]
            xs = x.std(); ys = y.std()
            if xs == 0 or ys == 0:
                continue
            xm = x - x.mean(); ym = y - y.mean()
            denom = np.sqrt((xm * xm).sum() * (ym * ym).sum())
            if denom == 0:
                continue
            r = float((xm * ym).sum() / denom)
            r2 = r * r
            cnt += 1
            if r2 > cur_max:
                cur_max = r2
                cur_pos = int(n_pos)
                cur_dist = int(n_pos - a_pos)
        n_in_win[i] = cnt
        if cur_max >= 0:
            max_r2[i] = cur_max
            best_pos[i] = cur_pos
            best_dist[i] = cur_dist
        if (i + 1) % 2000 == 0:
            print(f"  ... SV {i+1:,}/{n_sv:,}  elapsed {time.time()-t0:.0f}s",
                  file=sys.stderr)
    return max_r2, best_pos, best_dist, n_in_win


def get_founder_order(vcf_path):
    v = cyvcf2.VCF(vcf_path)
    return [s for s in v.samples if s.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prod-vcf", required=True,
                    help="Production 231-founder VCF (founders_231_chr.vcf.gz)")
    ap.add_argument("--mode", choices=["within", "cross"], required=True)
    ap.add_argument("--snp-vcf",
                    help="GrENE-Net SNP-only VCF (cross mode); ignored otherwise")
    ap.add_argument("--prod-chrom", required=True, help="e.g. Chr1")
    ap.add_argument("--snp-chrom",  help="cross-mode chrom label, e.g. 1")
    ap.add_argument("--out", required=True)
    ap.add_argument("--window-bp",     type=int,   default=50_000)
    ap.add_argument("--min-maf",       type=float, default=0.05)
    ap.add_argument("--min-call-frac", type=float, default=0.7)
    ap.add_argument("--min-pair-frac", type=float, default=0.5)
    args = ap.parse_args()

    if args.mode == "cross" and (not args.snp_vcf or not args.snp_chrom):
        sys.exit("--snp-vcf and --snp-chrom required for --mode cross")

    # Anchor founder order is the production VCF's sample order.
    founder_order = get_founder_order(args.prod_vcf)
    print(f"[init] founder_order: {len(founder_order)} samples", file=sys.stderr)

    # SVs from production VCF
    print(f"[stream-prod-svs] {args.prod_vcf} chrom={args.prod_chrom}", file=sys.stderr)
    sv = stream_vcf_decomposed(
        args.prod_vcf, args.prod_chrom, founder_order=None,
        min_maf=args.min_maf, min_call_frac=args.min_call_frac,
        sv_only=True, snp_only=False,
    )
    print(f"[svs] {sv['pos'].size:,} SVs after MAF + singleton filter "
          f"(class breakdown {dict((c, int((sv['cls']==c).sum())) for c in SV_CLASSES)})",
          file=sys.stderr)
    if sv["pos"].size == 0:
        sys.exit("no SVs to process")

    # SNP neighbors
    if args.mode == "within":
        print(f"[stream-prod-snps] {args.prod_vcf} chrom={args.prod_chrom}", file=sys.stderr)
        snp = stream_vcf_decomposed(
            args.prod_vcf, args.prod_chrom, founder_order=None,
            min_maf=args.min_maf, min_call_frac=args.min_call_frac,
            sv_only=False, snp_only=True,
        )
    else:
        print(f"[stream-grenenet-snps] {args.snp_vcf} chrom={args.snp_chrom}", file=sys.stderr)
        snp = stream_vcf_decomposed(
            args.snp_vcf, args.snp_chrom, founder_order=founder_order,
            min_maf=args.min_maf, min_call_frac=args.min_call_frac,
            sv_only=False, snp_only=True,
        )
    print(f"[snps] {snp['pos'].size:,} SNPs", file=sys.stderr)

    print(f"[scan] window=±{args.window_bp:,} bp", file=sys.stderr)
    max_r2, best_pos, best_dist, n_in_win = compute_max_r2(
        sv["pos"], sv["dose"], snp["pos"], snp["dose"],
        args.window_bp, args.min_pair_frac,
    )

    with open(args.out, "w") as f:
        f.write("chrom\tsv_pos\tsv_class\tsv_size\tsv_ac\tmax_r2\t"
                "best_tag_snp_pos\tbest_tag_snp_dist_bp\tn_snps_in_window\n")
        for i in range(sv["pos"].size):
            r2 = max_r2[i]
            r2_str = f"{r2:.6f}" if not np.isnan(r2) else "NA"
            f.write(f"{args.prod_chrom}\t{sv['pos'][i]}\t{sv['cls'][i]}\t"
                    f"{sv['size'][i]}\t{sv['ac'][i]}\t{r2_str}\t"
                    f"{best_pos[i] if best_pos[i] >= 0 else 'NA'}\t"
                    f"{best_dist[i]}\t{n_in_win[i]}\n")
    n_tagged = int(np.sum(~np.isnan(max_r2)))
    print(f"[done] wrote {sv['pos'].size:,} rows ({n_tagged:,} with at least one tag) -> {args.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()

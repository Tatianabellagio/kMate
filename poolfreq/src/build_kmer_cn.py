"""
Build the founder × k-mer copy-number matrix from PanGenie-index output + the
input multi-allelic VCF.

PanGenie merges nearby bubbles (within k-mer distance) into a single "merged
bubble". For each merged bubble at [start, end], we:
  1. Look up all VCF records in [start, end]
  2. For each founder, reconstruct their haplotype sequence across the region
     by applying their chosen alt for each VCF record
  3. Compute the canonical k-mer set of each founder's reconstructed sequence
  4. cn[f, k] = 1 iff k-mer k (from the bubble's unique-kmers list) is in
     founder f's k-mer set, else 0

Output is a sparse F × K matrix (CSR), plus arrays mapping each k-mer to its
bubble id and chromosome.
"""
from __future__ import annotations
import gzip
import numpy as np
from pathlib import Path
from scipy.sparse import csr_matrix, save_npz
import pysam


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT)[::-1]

def canonical(kmer: str) -> str:
    rc = revcomp(kmer)
    return kmer if kmer < rc else rc

def canonical_kmer_set(seq: str, k: int) -> set[str]:
    """Set of canonical k-mers in seq (excluding any with N)."""
    s = set()
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        km = seq[i:i+k]
        if "N" in km: continue
        s.add(canonical(km))
    return s


def reconstruct_haplotype(
    ref_seq: str,        # reference sequence over bubble region
    region_start: int,   # genomic start of ref_seq (1-based)
    vcf_records: list,   # list of (pos, ref, alts_tuple, founder_gt_idx)
    flank_left: str = "",
    flank_right: str = "",
) -> str:
    """Apply variants right-to-left to ref_seq, return the reconstructed sequence.

    vcf_records: list of tuples (pos, ref, [alt1, alt2, ...], gt_alt_index)
                 gt_alt_index: 0 = ref, 1 = first alt, etc.
                               -1 = missing (./.) — write N over the REF span
    """
    seq = ref_seq
    for pos, ref, alts, gt in sorted(vcf_records, key=lambda r: r[0], reverse=True):
        if gt == 0:
            continue
        offset = pos - region_start
        if offset < 0 or offset > len(seq):
            continue  # out of range — defensive
        rl = len(ref)
        if gt == -1:
            # Missing GT: mark the REF span with N. canonical_kmer_set skips
            # k-mers containing N, so the founder neither gets REF k-mer credit
            # nor ALT k-mer credit for k-mers spanning this position.
            seq = seq[:offset] + ("N" * rl) + seq[offset + rl:]
            continue
        if gt - 1 >= len(alts):
            continue  # missing alt index
        alt = alts[gt - 1]
        seq = seq[:offset] + alt + seq[offset + rl:]
    return flank_left + seq + flank_right


def parse_kmer_file(kmers_tsv_gz: str) -> list[dict]:
    """Parse a PanGenie kmers.tsv.gz file. Returns list of bubble dicts."""
    bubbles = []
    with gzip.open(kmers_tsv_gz, "rt") as f:
        next(f)  # header
        for line in f:
            p = line.rstrip("\n").split("\t")
            chrom, start, end = p[0], int(p[1]), int(p[2])
            kmers = p[3].split(",") if p[3] else []
            overhang = p[4] if len(p) > 4 else ""
            bubbles.append({
                "chrom": chrom, "start": start, "end": end,
                "kmers": kmers, "overhang": overhang,
            })
    return bubbles


def build_cn_for_chrom(
    kmers_tsv_gz: str,
    vcf_path: str,
    ref_fasta: str,
    chrom: str,
    k: int = 31,
    flank: int = 100,           # extra ref bases on each side for k-mer scanning
    max_bubbles: int | None = None,
    verbose: bool = True,
    treat_missing_as_n: bool = False,  # if True: ./. → N (no evidence). If False: ./. → REF (legacy v3 behavior).
):
    """Build the founder × k-mer matrix for one chromosome.

    Returns:
        cn:          scipy.sparse.csr_matrix of shape (F, K_total)
        kmer_index:  list of canonical k-mer strings (length K_total)
        bubble_id:   K_total-vector mapping each k-mer to its bubble idx (0..n_bubbles-1)
        bubble_meta: list of (chrom, start, end) per bubble (length n_bubbles)
        founders:    list of founder names (length F)
    """
    bubbles = parse_kmer_file(kmers_tsv_gz)
    bubbles = [b for b in bubbles if b["chrom"] == chrom]
    if max_bubbles is not None:
        bubbles = bubbles[:max_bubbles]
    if verbose:
        print(f"[build_cn] {chrom}: {len(bubbles):,} bubbles")

    # Open reference and VCF
    fasta = pysam.FastaFile(ref_fasta)
    vcf = pysam.VariantFile(vcf_path)
    founders = list(vcf.header.samples)
    F = len(founders)
    if verbose:
        print(f"[build_cn] founders: F={F}")

    # First pass: count total k-mers and build kmer-index per bubble
    bubble_kmer_canon: list[list[str]] = []
    total_K = 0
    for b in bubbles:
        canon = [canonical(km) for km in b["kmers"]]
        bubble_kmer_canon.append(canon)
        total_K += len(canon)
    if verbose:
        print(f"[build_cn] total k-mers across bubbles: {total_K:,}")

    # Pre-allocate cn as dense int8 (F × total_K). For Chr1: 82 × ~6M = ~500MB.
    # Use sparse builders if memory tight.
    cn_rows = []
    cn_cols = []
    bubble_id = np.zeros(total_K, dtype=np.int32)
    kmer_index = [""] * total_K

    k_offset = 0
    for b_idx, b in enumerate(bubbles):
        canon = bubble_kmer_canon[b_idx]
        K_b = len(canon)
        if K_b == 0:
            continue

        # Fetch ref sequence with flanks
        rs = max(0, b["start"] - flank)
        re_ = b["end"] + flank
        try:
            full_ref = fasta.fetch(chrom, rs, re_).upper()
        except Exception as e:
            if verbose:
                print(f"  WARN: fetch failed for {chrom}:{rs}-{re_}: {e}")
            k_offset += K_b
            continue

        # Find VCF records in [start, end] (1-based, inclusive)
        try:
            vcf_records = []
            for rec in vcf.fetch(chrom, b["start"] - 1, b["end"]):
                vcf_records.append(rec)
        except Exception:
            k_offset += K_b
            continue

        # Record kmer index + bubble id
        for j, km in enumerate(canon):
            bubble_id[k_offset + j] = b_idx
            kmer_index[k_offset + j] = km
        canon_set = set(canon)

        # For each founder, reconstruct haplotype and check kmer membership
        for f_idx, founder in enumerate(founders):
            # Build vcf_record_summary: (pos, ref, alts, gt_alt_index_for_this_founder)
            rec_summary = []
            for rec in vcf_records:
                gt = rec.samples[founder]["GT"]
                if gt is None or len(gt) == 0:
                    continue
                # Panel must be pre-haploidized. cn_full reconstructs ONE
                # sequence per founder; with diploid GTs (e.g. 0/1) we'd
                # silently use only gt[0] and disagree with cn_var's
                # any(a>0) carrier rule. Fail loudly instead.
                if len(gt) != 1:
                    raise ValueError(
                        f"Panel VCF must be haploid (got len(GT)={len(gt)} for "
                        f"{founder} at {rec.chrom}:{rec.pos}). Haploidize the "
                        f"VCF before building cn_full."
                    )
                if gt[0] is None:
                    # Missing GT. Behavior depends on treat_missing_as_n flag.
                    gt_idx = -1 if treat_missing_as_n else 0
                else:
                    gt_idx = gt[0]
                rec_summary.append((rec.pos, rec.ref, rec.alts, gt_idx))

            # Reconstruct: pad ref with left+right flanks already present in full_ref
            # Compute offset of bubble start within full_ref
            left_flank = full_ref[:b["start"] - 1 - rs]
            ref_in_bubble = full_ref[b["start"] - 1 - rs : b["end"] - rs]
            right_flank = full_ref[b["end"] - rs:]

            hap = reconstruct_haplotype(
                ref_in_bubble, b["start"], rec_summary,
                flank_left=left_flank, flank_right=right_flank,
            )
            hap_kmers = canonical_kmer_set(hap, k)

            # Check each canonical kmer
            for j, km in enumerate(canon):
                if km in hap_kmers:
                    cn_rows.append(f_idx)
                    cn_cols.append(k_offset + j)

        k_offset += K_b
        if verbose and (b_idx + 1) % 1000 == 0:
            print(f"  ... processed {b_idx+1:,}/{len(bubbles):,} bubbles, "
                  f"current cn nnz={len(cn_rows):,}")

    # Build sparse matrix
    cn = csr_matrix(
        (np.ones(len(cn_rows), dtype=np.int8), (cn_rows, cn_cols)),
        shape=(F, total_K),
        dtype=np.int8,
    )
    bubble_meta = [(b["chrom"], b["start"], b["end"]) for b in bubbles]
    if verbose:
        print(f"[build_cn] DONE. cn shape={cn.shape}, nnz={cn.nnz:,}, "
              f"density={cn.nnz / (cn.shape[0]*cn.shape[1]):.4%}")
    return cn, kmer_index, bubble_id, bubble_meta, founders


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmers", required=True, help="PanGenie kmers.tsv.gz")
    ap.add_argument("--vcf",   required=True, help="PanGenie input VCF (diploid, indexed)")
    ap.add_argument("--ref",   required=True, help="Reference FASTA")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out",   required=True, help="Output prefix")
    ap.add_argument("--max-bubbles", type=int, default=None)
    ap.add_argument("--treat-missing-as-n", action="store_true",
                    help="Write N at positions where founder has ./. GT, instead "
                         "of defaulting to REF. Drops k-mers spanning N from cn_full. "
                         "Use for v3qc-v2 and onwards. v3/v3qc were built without this.")
    args = ap.parse_args()

    if args.treat_missing_as_n:
        print("[build_cn] --treat-missing-as-n ON: ./. → N in haplotype")

    cn, kmer_index, bubble_id, bubble_meta, founders = build_cn_for_chrom(
        args.kmers, args.vcf, args.ref, args.chrom,
        max_bubbles=args.max_bubbles,
        treat_missing_as_n=args.treat_missing_as_n,
    )
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    save_npz(args.out + ".cn.npz", cn)
    np.savez(args.out + ".meta.npz",
             kmer_index=np.array(kmer_index),
             bubble_id=bubble_id,
             bubble_chrom=np.array([m[0] for m in bubble_meta]),
             bubble_start=np.array([m[1] for m in bubble_meta], dtype=np.int64),
             bubble_end=np.array([m[2] for m in bubble_meta], dtype=np.int64),
             founders=np.array(founders))
    print(f"Wrote {args.out}.cn.npz and {args.out}.meta.npz")

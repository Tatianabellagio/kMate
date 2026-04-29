"""
Build cn_var: founder × biallelic-VCF-record matrix.

For each row of the biallelic.norm VCF (one specific ALT of one bubble),
encode whether each founder carries that ALT in its phased haplotype.

cn_var[f, r] = 1 if founder f's GT at record r is the alt allele (any GT > 0),
               0 otherwise.

Used to project founder frequencies → per-record alt-allele frequencies:
    f_alt(r, sample i) = h_i · cn_var[:, r]

Output: scipy sparse CSR (F × N_records) + metadata (chrom, pos, ref_len, alt_len).
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import numpy as np
import pysam
from scipy.sparse import csr_matrix, save_npz


def build_cn_var(vcf_path: str, out_prefix: str, chrom_filter: str | None = None,
                 verbose: bool = True):
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    F = len(samples)
    print(f"Founders: F={F}")

    # First pass: count records (for sparse matrix sizing)
    chrom_arr, pos_arr, ref_len_arr, alt_len_arr = [], [], [], []
    rows, cols = [], []   # for sparse matrix
    n_kept = 0
    iterator = vcf.fetch(chrom_filter) if chrom_filter else vcf.fetch()

    t0 = time.time()
    for r_idx, rec in enumerate(iterator):
        # Each row in biallelic.norm VCF has exactly one ALT
        if not rec.alts:
            continue
        chrom_arr.append(rec.chrom)
        pos_arr.append(rec.pos)
        ref_len_arr.append(len(rec.ref))
        alt_len_arr.append(len(rec.alts[0]))

        for f_idx, sname in enumerate(samples):
            gt = rec.samples[sname]["GT"]
            if gt is None or len(gt) == 0:
                continue
            # GT is (a1, a2) for diploid; take any non-zero allele
            has_alt = any(a is not None and a > 0 for a in gt)
            if has_alt:
                rows.append(f_idx)
                cols.append(n_kept)
        n_kept += 1
        if verbose and n_kept % 100000 == 0:
            print(f"  {n_kept:,} records processed, nnz={len(rows):,}, "
                  f"elapsed={time.time()-t0:.0f}s")

    print(f"\nTotal records: {n_kept:,}")
    print(f"  nnz: {len(rows):,}, density: {len(rows)/(F*n_kept):.4f}")

    # Build sparse matrix
    cn_var = csr_matrix(
        (np.ones(len(rows), dtype=np.int8), (rows, cols)),
        shape=(F, n_kept),
        dtype=np.int8,
    )

    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    save_npz(out_prefix + ".cn_var.npz", cn_var)
    np.savez(
        out_prefix + ".meta.npz",
        founders=np.array(samples),
        chrom=np.array(chrom_arr),
        pos=np.array(pos_arr, dtype=np.int64),
        ref_len=np.array(ref_len_arr),
        alt_len=np.array(alt_len_arr),
    )
    print(f"Wrote {out_prefix}.cn_var.npz and {out_prefix}.meta.npz")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="biallelic.norm.vcf.gz path")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--chrom", default=None, help="optional: limit to one chrom")
    args = ap.parse_args()
    build_cn_var(args.vcf, args.out, args.chrom)

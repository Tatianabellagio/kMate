"""
Build var_pa: founder × biallelic-VCF-record matrices for AF projection.

Two parallel sparse matrices per record set:

  var_pa[f, r]        = 1 if founder f confirmed-carries the ALT at record r (GT 1/1 or any >0)
                        0 if confirmed REF (0/0) OR missing (./.)

  var_called[f, r] = 1 if founder f's GT is NOT ./. at record r (called genotype)
                        0 if GT is ./.

Projection convention (handles missing properly — see SESSION_2026-05-18.md var_pa bug):
    AF_est[r] = (h @ var_pa)[r] / (h @ var_called)[r]

For uniform h=1/F this reduces to AC/AN — the standard "AF among called samples".
The pre-fix convention (h @ var_pa, divide by 1) treats ./. as REF and systematically
under-counts AF at records with high F_MISSING (e.g. cactus-only or PG-only records
after the cactus_78 + PG_153 merge, which have F_MISSING up to 0.66).

Output: var_pa.npz, var_called.npz, meta.npz (with chrom/pos/ref/alt strings).
"""
from __future__ import annotations
import argparse, time
from pathlib import Path
import numpy as np
import pysam
from scipy.sparse import csr_matrix, save_npz


def build_var_pa(vcf_path: str, out_prefix: str, chrom_filter: str | None = None,
                 verbose: bool = True):
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    F = len(samples)
    print(f"Founders: F={F}")

    chrom_arr, pos_arr, ref_arr, alt_arr = [], [], [], []
    ref_len_arr, alt_len_arr = [], []
    carrier_rows, carrier_cols = [], []
    called_rows,  called_cols  = [], []
    n_kept = 0
    iterator = vcf.fetch(chrom_filter) if chrom_filter else vcf.fetch()

    t0 = time.time()
    for r_idx, rec in enumerate(iterator):
        if not rec.alts:
            continue
        chrom_arr.append(rec.chrom)
        pos_arr.append(rec.pos)
        ref_arr.append(rec.ref)
        alt_arr.append(rec.alts[0])
        ref_len_arr.append(len(rec.ref))
        alt_len_arr.append(len(rec.alts[0]))

        for f_idx, sname in enumerate(samples):
            gt = rec.samples[sname]["GT"]
            if gt is None or len(gt) == 0:
                # treat as ./.
                continue
            # Panel must be pre-haploidized; see build_kmer_pa.py for rationale.
            # Crash early on diploid GTs so var_pa doesn't silently disagree
            # with kmer_pa's GT[0]-only reconstruction.
            if len(gt) != 1:
                raise ValueError(
                    f"Panel VCF must be haploid (got len(GT)={len(gt)} for "
                    f"{sname} at {rec.chrom}:{rec.pos}). Haploidize the VCF "
                    f"before building var_pa."
                )
            # any None allele in GT marks the genotype as missing
            is_missing = all(a is None for a in gt)
            if is_missing:
                continue
            # genotype is called — mark in var_called
            called_rows.append(f_idx)
            called_cols.append(n_kept)
            # check ALT carriage
            has_alt = any(a is not None and a > 0 for a in gt)
            if has_alt:
                carrier_rows.append(f_idx)
                carrier_cols.append(n_kept)
        n_kept += 1
        if verbose and n_kept % 100000 == 0:
            print(f"  {n_kept:,} records, carriers={len(carrier_rows):,}, "
                  f"calls={len(called_rows):,}, elapsed={time.time()-t0:.0f}s")

    print(f"\nTotal records: {n_kept:,}")
    print(f"  var_pa carrier nnz: {len(carrier_rows):,}, "
          f"density: {len(carrier_rows)/(F*n_kept):.4f}")
    print(f"  var_called nnz: {len(called_rows):,}, "
          f"density: {len(called_rows)/(F*n_kept):.4f}")

    var_pa = csr_matrix(
        (np.ones(len(carrier_rows), dtype=np.int8), (carrier_rows, carrier_cols)),
        shape=(F, n_kept), dtype=np.int8,
    )
    var_called = csr_matrix(
        (np.ones(len(called_rows), dtype=np.int8), (called_rows, called_cols)),
        shape=(F, n_kept), dtype=np.int8,
    )

    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    save_npz(out_prefix + ".var_pa.npz", var_pa)
    save_npz(out_prefix + ".var_called.npz", var_called)
    np.savez(
        out_prefix + ".meta.npz",
        founders=np.array(samples),
        chrom=np.array(chrom_arr),
        pos=np.array(pos_arr, dtype=np.int64),
        # ref/alt use dtype=object so long SV alleles don't trigger fixed-width
        # numpy strings (some SV REFs are 100K+ bp and would otherwise allocate
        # 2.2 TB for the array).
        ref=np.array(ref_arr, dtype=object),
        alt=np.array(alt_arr, dtype=object),
        ref_len=np.array(ref_len_arr),
        alt_len=np.array(alt_len_arr),
    )
    print(f"Wrote {out_prefix}.var_pa.npz, {out_prefix}.var_called.npz, "
          f"{out_prefix}.meta.npz (with REF/ALT base strings)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="biallelic.norm.vcf.gz path")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--chrom", default=None, help="optional: limit to one chrom")
    args = ap.parse_args()
    build_var_pa(args.vcf, args.out, args.chrom)


if __name__ == "__main__":
    main()

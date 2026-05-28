"""
Build an ATOMIZED cn_var: per-base SNP catalog with founder carriers UNIONed across
all source biallelic VCF records that imply each per-base substitution.

For each source record (pos, ref, alt) with carrier set C and called set K:
  L = min(len(ref), len(alt))   # length of the aligned overlap region
  for i in 0..L-1:
      ref_b = ref[i]; alt_b = alt[i]; atom_pos = pos + i
      called_at_pos[(atom_pos, ref_b)] |= K          # this record can call REF base at atom_pos
      if ref_b != alt_b and ref_b in ACGT and alt_b in ACGT:
          carriers[(atom_pos, ref_b, alt_b)] |= C    # this record's ALT carriers carry the atomic substitution
          called_at_pos[(atom_pos, ref_b)] |= C      # defensive: carriers must be in called set

Output rows = unique (pos, ref_b, alt_b) tuples; per-row carrier and called sets are unioned.
Insertions/deletions beyond the aligned overlap region contribute NO atomized records (no
per-base SNP equivalent). Pure INS/DEL beyond REF/ALT[0] still don't atomize; they remain
in the path-aware raw cn_var for downstream SV-aware analyses.

N alleles (including IUPAC ambiguity codes that cactus normalizes to N) are skipped.
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
from collections import defaultdict
import numpy as np
import pysam
from scipy.sparse import csc_matrix, save_npz


def atomize_vcf(vcf_path: str, out_prefix: str, chrom_filter: str | None = None):
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    F = len(samples)
    print(f"Founders: F={F}", flush=True)

    # bool array (F,) per key — much cheaper than Python set[int] (~24 B/int).
    # F=231 bools = 231 bytes per array; for 7.5M keys → ~1.7 GB total carriers
    # plus ~1 GB called_at_pos. Hash-table overhead per dict entry is the dominant cost.
    def new_mask():
        return np.zeros(F, dtype=bool)
    carriers: dict[tuple, np.ndarray] = defaultdict(new_mask)
    called_at_pos: dict[tuple, np.ndarray] = defaultdict(new_mask)

    ACGT = {"A", "C", "G", "T"}
    iterator = vcf.fetch(chrom_filter) if chrom_filter else vcf.fetch()

    t0 = time.time()
    n_records = 0
    n_atomic_emit = 0
    n_skipped_n_allele = 0
    chrom_seen: set[str] = set()

    rec_carriers = np.zeros(F, dtype=bool)
    rec_called = np.zeros(F, dtype=bool)

    for rec in iterator:
        if not rec.alts:
            continue
        ref = rec.ref.upper()
        alt = rec.alts[0].upper()
        pos = rec.pos  # 1-based per pysam
        chrom_seen.add(rec.chrom)

        # Build per-record carrier + called bool masks (reuse scratch arrays)
        rec_carriers[:] = False
        rec_called[:] = False
        for f_idx, sname in enumerate(samples):
            gt = rec.samples[sname]["GT"]
            if gt is None or len(gt) == 0:
                continue
            if all(a is None for a in gt):
                continue
            rec_called[f_idx] = True
            if any(a is not None and a > 0 for a in gt):
                rec_carriers[f_idx] = True

        # Iterate over the aligned overlap region (positions covered by BOTH ref and alt)
        L = min(len(ref), len(alt))
        for i in range(L):
            ref_b = ref[i]
            alt_b = alt[i]
            if ref_b not in ACGT or alt_b not in ACGT:
                n_skipped_n_allele += 1
                continue
            atom_pos = pos + i

            # Record that this source allows us to call the REF base at atom_pos
            np.logical_or(called_at_pos[(atom_pos, ref_b)], rec_called,
                          out=called_at_pos[(atom_pos, ref_b)])

            if ref_b != alt_b:
                np.logical_or(carriers[(atom_pos, ref_b, alt_b)], rec_carriers,
                              out=carriers[(atom_pos, ref_b, alt_b)])
                # Carriers are by definition called at the REF position
                np.logical_or(called_at_pos[(atom_pos, ref_b)], rec_carriers,
                              out=called_at_pos[(atom_pos, ref_b)])
                n_atomic_emit += 1

        n_records += 1
        if n_records % 200_000 == 0:
            print(f"  {n_records:,} source records,  {len(carriers):,} unique atoms,  "
                  f"{n_atomic_emit:,} emissions,  elapsed={time.time()-t0:.0f}s", flush=True)

    print(f"\nProcessed {n_records:,} biallelic source records", flush=True)
    print(f"Unique atomized substitutions: {len(carriers):,}", flush=True)
    print(f"Skipped emissions due to N/ambiguous alleles: {n_skipped_n_allele:,}", flush=True)
    print(f"Source-record emissions to atomic sites: {n_atomic_emit:,}", flush=True)
    print(f"Chromosomes seen: {sorted(chrom_seen)}", flush=True)

    # Build CSC sparse matrices directly from per-column bool masks.
    # Avoids materializing 35M-element (row, col) lists.
    print("\nBuilding CSC sparse matrices column-by-column ...", flush=True)
    keys = sorted(carriers.keys())  # canonical order: by (pos, ref, alt)
    Ncols = len(keys)
    pos_arr = np.empty(Ncols, dtype=np.int64)
    chrom_arr = np.empty(Ncols, dtype=object)
    ref_arr = np.empty(Ncols, dtype=object)
    alt_arr = np.empty(Ncols, dtype=object)
    ref_len_arr = np.ones(Ncols, dtype=np.int64)
    alt_len_arr = np.ones(Ncols, dtype=np.int64)

    chrom_default = sorted(chrom_seen)[0] if chrom_seen else "Chr1"

    # CSC: indptr[N+1], indices[total_nnz], data[total_nnz]
    indptr_c = np.zeros(Ncols + 1, dtype=np.int64)
    indptr_called = np.zeros(Ncols + 1, dtype=np.int64)
    # We do not know nnz upfront; first pass to compute, second pass to fill.
    nnz_c = 0; nnz_called = 0
    for col_idx, key in enumerate(keys):
        atom_pos, ref_b, alt_b = key
        nnz_c += int(carriers[key].sum())
        nnz_called += int(called_at_pos[(atom_pos, ref_b)].sum())
    print(f"  cn_var nnz total: {nnz_c:,}", flush=True)
    print(f"  cn_var_called nnz total: {nnz_called:,}", flush=True)

    indices_c = np.empty(nnz_c, dtype=np.int32)
    indices_called = np.empty(nnz_called, dtype=np.int32)
    cur_c = 0
    cur_called = 0
    for col_idx, key in enumerate(keys):
        atom_pos, ref_b, alt_b = key
        pos_arr[col_idx] = atom_pos
        chrom_arr[col_idx] = chrom_default
        ref_arr[col_idx] = ref_b
        alt_arr[col_idx] = alt_b

        cmask = carriers[key]
        kmask = called_at_pos[(atom_pos, ref_b)]
        c_idx = np.nonzero(cmask)[0]
        k_idx = np.nonzero(kmask)[0]
        n_c = c_idx.shape[0]; n_k = k_idx.shape[0]
        indices_c[cur_c:cur_c + n_c] = c_idx
        indices_called[cur_called:cur_called + n_k] = k_idx
        cur_c += n_c; cur_called += n_k
        indptr_c[col_idx + 1] = cur_c
        indptr_called[col_idx + 1] = cur_called
        if (col_idx + 1) % 500_000 == 0:
            print(f"  {col_idx+1:,}/{Ncols:,} columns filled, elapsed={time.time()-t0:.0f}s", flush=True)

    data_c = np.ones(nnz_c, dtype=np.int8)
    data_called = np.ones(nnz_called, dtype=np.int8)

    cn_var = csc_matrix((data_c, indices_c, indptr_c), shape=(F, Ncols), dtype=np.int8)
    cn_var_called = csc_matrix((data_called, indices_called, indptr_called), shape=(F, Ncols), dtype=np.int8)

    print(f"cn_var: {cn_var.shape}  {cn_var.nnz:,} nnz  density={cn_var.nnz/(F*Ncols):.4f}", flush=True)
    print(f"cn_var_called: {cn_var_called.shape}  {cn_var_called.nnz:,} nnz  density={cn_var_called.nnz/(F*Ncols):.4f}", flush=True)

    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    save_npz(out_prefix + ".cn_var.npz", cn_var)
    save_npz(out_prefix + ".cn_var_called.npz", cn_var_called)
    np.savez(
        out_prefix + ".meta.npz",
        founders=np.array(samples),
        chrom=chrom_arr,
        pos=pos_arr,
        ref=ref_arr,
        alt=alt_arr,
        ref_len=ref_len_arr,
        alt_len=alt_len_arr,
    )
    print(f"\nWrote: {out_prefix}.{{cn_var,cn_var_called,meta}}.npz", flush=True)
    print(f"Total time: {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="biallelic VCF (e.g., merged_231_chr1_final.vcf.gz)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--chrom", default=None, help="optional: limit to one chromosome")
    args = ap.parse_args()
    atomize_vcf(args.vcf, args.out, args.chrom)

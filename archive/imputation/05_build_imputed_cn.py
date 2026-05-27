"""
Convert imputed VCF (231 GrENE-Net founders × cactus SVs) into a
231 × K_sv cn matrix that plugs into the pool-seq pipeline.

Steps:
1. Read imputed_151.vcf.gz (151 GrENE-only founders with imputed SV genotypes)
2. Read cactus_svs_renamed.vcf.gz (80 cactus founders with truth SV genotypes)
3. Merge into one 231-founder × K_sv matrix
4. Save as scipy sparse CSR (`cn_imputed_231.npz`) + metadata
"""
from __future__ import annotations
import sys, os
import numpy as np
import pysam
from scipy.sparse import csr_matrix, save_npz


WORK = "/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work"
OUT = "/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/cn_imputed_231"


def main():
    cactus = pysam.VariantFile(os.path.join(WORK, "cactus_svs_renamed.vcf.gz"))
    imputed = pysam.VariantFile(os.path.join(WORK, "imputed_151.vcf.gz"))

    cactus_samples = list(cactus.header.samples)  # 80
    imputed_samples = list(imputed.header.samples)  # 151
    all_samples = list(cactus_samples) + list(imputed_samples)  # 231
    F = len(all_samples)
    print(f"Total founders: {F}  (cactus={len(cactus_samples)}, imputed={len(imputed_samples)})")

    # Index variant records by (chrom, pos, REF, ALT) to match across VCFs
    # Note: Beagle output may differ in REF/ALT representation for normalized SVs.
    # We match by (chrom, pos) since both are biallelic-normalized at the same sites.

    # Build records list from cactus (the ground truth subset of variants)
    records = []
    for rec in cactus.fetch():
        # SV filter (matching the prep step)
        sv_size = abs(len(rec.ref) - len(rec.alts[0]))
        if sv_size < 50 or sv_size > 50000:
            continue
        records.append({
            "chrom": rec.chrom,
            "pos": rec.pos,
            "ref": rec.ref,
            "alt": rec.alts[0],
            "size": sv_size,
        })
    K = len(records)
    print(f"SV records: {K:,}")

    # Build cn matrix: F × K, dtype int8
    # For cactus samples: known dosage (0/1, since haploid → diploid is 0|0 or 1|1)
    # For imputed samples: dosage from Beagle GT or DS
    cn = np.zeros((F, K), dtype=np.int8)

    # Index records for fast lookup
    rec_index = {(r["chrom"], r["pos"]): i for i, r in enumerate(records)}
    sample_index = {s: i for i, s in enumerate(all_samples)}

    # Fill cactus rows (truth)
    cactus = pysam.VariantFile(os.path.join(WORK, "cactus_svs_renamed.vcf.gz"))
    for rec in cactus.fetch():
        sv_size = abs(len(rec.ref) - len(rec.alts[0]))
        if sv_size < 50 or sv_size > 50000:
            continue
        k_idx = rec_index.get((rec.chrom, rec.pos))
        if k_idx is None: continue
        for s in cactus_samples:
            gt = rec.samples[s]["GT"]
            if gt is None or len(gt) == 0: continue
            # Inbred → 0|0 or 1|1; we record 1 if any alt allele present
            dose = sum(1 for g in gt if g is not None and g > 0) // max(1, len(gt))
            cn[sample_index[s], k_idx] = 1 if dose > 0 else 0

    # Fill imputed rows
    imputed = pysam.VariantFile(os.path.join(WORK, "imputed_151.vcf.gz"))
    for rec in imputed.fetch():
        # Same SV filter
        try:
            sv_size = abs(len(rec.ref) - len(rec.alts[0]))
        except Exception:
            continue
        if sv_size < 50 or sv_size > 50000:
            continue
        k_idx = rec_index.get((rec.chrom, rec.pos))
        if k_idx is None: continue
        for s in imputed_samples:
            gt = rec.samples[s]["GT"]
            if gt is None or len(gt) == 0: continue
            # For imputed samples, use the rounded dosage (0/1 since inbred)
            dose = sum(1 for g in gt if g is not None and g > 0)
            # If both alleles agree (0|0 or 1|1) → 0 or 1
            # If discordant (0|1 — Beagle posterior may produce this from imputation
            # uncertainty), pick majority
            cn[sample_index[s], k_idx] = 1 if dose > 0 else 0

    # Save
    print(f"cn shape: {cn.shape}, nnz: {(cn>0).sum():,}")
    cn_sparse = csr_matrix(cn)
    save_npz(OUT + ".cn.npz", cn_sparse)
    np.savez(OUT + ".meta.npz",
             founders=np.array(all_samples),
             chrom=np.array([r["chrom"] for r in records]),
             pos=np.array([r["pos"] for r in records], dtype=np.int64),
             ref_len=np.array([len(r["ref"]) for r in records]),
             alt_len=np.array([len(r["alt"]) for r in records]),
             size=np.array([r["size"] for r in records]))
    print(f"Wrote {OUT}.cn.npz and {OUT}.meta.npz")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Build founder × SV genotype matrix from the SV panel VCF.

Output:
  data/founder_sv_matrix.parquet
    - rows = SV records keyed by (chrom, pos, alt_idx, ref_len, alt_len, var_type, sv_size)
    - cols = ecotype IDs (mapped from Assembly_ID via ASSEMBLIES_Best_version_of_dataset.csv)
    - values = 0/1 dosage of ALT allele (1 = haploid carrier, 0 = REF)

Restrict to SVs (var_type DEL/INS, |ref_len - alt_len| >= 50 bp). For INS, multiple
SV-sized ALTs at the same position are kept as separate rows (alt_idx 1, 2, ...).

Genotypes in the panel are biallelic 0/0 or 1/1 per accession (haploid model). When
bcftools merge created a multi-allelic site, an accession's GT for ALT_i is 1 iff that
accession's per-sample VCF carried that exact ALT_i.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

SV_MIN = 50  # ≥ 50 bp


def load_assembly_to_accession() -> dict:
    """Assembly_ID -> Accession_ID (string, both numeric)."""
    df = pd.read_csv(
        "/home/tbellagio/scratch/pang/long_read_seq_ara/ASSEMBLIES_Best_version_of_dataset.csv",
        dtype=str,
    )
    df = df.dropna(subset=["Assembly_ID", "Accession_ID"])
    df["Accession_ID"] = df["Accession_ID"].astype(str).str.split(".").str[0]
    return dict(zip(df["Assembly_ID"], df["Accession_ID"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--vcf",
        default="/home/tbellagio/scratch/pang/sv_panel/merge_vcfs_bcftools/panel.snp_ins_del.merged.no_sv_singletons.vcf.gz",
    )
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--out-meta", required=True, type=Path)
    args = ap.parse_args()

    asm_to_acc = load_assembly_to_accession()
    print(f"Loaded {len(asm_to_acc):,} Assembly_ID → Accession_ID mappings")

    vcf = pysam.VariantFile(str(args.vcf))
    panel_assemblies = list(vcf.header.samples)
    print(f"SV panel has {len(panel_assemblies)} accessions (Assembly_IDs)")

    panel_to_acc = {asm: asm_to_acc[asm] for asm in panel_assemblies if asm in asm_to_acc}
    panel_acc_ids = list(panel_to_acc.values())
    seen = set()
    dup = [a for a in panel_acc_ids if a in seen or seen.add(a)]
    if dup:
        print(f"WARNING: duplicate accession IDs in panel after mapping: {dup}")

    panel_assemblies_kept = [a for a in panel_assemblies if a in panel_to_acc]
    panel_accessions      = [panel_to_acc[a] for a in panel_assemblies_kept]
    print(f"Mapped {len(panel_accessions)} of {len(panel_assemblies)} panel accessions")

    # Asyc-check mapping consistency
    print("First 5 mappings:")
    for a in panel_assemblies_kept[:5]:
        print(f"  {a} -> {panel_to_acc[a]}")

    # Iterate VCF and select SV-sized records (any size diff ≥ 50 between REF and ALT_i)
    rows = []
    cols = panel_accessions  # final ecotype IDs as columns
    asm_idx_lookup = {s: i for i, s in enumerate(panel_assemblies)}
    asm_idx = [asm_idx_lookup[a] for a in panel_assemblies_kept]
    print(f"Iterating VCF: {args.vcf}")
    t0 = time.time()
    n_records = 0
    n_sv_alts = 0
    geno_chunks = []
    meta_chunks = []
    for rec in vcf.fetch():
        n_records += 1
        if n_records % 250_000 == 0:
            print(f"  {n_records:,} records, {n_sv_alts:,} SV ALTs, {time.time()-t0:.0f}s")
        ref = rec.ref
        alts = rec.alts or ()
        ref_len = len(ref)
        for alt_idx, alt in enumerate(alts, start=1):
            if alt is None or alt == "*":
                continue
            alt_len = len(alt)
            size_diff = abs(alt_len - ref_len)
            if size_diff < SV_MIN:
                continue
            var_type = "INS" if alt_len > ref_len else ("DEL" if alt_len < ref_len else "OTHER")
            if var_type == "OTHER":
                continue
            # extract genotypes
            geno = np.zeros(len(asm_idx), dtype=np.int8)
            for j, sidx in enumerate(asm_idx):
                gt = rec.samples[sidx].get("GT", (None,))
                if gt is None or len(gt) == 0:
                    continue
                # haploid model — many accessions report (0,) or (1,)
                # bcftools merge can give (0, 0) or (1, 1) — treat any 1 as ALT carrier of that ALT_i
                # alt_idx 1 corresponds to GT value 1, alt_idx 2 to GT value 2, etc.
                # so check if any GT element equals alt_idx
                if any(a == alt_idx for a in gt if a is not None):
                    geno[j] = 1
            ac = int(geno.sum())
            if ac == 0:
                # nobody in mapped subset carries this ALT — drop (no information)
                continue
            geno_chunks.append(geno)
            meta_chunks.append({
                "chrom": rec.contig,
                "pos": rec.pos,
                "alt_idx": alt_idx,
                "ref_len": ref_len,
                "alt_len": alt_len,
                "var_type": var_type,
                "sv_size": size_diff,
                "ac_in_panel": ac,
            })
            n_sv_alts += 1

    print(f"Total: {n_records:,} records → {n_sv_alts:,} SV ALTs (≥{SV_MIN}bp) with ≥1 carrier in mapped subset")
    if not geno_chunks:
        print("No SV ALTs survived. Exiting.")
        sys.exit(1)

    G = np.vstack(geno_chunks).astype(np.int8)
    meta = pd.DataFrame(meta_chunks)
    print(f"G shape: {G.shape}  (SV ALTs × ecotypes)")

    # write founder × SV matrix as parquet (long form: row_idx, ecotype, dosage)
    # but parquet handles wide better — store as dense matrix with ecotype columns
    G_df = pd.DataFrame(G, columns=cols)
    G_df.insert(0, "chrom", meta["chrom"].astype("category"))
    G_df.insert(1, "pos", meta["pos"].astype(np.int32))
    G_df.insert(2, "alt_idx", meta["alt_idx"].astype(np.int8))
    G_df.to_parquet(args.out, compression="zstd")
    meta.to_parquet(args.out_meta, compression="zstd")
    print(f"wrote: {args.out}  ({args.out.stat().st_size/1e6:.1f} MB)")
    print(f"wrote: {args.out_meta}")

    # quick stats
    print("\nVariant type distribution:")
    print(meta["var_type"].value_counts())
    print("\nSV size quantiles:")
    print(meta["sv_size"].describe())
    print("\nAC in mapped panel quantiles:")
    print(meta["ac_in_panel"].describe())


if __name__ == "__main__":
    main()

#!/usr/bin/env python -u
"""
Tier 2 (gtcheck substitute): pairwise SNP-genotype identity, streamed in chunks
so memory stays bounded. Works on haploid input directly (no diploid conversion).

For each pair (extras_i, cactus_j), compute identity = #records where both called
AND have the same allele / #records both called. Outputs 53×82 matrix.
"""
import gzip
import subprocess
import sys
import numpy as np
import pandas as pd
from pathlib import Path

WORK = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_135_vs_82")
SNPS = WORK / "data" / "pang135_biallelic_snps.vcf.gz"   # haploid biallelic SNPs (already built)
SAMPLES_135 = WORK / "data" / "pang_135_samples.txt"
GROUP82 = WORK / "data" / "grenenet_82_in_pang135.txt"
GROUP53 = WORK / "data" / "extras_53.txt"
OUT_MAT = WORK / "results" / "tier2_snp_identity_matrix.tsv"
OUT_DUPS = WORK / "results" / "tier2_quasi_duplicates.tsv"
THRESH = 0.99

def main():
    samples = [s.strip() for s in open(SAMPLES_135) if s.strip()]
    group82 = set(s.strip() for s in open(GROUP82) if s.strip())
    group53 = set(s.strip() for s in open(GROUP53) if s.strip())
    n = len(samples)
    idx82 = [i for i, s in enumerate(samples) if s in group82]
    idx53 = [i for i, s in enumerate(samples) if s in group53]
    assert len(idx82) == 82
    assert len(idx53) == 53
    print(f"[init] 82 cactus indices, 53 extras among {n} samples")

    # Streaming bcftools query
    cmd = ["bcftools", "query", "-f", "[%GT\t]\n", str(SNPS)]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=8192)

    # Running totals for pairs (53 x 82)
    matches = np.zeros((53, 82), dtype=np.int64)
    both_called = np.zeros((53, 82), dtype=np.int64)

    CHUNK = 100_000
    chunk = []
    nrec = 0

    def process_chunk(chunk):
        nonlocal matches, both_called
        m = len(chunk)
        # Build (m, n) int8 matrix: 0 → REF, 1 → ALT, -1 → missing
        gt = np.full((m, n), -1, dtype=np.int8)
        for i, tokens in enumerate(chunk):
            for j, tok in enumerate(tokens[:n]):
                if tok == "0":   gt[i, j] = 0
                elif tok == "1": gt[i, j] = 1
                # else stays -1
        A = gt[:, idx53]   # (m, 53)
        B = gt[:, idx82]   # (m, 82)
        # masks of called
        mA = (A != -1).astype(np.int32)
        mB = (B != -1).astype(np.int32)
        both_called += mA.T @ mB
        # matches: (A==0 ∧ B==0) + (A==1 ∧ B==1)
        A0 = (A == 0).astype(np.int32); A1 = (A == 1).astype(np.int32)
        B0 = (B == 0).astype(np.int32); B1 = (B == 1).astype(np.int32)
        matches += A0.T @ B0 + A1.T @ B1

    for line in p.stdout:
        tokens = line.rstrip("\n").split("\t")
        chunk.append(tokens)
        nrec += 1
        if len(chunk) >= CHUNK:
            process_chunk(chunk)
            chunk = []
            print(f"[stream] {nrec:,} biallelic SNPs processed")
    if chunk:
        process_chunk(chunk)
    p.wait()
    print(f"[done] {nrec:,} biallelic SNPs streamed")

    with np.errstate(divide='ignore', invalid='ignore'):
        identity = np.where(both_called > 0, matches / both_called, np.nan)

    extras = [samples[i] for i in idx53]
    cactus = [samples[i] for i in idx82]

    df = pd.DataFrame(identity, index=extras, columns=cactus)
    df.index.name = "extras_id"
    df.to_csv(OUT_MAT, sep="\t")
    print(f"[out] identity matrix {df.shape} → {OUT_MAT}")

    pairs = []
    for i, e in enumerate(extras):
        for j, c in enumerate(cactus):
            v = identity[i, j]
            if not np.isnan(v) and v >= THRESH:
                pairs.append((e, c, float(v), int(both_called[i, j])))
    dups = pd.DataFrame(pairs, columns=["extras_id", "cactus_id", "snp_identity", "n_called_snps"])
    dups = dups.sort_values("snp_identity", ascending=False)
    dups.to_csv(OUT_DUPS, sep="\t", index=False)
    print(f"[out] {len(dups)} pairs with identity >= {THRESH} → {OUT_DUPS}")

    print("\n[summary]")
    print(f"  records: {nrec:,} biallelic SNPs")
    print(f"  median identity (53×82): {np.nanmedian(identity):.4f}")
    print(f"  max identity per extras (top 10):")
    max_per_extras = pd.Series(np.nanmax(identity, axis=1), index=extras).sort_values(ascending=False)
    print(max_per_extras.head(10).to_string())
    print(f"\n  max identity per cactus (top 5):")
    max_per_cactus = pd.Series(np.nanmax(identity, axis=0), index=cactus).sort_values(ascending=False)
    print(max_per_cactus.head(5).to_string())

if __name__ == "__main__":
    main()

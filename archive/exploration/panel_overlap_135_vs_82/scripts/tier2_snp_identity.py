#!/usr/bin/env python
"""
Tier 2 (substitute for k-mer Jaccard):
Pairwise SNP-genotype identity between the 53 extras and 82 GrENE-Net cactus,
computed on biallelic SNPs from pang_135 raw VCF.

This is the same metric that caught the 5772/6150 cactus mislabel (99.81% SNP
identity). Quasi-duplicates of any 82-cactus assembly contribute zero new
genotyping information to the 153 PG founders, so should be flagged.

Output:
  panel_overlap_135_vs_82/results/tier2_snp_identity_matrix.tsv   (53 x 82)
  panel_overlap_135_vs_82/results/tier2_quasi_duplicates.tsv      (pairs >= threshold)
"""
import gzip
import subprocess
import sys
import numpy as np
import pandas as pd
from pathlib import Path

WORK = Path(__file__).resolve().parents[1]
RAW = "/global/home/users/tbellg/scratch/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.raw.vcf.gz"
SAMPLES_135 = WORK / "data" / "pang_135_samples.txt"
GROUP82 = WORK / "data" / "grenenet_82_in_pang135.txt"
GROUP53 = WORK / "data" / "extras_53.txt"
OUT_MAT = WORK / "results" / "tier2_snp_identity_matrix.tsv"
OUT_DUPS = WORK / "results" / "tier2_quasi_duplicates.tsv"
THRESH = 0.99  # flag identity >= this

def main():
    samples = [s.strip() for s in open(SAMPLES_135) if s.strip()]
    group82 = set(s.strip() for s in open(GROUP82) if s.strip())
    group53 = set(s.strip() for s in open(GROUP53) if s.strip())
    n = len(samples)
    idx82 = [i for i, s in enumerate(samples) if s in group82]
    idx53 = [i for i, s in enumerate(samples) if s in group53]
    assert len(idx82) == 82, f"expected 82, got {len(idx82)}"
    assert len(idx53) == 53, f"expected 53, got {len(idx53)}"
    print(f"[init] 82 cactus indices, 53 extras indices among {n} samples")

    # Pull biallelic SNPs only — fast. Haploid cactus output: each GT cell is "0" or "1" (or ".").
    # bcftools view -v snps -m2 -M2 filters to biallelic SNPs.
    # Stream GT field for all samples.
    cmd = [
        "bcftools", "view", "--threads", "4", "-v", "snps", "-m2", "-M2",
        "-r", "Chr1,Chr2,Chr3,Chr4,Chr5", RAW
    ]
    cmd2 = [
        "bcftools", "query", "-f", "[%GT\t]\n"
    ]
    # build via Popen pipe
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
    p1.stdout.close()

    # Encode GTs as int8:  '0' → 0, '1' → 1, anything else → -1 (missing/non-call)
    gt_mat_chunks = []  # list of int8 arrays per chunk
    nrec = 0
    chunk = []
    CHUNK = 200_000
    for line in p2.stdout:
        tokens = line.rstrip("\n").split("\t")[:n]
        row = np.empty(n, dtype=np.int8)
        for i, tok in enumerate(tokens):
            # haploid: tok is "0" or "1" or "."
            if tok == "0":
                row[i] = 0
            elif tok == "1":
                row[i] = 1
            else:
                row[i] = -1
        chunk.append(row)
        nrec += 1
        if len(chunk) == CHUNK:
            gt_mat_chunks.append(np.stack(chunk))
            chunk = []
            print(f"[stream] {nrec} biallelic SNPs read")
    if chunk:
        gt_mat_chunks.append(np.stack(chunk))
    p2.wait()
    p1.wait()

    if not gt_mat_chunks:
        print("ERROR: no records parsed", file=sys.stderr)
        sys.exit(2)
    gt = np.concatenate(gt_mat_chunks, axis=0)  # (nrec, n_samples)
    print(f"[matrix] gt shape {gt.shape}")

    # Compute pairwise identity for (53 extras) x (82 cactus):
    # identity(i, j) = #records where both non-missing AND gt[*,i] == gt[*,j]  /  #records both non-missing
    A = gt[:, idx53]  # (nrec, 53)
    B = gt[:, idx82]  # (nrec, 82)

    # Use int8 carefully. Build masks of non-missing.
    mask_A = (A != -1).astype(np.int32)  # (nrec, 53)
    mask_B = (B != -1).astype(np.int32)  # (nrec, 82)

    # Pairwise mask sum: how many records both non-missing
    both_called = mask_A.T @ mask_B  # (53, 82)

    # Pairwise matches: for non-missing positions, count where genotypes equal.
    # Trick: at non-missing positions both A and B in {0,1}. Match count =
    #   (A==1 ∧ B==1) + (A==0 ∧ B==0). With masks, equivalent to:
    #   (A_masked.T @ B_masked) + ((1-A)_masked.T @ (1-B)_masked)
    # where _masked means we zero out the missing positions.
    A0 = np.where(A == 0, 1, 0).astype(np.int32)  # 1 where call is REF
    A1 = np.where(A == 1, 1, 0).astype(np.int32)  # 1 where call is ALT
    B0 = np.where(B == 0, 1, 0).astype(np.int32)
    B1 = np.where(B == 1, 1, 0).astype(np.int32)
    matches = A0.T @ B0 + A1.T @ B1  # (53, 82)

    with np.errstate(divide='ignore', invalid='ignore'):
        identity = np.where(both_called > 0, matches / both_called, np.nan)

    extras = [samples[i] for i in idx53]
    cactus = [samples[i] for i in idx82]

    df = pd.DataFrame(identity, index=extras, columns=cactus)
    df.index.name = "extras_id"
    df.to_csv(OUT_MAT, sep="\t")
    print(f"[out] identity matrix {df.shape} → {OUT_MAT}")

    # Long-form quasi-duplicates table
    pairs = []
    for i, e in enumerate(extras):
        for j, c in enumerate(cactus):
            v = identity[i, j]
            if np.isnan(v):
                continue
            if v >= THRESH:
                pairs.append((e, c, float(v), int(both_called[i, j])))
    dups = pd.DataFrame(pairs, columns=["extras_id", "cactus_id", "snp_identity", "n_called_snps"])
    dups = dups.sort_values("snp_identity", ascending=False)
    dups.to_csv(OUT_DUPS, sep="\t", index=False)
    print(f"[out] {len(dups)} pairs with identity >= {THRESH} → {OUT_DUPS}")

    # Headline summary
    print(f"\n[summary]")
    print(f"  records used: {gt.shape[0]} biallelic SNPs")
    print(f"  median identity (53×82): {np.nanmedian(identity):.4f}")
    print(f"  max identity per extras assembly: top 10")
    max_per_extras = pd.Series(np.nanmax(identity, axis=1), index=extras).sort_values(ascending=False)
    print(max_per_extras.head(10).to_string())

if __name__ == "__main__":
    main()

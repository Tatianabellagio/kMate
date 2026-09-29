#!/usr/bin/env python
"""Genomic length share of each context class -- the denominator for "where are they found" (1c of
sv_selection_key_results). Exact base-pair coverage of Chr1-5 under the SAME first-match rule as
`_sv_sfs_context.context()`: CDS > UTR > intron (gene & ~exon) > TE > promoter (2 kb upstream,
strand-aware) > intergenic. Chromosome lengths from the TAIR10 reference used for degeneracy.
Env: kmate. Writes results/sv_adaptive/context_genome_share.csv (context, bp, share).
"""
import os, sys, gzip, importlib.util
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
_H = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("sfsctx", os.path.join(_H, "_sv_sfs_context.py"))
sc = importlib.util.module_from_spec(_s); _s.loader.exec_module(sc)
OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
FA = f"{OUT.rsplit('/', 1)[0]}/ref/TAIR10.fa.gz"


def chrom_lengths():
    L, name, n = {}, None, 0
    with gzip.open(FA, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                if name: L[name] = n
                name = line[1:].split()[0]; n = 0
            else:
                n += len(line.strip())
    L[name] = n
    return {f"Chr{i}": L[str(i)] for i in range(1, 6)}


def mark(arr, iv, c, val):
    if c in iv:
        for s, e in zip(*iv[c]):
            arr[max(s, 1):min(e, arr.size - 1) + 1] = val


def main():
    IV = sc.layers(); L = chrom_lengths(); bp = np.zeros(len(sc.CTX), np.int64)
    for c in lib.CHROMS:
        a = np.full(L[c] + 1, 5, np.int8)                       # 1-based; index 0 unused
        mark(a, IV["promoter"], c, 4); mark(a, IV["TE"], c, 3)
        g = np.zeros(L[c] + 1, bool); mark(g, IV["gene"], c, True)
        x = np.zeros(L[c] + 1, bool); mark(x, IV["exon"], c, True)
        a[g & ~x] = 2
        mark(a, IV["UTR"], c, 1); mark(a, IV["CDS"], c, 0)
        bp += np.bincount(a[1:], minlength=6)
        print(f"[{c}] {L[c]:,} bp", flush=True)
    T = pd.DataFrame({"context": sc.CTX, "bp": bp, "share": bp / bp.sum()})
    T.to_csv(f"{OUT}/context_genome_share.csv", index=False)
    print(T.to_string(index=False)); print(f"[wrote] {OUT}/context_genome_share.csv")


if __name__ == "__main__":
    main()

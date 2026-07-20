#!/usr/bin/env python
"""How multi-haplotype are the recomputed (CLQcut=0.9) blocks, on the FOUNDERS?

For each block, restrict the 231-founder all-class genotype matrix (same MAF/called
filter + unique positions used to DEFINE the blocks) to that block's variants, and
count haplotypes:
  - n_distinct : number of unique founder haplotype rows (missing -> per-variant major)
  - n_eff      : effective # haplotypes = 1 / sum(freq^2)  (inverse Simpson)
  - n_ge2      : # haplotypes carried by >=2 founders (drops private singletons)
Summarize the distribution, and how it relates to HapFM's hardcoded k=7 cluster
trigger and to "essentially biallelic" (n_eff<=2).

Run in the kmate env. Reuses recompute_blocks.build_common_matrix for an exact match
to how the blocks were defined.
"""
import sys, os
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from recompute_blocks import build_common_matrix

BR = "analysis/grenenet_gea/blocks_recompute"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
MAF, MINCF = 0.05, 0.5


def hap_counts(sub):
    """sub: F x m int 0/1 founder matrix -> (n_distinct, n_eff, n_ge2)."""
    F = sub.shape[0]
    # unique rows + their multiplicity among the F founders
    _, counts = np.unique(sub, axis=0, return_counts=True)
    p = counts / F
    n_distinct = len(counts)
    n_eff = 1.0 / np.sum(p ** 2)
    n_ge2 = int(np.sum(counts >= 2))
    return n_distinct, n_eff, n_ge2


def main():
    allrows = []
    for chrom in CHROMS:
        chrlc = chrom.lower()
        pref = f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}"
        bf = f"{BR}/{chrlc}_clq0.9_blocks_clq0.9.tsv"
        print(f"[{chrom}] building founder matrix ...", flush=True)
        _, raw, positions = build_common_matrix(pref, MAF, MINCF)
        positions = np.asarray(positions)
        geno = (raw >= 0.5).astype(np.int8)   # missing was mean-imputed -> rounds to major
        bl = pd.read_csv(bf, sep="\t")
        rows = []
        for s, e in zip(bl["start_pos"], bl["end_pos"]):
            lo = int(np.searchsorted(positions, s))
            hi = int(np.searchsorted(positions, e, side="right"))
            if hi - lo < 2:
                continue
            nd, ne, ng = hap_counts(geno[:, lo:hi])
            rows.append((chrom, s, e, hi - lo, nd, ne, ng))
        df = pd.DataFrame(rows, columns=["chrom", "start", "end", "nvar",
                                         "n_distinct", "n_eff", "n_ge2"])
        allrows.append(df)
        print(f"[{chrom}] {len(df)} testable blocks (>=2 variants)")
    D = pd.concat(allrows, ignore_index=True)
    out = f"{BR}/block_haplotype_counts_clq0.9.csv"
    D.to_csv(out, index=False)

    n = len(D)
    print(f"\n=== ALL CHROMS: {n} testable blocks (>=2 variants) ===")
    print(f"variants/block         median {D.nvar.median():.0f}  "
          f"q25 {D.nvar.quantile(.25):.0f}  q75 {D.nvar.quantile(.75):.0f}  max {D.nvar.max()}")
    for col, lab in [("n_distinct", "distinct haps "),
                     ("n_eff", "EFFECTIVE haps"),
                     ("n_ge2", "haps (>=2 fdr)")]:
        v = D[col]
        print(f"{lab}  median {v.median():.1f}  q25 {v.quantile(.25):.1f}  "
              f"q75 {v.quantile(.75):.1f}  q95 {v.quantile(.95):.1f}  max {v.max():.0f}")
    print("\n--- 'how multiallelic?' using EFFECTIVE # haplotypes (inverse Simpson) ---")
    for lo, hi, lab in [(0, 2, "~biallelic   n_eff<=2 "),
                        (2, 4, "mild         2<n_eff<=4"),
                        (4, 7, "moderate     4<n_eff<=7"),
                        (7, 12, "high         7<n_eff<=12"),
                        (12, 1e9, "very high    n_eff>12 ")]:
        m = (D.n_eff > lo) & (D.n_eff <= hi)
        print(f"  {lab}: {m.sum():6d}  ({100*m.mean():.1f}%)")
    print("\n--- HapFM clustering trigger (uses DISTINCT haps > 7) ---")
    trig = D.n_distinct > 7
    print(f"  blocks with >7 distinct founder haplotypes (HapFM would cluster): "
          f"{trig.sum()} ({100*trig.mean():.1f}%)")
    print(f"  blocks with <=7 distinct (HapFM keeps each haplotype as-is):     "
          f"{(~trig).sum()} ({100*(~trig).mean():.1f}%)")
    print(f"\nwrote per-block table -> {out}")


if __name__ == "__main__":
    main()

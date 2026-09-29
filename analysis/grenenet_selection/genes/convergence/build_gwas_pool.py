#!/usr/bin/env python
"""GWAS candidate pool: every Bonferroni-significant NON-SNP marker, 30 per-garden scans.

What goes in
------------
`results/gwas_hits_allele_resolved_records.csv` (resolve_gwas_alleles.py) -- the same 30
per-garden GEMMA scans on the founder selection coefficient `s`, non-SNP classes only
(`nonsnp`, `sv`; `nonsnp` is smallindel+sv pooled, so a marker can appear in both), but
with the ALLELE of each hit recovered. The upstream `persite_bonferroni_hits.csv` is keyed
on position alone; taking the largest panel record at the position (what this script did
until 2026-09-17) named the wrong allele at 21 of 33 multiallelic markers, which is why
CYP71B4 / CML25 / EPFL5 were drawn with alleles significant in no garden.

**MAC floor stays at 5 (MAF 2.16%), per the user 2026-08-26.** The MAF-5% arm in
`gemma_gwas_mac12/` is deliberately NOT used. Raising the floor would not have bought
calibration anyway: the retained strata are still inflated 10.2x / 5.6x / 18.9x at
p<1e-6 (`gemma_gwas_mac12/mac_calibration.csv`), it only removes the worst stratum.
The MAC is carried per row instead so low-MAC candidates stay visible and triage-able.

Recurrence
----------
The GWAS analogue of the GEA's cross-axis recurrence is **how many gardens** a marker
is significant in. Unlike climate axes, gardens are genuinely independent replicates
(separate plantings), so `n_gardens` needs no cluster correction -- but the gardens are
not evenly powered: 6 of 30 gardens carry almost all hits, and `bio1_min/max` records
the climate span over which a marker recurs.

Variant size
------------
ref_len/alt_len now come with the hit, so `size` and `vclass` are the real allele's and
`size_inferred` is always False. `n_alleles_at_pos` / `pos_multiallelic` still flag
positions carrying more than one tested allele, and `key_rank` distinguishes same-length
ALT sequences (the second layer of `panel-multiallelic-pos-key-trap`): the pool is keyed
on (chrom, pos, ref_len, alt_len, key_rank), so two alleles at one position stay two rows.

Outputs -> results/
  gwas_pool.csv            one row per unique marker
  gwas_pool_records.csv    one row per (marker x class x garden) significant hit

env: kmate.  Run on a compute node.
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                    # noqa: E402

OUT = f"{HERE}/results"
HITS = f"{OUT}/gwas_hits_allele_resolved_records.csv"

CLASSES = ["nonsnp", "sv"]        # non-SNP layer only; the snp scan is not the pool


def main():
    os.makedirs(OUT, exist_ok=True)
    H = pd.read_csv(HITS)
    print(f"all Bonferroni hits: {len(H):,}")
    R = H[H.cls.isin(CLASSES)].copy()
    print(f"non-SNP classes only: {len(R):,} "
          f"({R.cls.value_counts().to_dict()})")

    key = ["chrom", "pos", "ref_len", "alt_len", "key_rank"]
    R = R.sort_values("nlp", ascending=False)
    g = R.groupby(key, sort=False)
    C = g.agg(
        best_nlp=("nlp", "max"),
        best_garden=("garden", "first"),
        mac=("mac", "min"),
        maf=("maf", "min"),
        bio1_min=("bio1", "min"),
        bio1_max=("bio1", "max"),
        n_records=("garden", "size"),
    ).reset_index()
    C["n_gardens"] = g["garden"].nunique().values
    C["gardens"] = g["garden"].apply(
        lambda s: ",".join(str(x) for x in dict.fromkeys(s))).values
    C["sig_classes"] = g["cls"].apply(lambda s: ",".join(sorted(set(s)))).values

    C["size"] = (C.alt_len - C.ref_len).abs()
    C["vclass"] = np.where(C["size"] > 50, "sv",
                           np.where(C["size"] > 0, "smallindel", "mnp"))
    # the allele travels with the hit now, so nothing about size is inferred
    C["size_inferred"] = False
    A = pd.read_csv(f"{OUT}/gwas_hits_allele_resolved.csv")
    nat = A.drop_duplicates(["chrom", "pos"]).set_index(["chrom", "pos"]).n_alleles_tested_at_pos
    C["n_alleles_at_pos"] = pd.MultiIndex.from_frame(C[["chrom", "pos"]]).map(nat)
    C["n_records_at_pos"] = C["n_alleles_at_pos"]
    C["pos_multiallelic"] = C.n_alleles_at_pos > 1
    C = C.sort_values(["n_gardens", "best_nlp"],
                      ascending=[False, False]).reset_index(drop=True)

    R.to_csv(f"{OUT}/gwas_pool_records.csv", index=False)
    C.to_csv(f"{OUT}/gwas_pool.csv", index=False)

    print(f"\n{len(R):,} significant hits -> {len(C):,} unique markers")
    print(C.vclass.value_counts(dropna=False).to_string())
    print("\nrecurrence over gardens:")
    print(C.n_gardens.value_counts().sort_index().to_string())
    print(f"\nMAC < 12 (would be dropped by the MAF-5% arm): "
          f"{int((C.mac < 12).sum()):,} / {len(C):,}")
    print(f"markers at a position carrying >1 tested allele (allele now resolved, "
          f"not guessed): {int(C.pos_multiallelic.fillna(False).sum()):,}")
    print(f"\nwrote {OUT}/gwas_pool.csv, {OUT}/gwas_pool_records.csv")


if __name__ == "__main__":
    main()

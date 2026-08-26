#!/usr/bin/env python
"""GWAS candidate pool: every Bonferroni-significant NON-SNP marker, 30 per-garden scans.

What goes in
------------
`r3_persite_gwas/results/gemma_gwas/persite_bonferroni_hits.csv` -- the 30 per-garden
GEMMA scans on the founder selection coefficient `s`, restricted to the non-SNP classes
(`nonsnp`, `sv`; `nonsnp` is smallindel+sv pooled, so a marker can appear in both).

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
The per-garden result arrays carry only chrom/pos, so ref_len/alt_len are recovered
from the arch3 panel meta. 2.14% of arch3 positions carry more than one biallelic
record (memory `panel-multiallelic-pos-key-trap`), so a position key can match the
wrong alt allele. Rather than hide that, the largest record at the position is taken
and `pos_multiallelic` flags every row where the position is ambiguous.

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

HITS = (f"{lib.GEA}/r3_persite_gwas/results/gemma_gwas/"
        "persite_bonferroni_hits.csv")
PANEL = "/global/scratch/projects/fc_moilab/tbellg/kmate/panel/arch3"
OUT = f"{HERE}/results"

CLASSES = ["nonsnp", "sv"]        # non-SNP layer only; the snp scan is not the pool


def panel_sizes(need: pd.DataFrame) -> pd.DataFrame:
    """chrom/pos -> ref_len/alt_len from the arch3 panel meta, with an ambiguity flag.

    Returns one row per (chrom,pos) present in `need`: the LARGEST record at that
    position, plus n_records_at_pos so multiallelic positions are visible.
    """
    # the GEMMA hit table writes chrom lowercase ("chr2"); lib.CHROMS is "Chr2"
    rows = []
    for ch in lib.CHROMS:
        key_ch = ch.lower()
        want = need.loc[need.chrom.str.lower() == key_ch, "pos"].astype(np.int64)
        if not len(want):
            continue
        f = f"{PANEL}/{key_ch}/var_pa_231_arch3_{key_ch}.meta.npz"
        if not os.path.exists(f):
            print(f"  MISSING {f}", flush=True)
            continue
        d = np.load(f, allow_pickle=True)
        pos = d["pos"]
        sel = np.isin(pos, want.to_numpy())
        if not sel.any():
            continue
        sub = pd.DataFrame({"chrom": key_ch, "pos": pos[sel],
                            "ref_len": d["ref_len"][sel], "alt_len": d["alt_len"][sel]})
        sub["size"] = (sub.alt_len - sub.ref_len).abs()
        n_at = sub.groupby("pos")["size"].size().rename("n_records_at_pos")
        best = sub.sort_values("size", ascending=False).groupby("pos").first()
        best = best.join(n_at).reset_index()
        best["chrom"] = key_ch
        rows.append(best)
        print(f"  {key_ch}: {len(best):,} of {len(set(want)):,} positions matched",
              flush=True)
    if not rows:
        raise SystemExit("no panel positions matched")
    M = pd.concat(rows, ignore_index=True)
    M["pos_multiallelic"] = M.n_records_at_pos > 1
    return M[["chrom", "pos", "ref_len", "alt_len", "size",
              "n_records_at_pos", "pos_multiallelic"]]


def main():
    os.makedirs(OUT, exist_ok=True)
    H = pd.read_csv(HITS)
    print(f"all Bonferroni hits: {len(H):,}")
    R = H[H.cls.isin(CLASSES)].copy()
    print(f"non-SNP classes only: {len(R):,} "
          f"({R.cls.value_counts().to_dict()})")

    key = ["chrom", "pos"]
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

    M = panel_sizes(C)
    n_before = len(C)
    C = C.merge(M, on=key, how="left")
    miss = int(C.ref_len.isna().sum())
    if miss:
        print(f"  WARNING: {miss} markers had no panel record (size left NaN)")
    assert len(C) == n_before, "panel merge changed row count"

    C["vclass"] = np.where(C["size"] > 50, "sv",
                           np.where(C["size"] > 0, "smallindel", "mnp"))
    # `size` is the LARGEST record at the position. Where the position carries more
    # than one record that is a GUESS, not the associated allele -- so vclass, any
    # frameshift call downstream, and any join keyed on size are all unreliable for
    # these rows. Flag them explicitly rather than let the guess look like data.
    C["size_inferred"] = C.pos_multiallelic.fillna(False)
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
    print(f"multiallelic positions (alt allele ambiguous): "
          f"{int(C.pos_multiallelic.fillna(False).sum()):,}")
    print(f"\nwrote {OUT}/gwas_pool.csv, {OUT}/gwas_pool_records.csv")


if __name__ == "__main__":
    main()

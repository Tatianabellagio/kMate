#!/usr/bin/env python
"""GEA candidate pool: every Bonferroni-significant NON-SNP record, 22 climate axes.

What goes in
------------
Union over cls in {sv, smallindel, nonsnp} x axis in {bio1..19, pc1, pc2, pc3} of the
raw per-record LFMM p-values in `wza_in_clq09_tile`, at per-(class x axis) genome-wide
Bonferroni (0.05 / n_records, MAF > 0.05).

Two deliberate choices, both set with the user 2026-08-26:

1. **All Bonferroni records, not one lead per block.** The older
   `dissection/screen_sig_blocks.py` kept only each significant block's strongest
   variant. That silently drops a genic variant whenever a marginally stronger
   intergenic neighbour shares its block -- the opposite of what the CARK lesson asks
   for (judge every candidate at its OWN position). Note the `nonsnp` class scan is
   smallindel+sv pooled, so a variant can be significant in both `sv` and `nonsnp`;
   `sig_classes` records which.

2. **Non-SNP hits only, SNP-shadowed or not.** The pool is the non-SNP layer. Whether a
   Bonferroni SNP sits nearby is carried as an annotation (`snp_cosig_2kb`, added by
   `dissection/augment_snp_cosig.py`) rather than used as a filter, so SNP-shared and
   SNP-unique candidates both stay in and can be split later.

Recurrence and inflation
------------------------
`n_axes` overcounts, because the bioclim axes are correlated (see `axis_clusters.py`).
So the pool also carries `n_clusters` -- recurrence over the 7 empirical correlation
clusters -- which is the number the convergence filter should use.

Genomic-control lambda for the axis is attached per record from
`raw_block_significance_lfmm.csv` (per class x axis). `min_lam` is the lowest lambda
among the axes a variant is significant on, and `has_clean_axis` marks variants with at
least one hit on an axis with lambda < LAM_CLEAN. bio15 (2.67) and pc3 (3.13) together
contribute ~58% of all Bonferroni records, so without this flag the pool is mostly the
two worst-calibrated axes.

**Caveat that travels with every row:** raw uncalibrated LFMM p. This is a candidate
net, not a calibrated hit list. See memory `raw-lfmm-over-wza-decision`.

Outputs -> results/
  gea_pool.csv            one row per unique variant
  gea_pool_records.csv    one row per (variant x class x axis) significant record

env: kmate.  Run on a compute node (reads ~2 GB of CSV).
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                    # noqa: E402
import axis_clusters as ac                                    # noqa: E402

WZAIN = (f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/"
         "wza_in_clq09_tile")
LAMTAB = (f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/"
          "raw_block_significance_lfmm.csv")
OUT = f"{HERE}/results"

AXES = ac.AXES                                   # bio1..19 + pc1,pc2,pc3
CLASSES = ["sv", "smallindel", "nonsnp"]         # non-SNP layer only
MAF_MIN = 0.05                                   # the LFMM record filter
LAM_CLEAN = 2.0                                  # "not one of the badly inflated axes"


def lambda_table() -> dict[tuple[str, str], float]:
    d = pd.read_csv(LAMTAB)
    d = d[d.model == "lfmm"]
    return {(r.cls, r.axis): float(r.lam) for r in d.itertuples()}


def sig_records() -> pd.DataFrame:
    """Every Bonferroni-significant record, one row per (variant, class, axis)."""
    lam = lambda_table()
    frames = []
    for cls in CLASSES:
        for ax in AXES:
            f = f"{WZAIN}/lfmm_{cls}_gen9_{ax}.csv"
            if not os.path.exists(f):
                print(f"  MISSING {os.path.basename(f)}", flush=True)
                continue
            d = pd.read_csv(f)
            d = d[(d.MAF > MAF_MIN) & d.block.notna() & (d.block != "")]
            if not len(d):
                continue
            thr = 0.05 / len(d)                          # genome-wide, this class x axis
            s = d[d.pval < thr].copy()
            print(f"  {cls:11s} {ax:5s}  n={len(d):>9,}  bonf<{thr:.3e}  "
                  f"sig={len(s):>5,}", flush=True)
            if not len(s):
                continue
            s["cls"] = cls
            s["axis"] = ax
            s["nlp"] = -np.log10(s.pval.clip(lower=1e-300))
            s["lam"] = lam.get((cls, ax), np.nan)
            frames.append(s[["cls", "axis", "block", "chrom", "pos", "ref_len",
                             "alt_len", "MAF", "nlp", "lam"]])
    if not frames:
        raise SystemExit("no significant records found")
    return pd.concat(frames, ignore_index=True)


def collapse(R: pd.DataFrame, cl: dict[str, str]) -> pd.DataFrame:
    """One row per unique variant, with recurrence over axes AND axis clusters."""
    R = R.copy()
    R["cluster"] = R.axis.map(cl)
    key = ["chrom", "pos", "ref_len", "alt_len"]
    R = R.sort_values("nlp", ascending=False)
    g = R.groupby(key, sort=False)

    out = g.agg(
        MAF=("MAF", "first"),
        best_axis=("axis", "first"),
        best_nlp=("nlp", "max"),
        best_lam=("lam", "first"),
        min_lam=("lam", "min"),
        block=("block", "first"),
        n_records=("axis", "size"),
    ).reset_index()

    out["n_axes"] = g["axis"].nunique().values
    out["sig_axes"] = g["axis"].apply(lambda s: ",".join(dict.fromkeys(s))).values
    out["n_clusters"] = g["cluster"].nunique().values
    out["sig_clusters"] = g["cluster"].apply(
        lambda s: ",".join(dict.fromkeys(s))).values
    out["sig_classes"] = g["cls"].apply(lambda s: ",".join(sorted(set(s)))).values

    # cleanest (least inflated) axis this variant is significant on
    idx_clean = g["lam"].idxmin()
    out["cleanest_axis"] = R.loc[idx_clean.values, "axis"].values
    out["has_clean_axis"] = out.min_lam < LAM_CLEAN

    out["size"] = (out.alt_len - out.ref_len).abs()
    # NB `dissection/screen_sig_blocks.py` labels size==0 as "snp". That is wrong here:
    # these scans contain no SNPs, and every size==0 record has ref_len==alt_len>=2,
    # i.e. a multi-nucleotide substitution. Label them mnp, not snp.
    out["vclass"] = np.where(out["size"] > 50, "sv",
                             np.where(out["size"] > 0, "smallindel", "mnp"))
    return out.sort_values(["n_clusters", "best_nlp"],
                           ascending=[False, False]).reset_index(drop=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    cl = ac.clusters()
    print(f"axis clusters ({len(set(cl.values()))}): "
          f"{sorted(set(cl.values()))}\n")

    R = sig_records()
    C = collapse(R, cl)

    R.to_csv(f"{OUT}/gea_pool_records.csv", index=False)
    C.to_csv(f"{OUT}/gea_pool.csv", index=False)

    print(f"\n{len(R):,} significant records -> {len(C):,} unique variants")
    print(C.vclass.value_counts().to_string())
    print("\nrecurrence over axis clusters:")
    print(C.n_clusters.value_counts().sort_index().to_string())
    print(f"\nwith >=1 axis at lambda < {LAM_CLEAN}: "
          f"{int(C.has_clean_axis.sum()):,} / {len(C):,}")
    print(f"\nwrote {OUT}/gea_pool.csv, {OUT}/gea_pool_records.csv")


if __name__ == "__main__":
    main()

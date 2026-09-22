#!/usr/bin/env python
"""Does the promoter result depend on how "promoter" is defined? Same enrichment, four rules.

Rules (gea_promoter_strict.py): the classifier's own 1 kb overlap, the same window but requiring
the variant's whole span to be CONTAINED, a 500 bp overlap window, and 500 bp contained. Same
statistic as the notebook's section 2: odds of being a GEA hit for promoter vs non-promoter
variants of the same class and type, Mantel-Haenszel over size x MAF strata, 95% CI from
resampling LD blocks. env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gea_region_by_axis as GR                                  # noqa: E402
OUT = f"{HERE}/results"
K4 = ["chrom", "pos", "ref_len", "alt_len"]
RULES = ["tier_1kb", "tier_contained_1kb", "tier_overlap_500", "tier_contained_500"]

U = pd.read_csv(f"{OUT}/gea_hit_landscape.csv.gz")
U = U.sort_values(["hit", "min_p"], ascending=[False, True]).drop_duplicates(K4)
U = U[U.kind.isin(["deletion", "insertion"])]
R = pd.read_csv(f"{OUT}/gea_universe_regions_strict.csv.gz",
                usecols=["chrom", "pos", "ref_len"] + [r for r in RULES if r != "tier_1kb"])
U = U.merge(R, on=["chrom", "pos", "ref_len"], how="left")
U["stratum"] = (pd.cut(U["size"], GR.SIZE_BINS, labels=False).astype(str) + "|"
                + pd.cut(U.MAF, GR.MAF_BINS, labels=False, include_lowest=True).astype(str))
rng = np.random.default_rng(1); rows = []
for cls in ("sv", "smallindel"):
    for kind in ("deletion", "insertion"):
        d = U[(U.cls == cls) & (U.kind == kind)]
        for rule in RULES:
            e = d.assign(_x=d[rule].eq("3_promoter"))
            o, lo, hi, p = GR.or_with_ci(e, "hit", "_x", rng)
            rows.append(dict(cls=cls, kind=kind, rule=rule.replace("tier_", ""),
                             promoter_tested=int(e._x.sum()), promoter_hits=int((e.hit & e._x).sum()),
                             hit_share=e[e.hit]._x.mean(), tested_share=e._x.mean(),
                             OR=o, lo=lo, hi=hi, p=p))
T = pd.DataFrame(rows)
T.to_csv(f"{OUT}/gea_promoter_strict_enrichment.csv", index=False)
pd.set_option("display.width", 200)
print(T.assign(**{c: T[c].round(3) for c in ["hit_share", "tested_share", "OR", "lo", "hi", "p"]}).to_string(index=False))

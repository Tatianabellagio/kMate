#!/usr/bin/env python
"""Flag -- never drop -- records whose pool-AF vector is byte-identical to a distant record.

af_duplicate_vectors.py found 5,589 testable non-SNP records (0.74%) sharing an identical
gen9 pool-AF vector with a record >50 kb away. Their climate signal is not locus-specific,
but the record, its mechanism and its expression evidence are still data. So this is an
annotation, like `snp_blind` / `tag_untestable`: it travels with the candidate and is read
alongside the other lines of evidence, it does not remove anything.

Adds to results/candidate_evidence.csv and results/mode_of_action_scan.csv:
  af_shared_vector   True if the record's pool-AF vector is shared with a record >50 kb away
  af_group_size      how many records share it
  af_group_span_mb   genomic span of the group (or 'multi-chrom')
"""
import os
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = f"{HERE}/results"

G = pd.read_csv(f"{OUT}/af_duplicate_groups.csv")
g = G.groupby("h").agg(n=("row", "size"), nchr=("chrom", "nunique"),
                       span=("pos", lambda x: x.max() - x.min()))
G = G.join(g, on="h")
info = {r.row: (r.n, "multi-chrom" if r.nchr > 1 else round(r.span / 1e6, 2))
        for r in G.itertuples()}
for f in ("candidate_evidence.csv", "mode_of_action_scan.csv"):
    D = pd.read_csv(f"{OUT}/{f}")
    D["af_shared_vector"] = D.store_row.isin(info)
    D["af_group_size"] = D.store_row.map(lambda r: info.get(r, (1, ""))[0])
    D["af_group_span_mb"] = D.store_row.map(lambda r: info.get(r, (1, ""))[1])
    D.to_csv(f"{OUT}/{f}", index=False)
    print(f"{f}: {int(D.af_shared_vector.sum())} of {len(D)} flagged")

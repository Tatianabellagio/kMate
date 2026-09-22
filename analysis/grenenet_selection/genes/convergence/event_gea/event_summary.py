#!/usr/bin/env python
"""Flags each tested SV event with how the ALLELE-level GEA saw it, for the notebook's event section.

allele_tested  at least one member record cleared MAF 0.05 as a record (so the record-level GEA
               could test the event at all); False = newly testable, only as an event
allele_hit     at least one member record was a record-level GEA hit
Adds both to event_gea/event_sv_hits.csv. env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
GEA = os.path.abspath(os.path.join(CONV, "..", ".."))
sys.path.insert(0, GEA)
import lib                                                       # noqa: E402
K4 = ["chrom", "pos", "ref_len", "alt_len"]
U = np.load(f"{GEA}/r1_sv_negative_selection/results/sv_adaptive/nonsnp_event_units.npz")
ptr, ind = U["indptr"], U["indices"]
idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
R = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"], "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
L = pd.read_csv(f"{CONV}/results/gea_hit_landscape.csv.gz", usecols=K4 + ["hit", "cls"])
L = L[L.cls == "sv"]
tested = set(map(tuple, L[K4].values)); hit = set(map(tuple, L[L.hit][K4].values))
E = pd.read_csv(f"{HERE}/event_sv_hits.csv")
mem = [[tuple(x) for x in R.iloc[ind[ptr[u]:ptr[u + 1]]][K4].values] for u in E.unit]
E["allele_tested"] = [any(m in tested for m in ms) for ms in mem]
E["allele_hit"] = [any(m in hit for m in ms) for ms in mem]
E.to_csv(f"{HERE}/event_sv_hits.csv", index=False)
print(f"{len(E):,} events; newly testable {int((~E.allele_tested).sum()):,}; event hits {int(E.hit.sum())} "
      f"({int((E.hit & ~E.allele_tested).sum())} newly testable); allele-hit events {int(E.allele_hit.sum())}, "
      f"still hits {int((E.hit & E.allele_hit).sum())}")

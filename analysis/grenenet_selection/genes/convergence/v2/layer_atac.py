#!/usr/bin/env python
"""v2 layer -- open chromatin (ATAC) overlap of each GEA unit, keyed on `unit`.

Same overlap rule as atac_overlap.py (`peak_overlap`: the unit's representative-record REF footprint [pos, pos+ref_len-1] against the
multi-tissue ATAC union; flower / leaf / root / shoot; n_tissues = tissues in which the overlapped peak was called). For an SV event the
footprint is its representative member (most carriers). Also reports the overlap rate over EVERYTHING the GEA tested (tile_cr50 nonsnp, 623,764
units) as the background, by class -- a context number, not a significance test.
Output: results/candidate_list_v2_events_cr50/v2_layer_atac.csv   (one row per pooled unit)
env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
sys.path.insert(0, CONV); sys.path.insert(0, os.path.abspath(os.path.join(CONV, "..", "..")))
import atac_overlap as AO                                        # noqa: E402
import lib                                                       # noqa: E402
V2 = os.path.join(CONV, "results", "candidate_list_v2_events_cr50")

P = pd.read_csv(f"{V2}/gea_pool.csv")[["unit", "chrom", "pos", "ref_len", "alt_len", "vclass"]]
peaks = AO.load_peaks()
ov = AO.peak_overlap(peaks, P.chrom.to_numpy(str), P.pos.to_numpy(), P.ref_len.to_numpy())
L = pd.concat([P[["unit"]].reset_index(drop=True), ov.reset_index(drop=True)], axis=1)
L.to_csv(f"{V2}/v2_layer_atac.csv", index=False)

# background: every unit the production GEA tested (nonsnp scan, bio1 file lists them all)
B = pd.read_csv(f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/tile_cr50/lfmm_nonsnp_gen9_bio1.csv",
                usecols=["unit", "chrom", "pos", "ref_len", "alt_len"])
B["vclass"] = np.where((B.alt_len - B.ref_len).abs() > 50, "sv", np.where(B.alt_len == B.ref_len, "mnp", "smallindel"))
Bo = AO.peak_overlap(peaks, B.chrom.to_numpy(str), B.pos.to_numpy(), B.ref_len.to_numpy())
B["in_atac"] = Bo.in_atac.values
rows = []
for c in ("sv", "smallindel", "mnp"):
    cand = L[P.vclass.values == c]; bg = B[B.vclass == c]
    rows.append(dict(vclass=c, n_candidates=len(cand), cand_in_atac_pct=100 * cand.in_atac.mean(), n_tested=len(bg),
                     tested_in_atac_pct=100 * bg.in_atac.mean()))
print(pd.DataFrame(rows).round(1).to_string(index=False))
print(f"\nwrote {V2}/v2_layer_atac.csv ({len(L):,} units; in a peak: {int(L.in_atac.sum()):,})")

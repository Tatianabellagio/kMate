# Top convergence candidates — verdicts after locus dissection

The three genes hit by **both** scans, plus GSH1 (same locus as CRK18, three themes),
put through `dissection/plot_locus_combined.py` via `plot_top_loci.py`. Figures in
`plots/loci/`.

The number that decides it is the printed **LD-confirm**: founder r² between the lead
variant and the candidate gene's own variants, against the best r² the lead achieves
anywhere local. If the lead is better linked to something else than to the gene it was
attributed to, the attribution is CARK-type detachment and the candidate is not about
that gene.

| gene | theme | r²(lead ↔ gene) max / median | best local r² | carriers | verdict |
|---|---|---|---|---|---|
| **CRK18** (AT4G23260) | stress | **0.874** / 0.064 | 0.933 | 15 | **ROBUST** — the lead is nearly as well linked to CRK18 as to anything in the window. Cross-scan, frameshift, and the locus is the WZA honest-null survivor. The best candidate in the analysis. |
| **GSH1** (AT4G23100) | temperature, flowering, stress | 0.452 / 0.005 | 0.463 | 11 | **MODERATE** — lead is not better linked elsewhere (0.452 vs 0.463 local), so the attribution is not detached, but r² is middling and the median is ~0. Same locus as CRK18. |
| **EPFL5** (AT3G22820) | stress | 0.391 / 0.023 | 0.700 | 22 | **WEAK** — the lead is nearly twice as well linked to something else (0.700) as to EPFL5. 4 gardens is real recurrence, but the gene assignment is suspect. |
| **BT4** (AT5G67480) | stress | **0.069** / 0.019 | 0.163 | 11 | **FRAGILE — do not carry forward as a BT4 result.** The lead is essentially unlinked to BT4 (r² 0.069) and poorly linked to anything local. Textbook CARK-type detachment: cross-scan support does not survive asking *which gene*. |

## What this changes

Cross-scan agreement (`L_cross`) says a **locus** is seen by two independent scans. It
does not say **which gene** at that locus is responsible — that is a separate question,
and it is what the LD-confirm answers. Of the three cross-scan genes, one survives
cleanly (CRK18), one is suspect (EPFL5) and one fails (BT4).

This is the same triage rule the GEA shortlist arrived at independently
(`dissection/results/loci/CANDIDATE_VERDICTS.md`): robustness = recurrence + founder
frequency + **LD to the gene**, not p-value. Adding cross-scan evidence strengthens the
first term; it does nothing for the third.

## Caveats specific to these four

- **CRK18's frameshift call rests on an inferred allele.** The frameshift variant
  (Chr4:12,169,155) is GWAS-derived, and that position carries 8 panel records with
  length differences 0,1,1,1,1,1,1,2; `build_gwas_pool.py` took the largest. Seven of
  the eight would still shift the frame, so "frameshift" holds, but "2 bp" does not.
  The variant plotted as the lead is the *GEA* one at Chr4:12,170,418 (1 bp, promoter,
  bio6) — a different variant in the same gene.
- All four sit in loci contributing many candidate genes (CRK18/GSH1 locus: 19;
  EPFL5: 20; BT4: 12), so gene counts at these loci are LD pileups, not independent
  hits.
- The figures inherit `plot_locus_combined.py`'s layout, which crowds the annotation
  against the gene track at this window width. Cosmetic, not a data problem.

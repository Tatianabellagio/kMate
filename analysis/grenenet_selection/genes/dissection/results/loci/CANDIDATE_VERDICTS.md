# Per-candidate locus dissection — verdicts

Each candidate got the cam5/CARK-style combined figure (`loci/<sym>_combined.png`):
LFMM Manhattan (Δp hot−cold colour) + founder-haplotype panel (cold→hot home bio1)
+ founder-LD triangle, with a printed LD-confirm of the lead variant vs the gene's
own variants. Sorted robust → fragile.

| gene | variant | region | axis (n_axes) | nlp | founder carriers | r²(lead↔gene) | verdict |
|---|---|---|---|---|---|---|---|
| **EMB1241** (GrpE co-chaperone) | 2 bp indel | 3′UTR | bio15 (**7**) | 11.9 | 27 | **0.92** | **ROBUST — top pick**: most replicated, strongly tags gene haplotype, hot-assoc |
| **GPX6** (glutathione peroxidase) | 1.2 kb SV del | promoter | pc1 (**6**) | 10.2 | 27 | 0.24 | **ROBUST**: recurrent; promoter cis (LD to gene modest but physically placed) |
| **TPT** (photosynth. acclimation) | 428 bp SV ins | promoter | bio1 (2) | 7.5 | **106** | 0.57 | **ROBUST**: common variant, hot-assoc, well-placed |
| **AGL100** (MADS-box TF) | 373 bp SV del | promoter | bio15 (2) | 10.0 | 10 | 0.70 | **GOOD**: well-linked to gene, hot-assoc (SNP-shadowed) |
| **AIPP3/RVR1** (Repressor of Vernalization 1) | 4 bp indel | promoter | bio6 (4) | 10.9 | 11 | 0.43 | **MODERATE**: recurrent, moderate LD; in a Chr4 candidate cluster w/ GPX6 |
| **AFP2** (ABI5-binding, ABA) | 1 bp indel | promoter | bio15 (1) | 11.5 | 172 | 0.47 | **MODERATE**: common but single-axis, cold-assoc, low-diversity region |
| **AT3G22142** (2S-albumin/LTP) | 451 bp SV del | **CDS** | bio15 (2) | 11.6 | **3** | 0.24 | **FRAGILE**: rare; the structural-SV of interest but signal rests on 3 lineages |
| **FAMA** (bHLH, stomatal) | 1 bp indel | promoter | bio15 (1) | 13.2 | 17 | 0.30 | **FRAGILE**: highest p but an isolated spike, no LD support, single-axis |
| **HSP90-4** (heat-shock 90) | 4 bp indel | promoter | bio15 (1) | 9.8 | 6 | 0.09 | **FRAGILE**: rare, isolated, essentially unlinked to gene |
| **FRL2** (FRIGIDA-like 2) | 1 bp del | 3′UTR | bio15 (1) | 9.5 | 6 | 0.19 | **FRAGILE**: rare, isolated spike, low gene-LD |

## Triage rule that emerged
Robust vs fragile separates cleanly on three axes-independent features, NOT on p-value:
- **recurrence** (`n_axes` ≥ 2–7),
- **founder frequency** (common vs 3–6 carriers),
- **LD to the gene's own haplotype** (r² high vs ~0).

The fragile four (FAMA, HSP90-4, FRL2, AT3G22142) are all rare (3–17 carriers),
single-axis, isolated ultra-high-p spikes with near-zero gene-LD — the classic
rare-variant/structure-inflation profile. FAMA's nlp 13.2 is the cleanest warning
that **p-value alone is not evidence of a real, gene-linked locus here**.

## Where to go next
Top picks for follow-up: **EMB1241, GPX6, TPT, AGL100** (+ AIPP3/RVR1). All are
functionally placed (UTR/promoter), climate-recurrent or common, and gene-linked.
The Chr4 ~6.99–7.01 Mb cluster (AIPP3/RVR1 + GPX6 + CRK33/PUX8 nearby) may be one
adaptive haplotype worth dissecting as a block.

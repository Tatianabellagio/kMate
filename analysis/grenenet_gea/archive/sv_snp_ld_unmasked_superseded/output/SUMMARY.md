# SV–SNP tagging (LD) — how much of our SVs are captured by SNPs

Per-SV maximum r² with any SNP within ±50 kb, across the 231 founders, genome-wide
(226,468 true SVs, >50 bp). Two SNP panels:

- **our panel** — SNPs from the same `merged_231` assembly pangenome as the SVs (V_pa)
- **short-read** — GrENE-net short-read SNP panel `greneNet_final_v1.1` (Xing Wu)

Proper 3 tiers (well-tagged / partial / truly independent):

| Panel | **well-tagged** (r²≥0.8) | **partial** (0.2–0.8) | **independent** (r²<0.2) |
|---|---|---|---|
| our panel (assembly SNPs) | 85.8% | 13.7% | **0.4%** |
| short-read SNPs (GrENE-net) | 47.4% | 45.4% | **7.3%** |

**Headline:** with a *complete* assembly-based pangenome SNP set, SVs are largely
SNP-tagged (~86% well-tagged). With a *conventional short-read* SNP panel (what GrENE-net
and most studies use), only **~47% of SVs are well-tagged (r²≥0.8)** — the rest split into
**~45% partially tagged** (a SNP in moderate LD → SNP-GEA detects them with REDUCED power)
and **~7% truly independent** (no SNP in LD → entirely SNP-invisible).

Note: "independent" = best r² < 0.2 (only ~7%, NOT 53%). The "53% not well-tagged"
figure mixes partial + independent; most of it is partial. The defensible claims are:
(i) only ~47% well-tagged by short-read SNPs vs ~86% by the full pangenome; (ii) ~7% of
SVs have no SNP in LD at all. Consistent with Yan et al. 2021 (54% of adaptive human SVs
in strong SNP LD).

**Larger SVs are even less SNP-tagged by short reads** (>5 kb: ~43%) — they sit in
repeat-rich regions where short-read SNP calling fails — whereas our assembly panel tags
them uniformly (~85%). So the SVs most likely to be missed by SNP studies are the large ones.

Figure: `sv_snp_tagging.png` (ECDF + %tagged bars + tagging-vs-size).
Per-chrom per-SV r²: `sv_snp_ld_{panel,shortread}_Chr*.npz` (best_r2, best_snp_pos, sv_size).
Code: `analysis/grenenet_gea/build_sv_snp_ld.py`, `extract_shortread_geno.py`.

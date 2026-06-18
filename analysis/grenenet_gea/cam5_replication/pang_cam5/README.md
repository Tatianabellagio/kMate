# CAM5 / AT2G27030 — combined GEA + pangenome-haplotype figures

Final, talk-ready figures for **CAM5** (Chr2:11,531,800–11,534,333, + strand): the
3′ climate-GEA signal is a **mixed SNP + small-indel haplotype** (the 3 Bonferroni-
significant variants are in **r²≈1** — one tight haplotype), carried by rare,
warm-origin founders. (All exploratory plots from earlier sessions — tube-map,
variation tracks, bandage/vg graphs, standalone gggenomes/hap-zoom — were removed
2026-06-17; this dir is now just the final set + its pipeline.)

## Final figures
Each is a stacked, x-aligned LocusZoom-style panel (PNG + vector PDF):
1. **LFMM Manhattan** — −log10 p; colour = **Δp (p_final − p0, SEEDMIX)**; shape =
   SNP circle / indel square; **gene-level Bonferroni (0.05/71)** line; 3 sig
   variants marked with grey droplines. Frameless.
2. **Unique-haplotype pangenome** — gene model (isoform AT2G27030.3), per-founder
   variants (grey), √-scaled **abundance bars** coloured by mean home bio1 on a
   blue→grey→red diverging scale; reference haplotype labelled.
3. **Founder-LD triangle** — r² from the 231-founder genotypes, rank-spaced, white→
   purple; 3 vertical grey connectors from the sig variants.

- `cam5_combined_region{,_squares}.png/.pdf` — zoomed to the 3′ region (17 unique
  haplotypes; `_squares` = variants as SNP-circle/indel-square instead of type/size
  glyphs).
- `cam5_combined_manhattan_hap{,_squares}.png/.pdf` — whole gene (109 haplotypes).

## Pipeline (env: `gggenomes` mamba env)
1. **Data prep** (basic env / pysam):
   - `prep_gggenomes.py` → `gg_hapseqs/gg_hapfeats/gg_hapmeta/gg_model.tsv`
     (whole-gene haplotype collapse + gene model).
   - `prep_region_haps.py [LO HI]` → `gg_rhapseqs/gg_rhapfeats/gg_region.tsv`
     (3′-region collapse; default window 11,533,810–11,534,333, i.e. **starting after
     the near-fixed 5′ block**, see note).
   - `compute_ld.py {region,gene} [LO HI]` → `gg_ld_{region,gene}.tsv` (founder r²).
2. **Plot** (`~/miniforge3/envs/gggenomes/bin/Rscript`):
   - `Rscript plot_combined_region.R [squares]`
   - `Rscript plot_combined_manhattan_hap.R [squares]`

## Source data (kept)
- `cam5_region.vcf.gz` — 231-founder genotypes in the region (panel VCF).
- `cam5_features.tsv` — TAIR10 isoform/exon model.
- `cam5_variants_gen9.tsv` — per-variant GEA stats: `tau, kendall_p, cls, lfmm_p, dp`
  (dp = p_final − p0 from `results/grenenet_gea/group_means.npz`; p0 = SEEDMIX).

## Notes worth remembering
- **Bonferroni is gene-level** (0.05 / 71 tested variants in CAM5), NOT genome-wide.
  Best variant p≈3.6e-4, far from genome-wide 5e-8.
- **Panel vs Manhattan variant counts differ on purpose**: founder panel = 169
  variants in the gene; only ~71 are GEA-tested (Manhattan). The ~93 untested are
  rare (median founder MAF 0.9%, ~half singletons/doubletons) → no testable cohort
  AF trajectory, so no LFMM p. The haplotype panel shows the full genotype panel.
- **Near-fixed 5′ block** (11,533,741–11,533,796: 4 SNP + an 11 bp del) is carried by
  ~81% of founders because TAIR10 holds the MINOR allele there (polarization, not
  divergence). The region figure starts after it so the common haplotype collapses
  into one reference row (n=172).
- LD triangle is **rank-spaced** (equal slot per variant) → clean triangle but not
  bp-exact vs the panels above. High-LD r²≈1 between two *non-adjacent* variants shows
  up *deep* in the triangle, not at the apex between them.

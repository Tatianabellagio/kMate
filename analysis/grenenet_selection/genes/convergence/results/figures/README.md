# Candidate figures, one tree

Built by `organize_figures.py` from the review tables; re-run it after any re-render.
The raw render output stays in `../plots/` (screen/, screen_ownaxis/, loci/) -- this tree
is the curated view, and every file here matches the representative variant in the
CURRENT tables (allele-resolved for the GWAS set, GEA-significant records only).

    01_top20/       80 files -- the 20 genes to look at first, rank-ordered
    gea_own_axis/   844 files -- 211 GEA genes judged on their own climate axis
    gwas/           180 files -- 45 per-garden GWAS genes
    round1/         328 files -- the round-1 89, minus the 7 withdrawn on 2026-09-17
    INDEX.csv       gene, set, rank, verdict, grade, axis, store_row, file names

Two figure types per gene:

  `<sym>__grid.png`   one panel per garden, ordered by the climate axis; dots = pools,
                      bold line = garden mean, dashed = founding frequency, red frame +
                      star = Bonferroni-significant in that garden's per-garden GWAS.
  `<sym>__locus.png`  top: -log10 p around the gene (lead variant outlined, gene shaded,
                      dashed lines = SNP / SV Bonferroni); middle: founder genotypes
                      sorted by home temperature; bottom: founder LD.

Verdicts: A pattern + co-signal locus, B pattern + lone SNP-blind variant, C pattern but
the lead is a passenger, D partial, E none. All are single-pass visual calls on raw
uncalibrated LFMM / GEMMA p-values.

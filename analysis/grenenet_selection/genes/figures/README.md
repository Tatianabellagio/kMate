# Candidate figures -- one flat folder

Every figure for every gene under `genes/`, in one place. Built by
`convergence/organize_figures.py`; re-run it after any re-render. Entries are hard links,
so this costs no extra disk and the raw render output stays where the scripts put it
(`convergence/results/plots/{screen,screen_ownaxis,loci,atac}`,
`dissection/results/plots`). Every reviewed file matches the representative variant in the
CURRENT tables (allele-resolved for the GWAS set, GEA-significant records only).

FLAT -- one directory, 744 PNGs and 744 PDFs, no subfolders. A gene's figures
share its symbol prefix, so they sort together and `ls | grep GPX6` finds all of them.
Grouping that used to be folders (top-20 rank, which scan, review round) is in the
`set` and `rank` columns of INDEX.csv instead.

Three figure types per gene:

  `<sym>__grid.png`   one panel per garden, ordered by the climate axis; dots = pools,
                      bold line = garden mean, dashed = founding frequency, red frame +
                      star = Bonferroni-significant in that garden's per-garden GWAS.
                      This is the allele-frequency-over-time view.
  `<sym>__locus.png`  top: -log10 p around the gene (lead variant outlined, gene shaded,
                      dashed lines = SNP / SV Bonferroni); middle: founder genotypes
                      sorted by home temperature; bottom: founder LD.
  `<sym>__atac.png`   functional tracks from `plot_atac.py`: GEA Manhattan, TFBS turnover,
                      per-tissue open chromatin, gene models on a shared axis.

INDEX.csv -- gene, sym, set, rank, verdict, grade, axis, store_row, grid, locus, atac.

    set=GEA_ownaxis  211 GEA genes judged on their own climate axis
    set=GWAS          45 per-garden GWAS genes
    set=round1        82 round-1 genes, the 7 withdrawn on 2026-09-17 removed
    set=unreviewed    rendered but carried by no review table: genes withdrawn on audit
                      (AT5G44220), ones that failed the LD-confirm (BT4), the
                      dissection-era ROBUST calls (EMB1241, GSH1), the CARK figures
    set=atac          functional-track figures
    rank=1..20        the genes to look at first

Files whose names are not `<sym>__*` are `set=unreviewed`: they keep their original render
names because they are not grid/locus review pairs. The coverage check this tree must pass
is that every gene with a PNG in any source dir has a row in INDEX.csv -- `unreviewed` is
what makes that true, without it ~90 genes existed on disk and nowhere here.

Verdicts: A pattern + co-signal locus, B pattern + lone SNP-blind variant, C pattern but
the lead is a passenger, D partial, E none. All are single-pass visual calls on raw
uncalibrated LFMM / GEMMA p-values.

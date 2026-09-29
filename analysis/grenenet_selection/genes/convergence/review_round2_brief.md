# Round-2 visual review brief

Reviewer instructions for `results/screen_visual_review_round2.csv` (268 genes: 223 `GEA_ownaxis`,
45 `GWAS`). Same reading as round 1 (`results/screen_visual_review.csv`, use its `notes`,
`V_ownaxis`, `locus_check` columns as calibration), with the criteria split by scan.

## Figures (in `results/plots/screen/`)

- `<sym>_garden_trajectories_ownaxis.png`: one panel per garden, gardens ordered by the climate
  axis in the `axis` column (value printed in each panel). Dots = pools, bold line = garden mean,
  dashed = p0, red frame and star = Bonferroni-significant in that garden's per-site GWAS.
- `<sym>_combined.png` (GEA) or `<sym>_combined_gwas.png` (GWAS): top = -log10 p Manhattan
  (the lead variant is the large outlined marker on the red dashed line, gene shaded, dashed lines
  = SNP/SV Bonferroni), middle = founder genotypes sorted by home temperature, bottom = founder LD.
- If `sym` collides with a round-1 gene the file is `<sym>_<target_gene>_...`.
- A missing locus figure means the gene has no coordinates; say so and judge the grid alone.

## GEA_ownaxis: pattern on the axis the hit was found on

V_grid: is there a monotone or threshold relationship between the axis and the change from p0
(rise at one end, fall or flat at the other)? Name the gardens that move and their values.
Grade STRONG / GOOD / MODERATE / WEAK.
- Garden-4 founder-sweep pattern (garden 4 rises strongly, 43/45/32 modestly, everything else
  flat) is WEAK regardless of axis. `r_garden4_set` > 0.6 is the numeric flag; say whether the
  picture is broader than garden 4.
- pc3 is the most inflated axis (lambda ~3.2); require a visibly clean gradient to grade above
  MODERATE.
- A pattern carried by 1-2 gardens is at most MODERATE.

## GWAS: per-garden significance, climate is context only

Rare alleles cannot give a clean climate gradient, and that is fine. Judge:
- how many gardens are starred (`gwas_n_gardens`), and whether those gardens actually show
  the allele rising (or falling) consistently over generations, not a single-pool spike;
- whether starred gardens share a climate (local adaptation) or span climates (global
  advantage). Say which.
V_grid grades: STRONG = >= 3 starred gardens with a consistent rise; GOOD = 2 starred with clean
rise, or >= 3 with some noise; MODERATE = 2 starred, noisy; WEAK = 1 starred garden, or the stars
do not match visible movement.
- Also note gene interest from `protein` (stress, climate, flowering, ion, signalling, etc).

## V_locus (both sets)

Is the lead the local peak or a passenger? Co-signal (SNPs/indels at the same place) or a lone
SNP-blind marker? Inside the gene's LD block? Neighbouring genes that could carry the signal?
Many founder carriers or few? One sentence to three.

## verdict (both sets)

- `A`: pattern + clean co-signal locus attributable to this gene
- `B`: pattern, lone SNP-blind variant with good placement
- `C`: pattern, but lead is a passenger or attribution is unclear
- `D`: partial pattern
- `E`: no pattern, or garden-4 sweep only

GWAS genes with >= 2 starred gardens and a real rise in them, whose gene is plausibly relevant,
should be `A`-`C` even without a climate gradient.

## Output

One CSV per batch in `results/review_round2_parts/` with columns
`target_gene,sym,V_grid,V_locus,verdict,grade,notes` (grade = the STRONG..WEAK word from V_grid).
Quote every field.

# convergence — one candidate pool from both scans, ranked by lines of evidence

Built 2026-08-26. The third layer under `genes/`, alongside `attribution/` (block→gene
from the GWAS) and `dissection/` (per-locus validation from the GEA). Those two are
organised by *which scan they came from*; this one exists precisely to **cross them**.

## Why

Neither upstream scan is calibrated well enough to rank on p-value:

| scan | inflation | reference |
|---|---|---|
| GEA — raw per-record LFMM, 22 climate axes | λ ≈ 1.72 median, up to **3.13** (pc3) | memory `raw-lfmm-over-wza-decision` |
| GWAS — 30 per-garden GEMMA scans on `s` | λ ≈ 1.0 but the **tail is 5–60× inflated**, worst where the hits are | memory `persite-gwas-low-mac-tail-inflation` |

Both failure modes are structure artifacts, but they act on **different traits**
(site climate vs. founder selection coefficient), **different units** (352 pools vs.
231 founders) and **different software** (LFMM vs. GEMMA). So a locus seen by both is
not explained by either one's inflation, and recurrence *within* a scan — across
uncorrelated climate clusters, or across independent gardens — is not what a single bad
axis or a single bad garden produces.

**Rank by independent lines of evidence, not by p.** This is the same triage rule that
emerged from `dissection/` (recurrence + frequency + LD-to-gene, not p-value), applied
across both tracks instead of within one.

## What goes into the pool

Set with the user 2026-08-26:

- **GEA:** every Bonferroni-significant record — *not* one lead per block. A block's
  lead can be intergenic while a marginally weaker record in the same block sits in a
  CDS; leads-only silently drops the genic one, which is the opposite of what the CARK
  lesson asks for. 22 axes × {sv, smallindel, nonsnp}, MAF > 0.05, per-(class × axis)
  genome-wide Bonferroni.
- **GWAS:** every Bonferroni-significant marker at the **MAC ≥ 5 (MAF 2.16%) arm**.
  The MAF-5% arm in `gemma_gwas_mac12/` is deliberately not used — raising the floor
  buys no calibration (the retained strata are still inflated 10.2× / 5.6× / 18.9× at
  p<1e-6) and only removes the worst stratum. MAC is carried per row instead.
- **Non-SNP hits only**, SNP-shadowed or not. Whether a Bonferroni SNP sits nearby is
  an annotation, never a filter.

Every variant is assigned to a gene **at its own position** via the TAIR10 GFF, reusing
`dissection/screen_sig_blocks.py`'s classifier rather than reimplementing it.

## The two corrections that change the answer

**1. Climate axes are correlated, so `n_axes` overcounts.** A variant significant on
five temperature axes has one line of evidence, not five. `axis_clusters.py` derives
the correlation structure from the data — average-linkage on 1 − |r| over the 31 site
climate vectors — giving **7 clusters at |r| ≥ 0.6**:

| cluster | axes |
|---|---|
| bio1 | bio1, bio5, bio9, bio10, bio11 |
| bio7 | bio2, bio4, bio6, bio7, pc2 |
| bio16 | bio12, bio13, bio16, bio19 |
| bio17 | bio14, bio15, bio17, bio18, pc1 |
| bio3 | bio3 |
| bio8 | bio8 |
| pc3 | pc3 |

Recurrence is counted over clusters (`n_clusters`), not axes. Note **pc3 — the most
inflated axis — forms its own cluster**, so it cannot be propped up by a correlated
neighbour.

**2. E1 and E2 are not independent.** Among genes recurrent across ≥2 clusters, 48%
have support from an axis with λ < 2, against 10% of the rest — recurring across
uncorrelated clusters is most of what it takes to land on a well-behaved axis. Summing
the four raw flags scores a purely GEA-internal result as "2 lines of evidence". So the
two are combined into **one gated GEA line**, and independence is counted over the three
genuinely separate sources:

| line | means |
|---|---|
| `L_gea` | recurrent over ≥2 climate clusters **and** ≥1 axis with λ < 2 |
| `L_gwas` | Bonferroni in ≥2 of the 30 gardens |
| `L_cross` | the gene is hit by **both** the GEA and the GWAS |

`n_evidence` (the raw 4-flag sum) is kept in the output for auditing. **Do not rank on
it.**

## A gene list is not a locus list

Neighbouring genes in one LD block are hit by the same haplotype. `add_loci()` collapses
candidate genes within **200 kb** into one locus and reports `locus_n_genes`. The Chr4
6.98–7.21 Mb region alone contributes **8** of the shortlist's genes — it is the cluster
already flagged as "possibly one adaptive haplotype" in memory
`genes-expl-candidate-pipeline`. This is the LD-scale version of
`sv-block-merge-artifact`. **Quote loci, not gene counts.**

## Scripts

| script | what it does | writes |
|---|---|---|
| `axis_clusters.py` | empirical correlation clusters of the 22 climate axes | `results/axis_clusters.csv` |
| `build_gea_pool.py` | all Bonferroni non-SNP GEA records, 22 axes → unique variants with cluster recurrence + λ | `results/gea_pool{,_records}.csv` |
| `build_gwas_pool.py` | all Bonferroni non-SNP per-garden GEMMA markers → unique markers with garden recurrence + MAC | `results/gwas_pool{,_records}.csv` |
| `build_convergence.py` | classify both pools at variant position, join at gene level, count independent lines, annotate | `results/{variants_classified,candidate_genes,convergent_genes}.csv` |
| `build_functional_variants.py` | SV/indels landing inside a gene or just before/after it, tiered by mechanism | `results/functional_{variants,shortlist}.csv` |
| `build_snp_tagging.py` | SNP co-significance + founder-panel tagging r², added in place | (columns on the above) |
| `build_master_table.py` | **the master table** — one row per gene, every column, annotated | `results/master_candidate_genes.{csv,md}` |
| `theme_master.py` | stress / temperature / circadian-light / flowering tags | `results/themed_candidates.csv` |
| `plot_candidate_overview.py` | shortlist + theme composition + genome distribution | `results/plots/candidate_overview.{png,pdf}` |
| `plot_top_loci.py` | drives `dissection/plot_locus_combined.py` for the top candidates | `results/plots/loci/` |

Run order: `axis_clusters` → `build_gea_pool` → `build_gwas_pool` → `build_convergence`
→ `build_functional_variants` → `build_snp_tagging` → `build_master_table` →
`theme_master` → the two plot scripts.

### The master table

`master_candidate_genes.csv` — **1,158 genes × 64 columns**, one row per gene, TAIR GO +
UniProt annotated for all of them (1,117/1,158 mapped). Column groups: identity ·
locus · discovery (`found_by`, axes, clusters, gardens, class scans) · blocks ·
variants · mechanism · SNP visibility · evidence · annotation · caveats.
`master_candidate_genes.md` is the 77-gene shortlist in readable form.

Two things it surfaces that the per-scan tables could not:

- **`found_by`**: 1,110 GEA-only, 42 GWAS-only, **6 by both**.
- **`block_gene_mismatch`**: **252 of 1,158 genes (22%)** sit in a clq0.9 block whose
  attributed gene set does not include them — block-level attribution would have named
  a different gene. That is the CARK artifact measured across the whole candidate set,
  and it is why everything here is assigned at the variant's own position.

### Themes

`theme_master.py` imports the regexes from `dissection/theme_filter.py` (one definition
of "flowering gene") but matches them against the **richer** text the master table
carries — UniProt free-text FUNCTION, UniProt keywords and GO biological process, not
just symbol and protein name. That is worth a lot: of 312 themed genes, **232 matched
only on free text** and would have been missed by name matching alone.

| theme | all candidates | in shortlist |
|---|---|---|
| stress | 200 | 16 |
| temperature | 96 | 8 |
| circadian / light | 47 | 3 |
| flowering | 63 | 2 |

Gating on evidence enriches for stress (17% → 21%) and temperature (8% → 10%), leaves
circadian flat (4% → 4%) and *depletes* flowering (5% → 3%) — so the shortlist is
tilted toward stress/temperature biology, but none of these shifts is large.

### Verdicts

`results/TOP_CANDIDATE_VERDICTS.md` — the cross-scan genes after locus dissection.
Short version: cross-scan agreement identifies a **locus**, not a gene, and the founder
LD-confirm decides which gene. CRK18 survives (r² 0.874 to its gene), GSH1 is moderate
(0.452), EPFL5 is suspect (0.391 vs 0.700 better elsewhere), **BT4 fails** (0.069 —
CARK-type detachment).

### Functional tiering (`build_functional_variants.py`)

"Which of these could actually *do* something?" — restricted to `size > 0` (SVs and
small indels; MNPs are substitutions, not length changes). Two things it adds over
`variants_classified.csv`:

- **A downstream flank tier.** The inherited classifier has a strand-aware 1 kb
  *upstream* promoter tier but no downstream one, so variants just past a gene's 3' end
  were falling into generic `7_proximal_intergenic`. Adding the mirror category recovers
  **272 variants** that had no gene handle before.
- **Frameshift status for coding indels** — `size % 3 != 0` inside CDS, the strongest
  functional prior available without expression data. SVs overlapping CDS are flagged
  separately (`cds_sv`), since a large deletion is likelier to remove exons than to
  shift a frame. Caveat: the modulo test assumes the indel lies wholly within coding
  sequence; one straddling an exon boundary is not handled.

Tiers, strongest mechanism first: `F1_CDS_frameshift` > `F2_CDS_inframe` > `F3_UTR` >
`F4_promoter` > `F5_downstream` > `F6_intron`.

Of 2,468 length-changing pooled variants, 1,715 have a gene handle (753 are TE or
intergenic with none): 24 CDS frameshift indels, 32 SVs overlapping CDS, 163 UTR, 742
promoter, 272 downstream, 445 intron. Restricted to genes carrying ≥1 gated line of
evidence: **160 variants in 72 genes / 39 loci**.

⚠ `locus_n_genes` in `functional_variants.csv` counts genes in the **pool**, not in the
shortlist, so it reads much higher there (57 for the Chr4 ~7 Mb locus). It is a pileup
warning, not a shortlist size.

### Are these visible to SNPs? (`build_snp_tagging.py`)

The pool keeps SNP-shadowed and SNP-unique candidates alike; these columns let you split
them afterwards. **They are annotations, never filters.** Two questions that are easy to
conflate and give opposite answers here:

| column | question | answer over the GEA-derived pool |
|---|---|---|
| `snp_cosig_2kb` | is the *association* also seen in SNPs? (a Bonferroni SNP within 2 kb, any of 22 axes) | **59% shadowed / 41% SNP-unique** |
| `best_r2_snp` | could a SNP *tag* this variant at all? (max founder-panel r² to any SNP within ±50 kb) | **median r² = 1.000**; only 3 of 1,421 testable are SNP-blind (r² < 0.2) |

So the non-SNP layer here is **not** finding variation SNPs cannot see — nearly every
candidate is in near-perfect LD with a nearby SNP. That is consistent with result #2 in
the tree README ("no kMate gain") and with `sv_indel_tagging_masked.ipynb`. The value of
this list is **causal candidacy**: an SV sitting in a promoter is a better mechanistic
hypothesis than the SNP tagging it, even though both mark the same haplotype.

r² is taken from the **masked** run (`sv_snp_ld_v2/tagging_{sv,indel}_panel_*`), not
`varexp/nonsnp_tagging_*` — the earlier build coded missingness as REF and had no MAC
floor, and panel SV records average ~60% missing. `tag_untestable` (no r² computable)
is reported separately from "untagged": they mean opposite things, and 242 of 1,663
GEA-derived functional variants are untestable, so absence of an r² is not evidence of
SNP-blindness.

### ⛔ `L_cargo` is withdrawn — a dinucleotide shuffle is not a null for long sequence

`cargo_line.py` turned the `cargo/` tree's TF-motif payload into an eighth line: what an
*insertion* brings with it, scored against each locus's own dinucleotide-shuffled null.
**It scores 0 everywhere and should not be revived by lowering a threshold.** Three steps
got there, and the third is the one that matters.

1. **The statistic counted motif hits, which are not independent.** One element is matched
   by every motif in its family at every offset — at Chr5:19,636,028 a single 44 bp ABRE
   draws 59 bZIP hits from 23 motifs — while a dinucleotide shuffle destroys clustered
   elements and rarely piles hits that way. One real element becomes dozens of counts
   against a null that cannot produce them.
2. **`cargo_sites_null.py`** redid all 193 loci with overlapping hits merged into distinct
   sites, everything else identical. 25/193 reached p ≤ 0.05 family-merged, median
   enrichment 1.06 — which *looks* calibrated. It is not: the rate is a function of
   **length**, 4.9% below 200 bp (nominal) against **63% above 2 kb**, at enrichments of
   only 1.14–1.19× with p pinned to the 1/101 floor. Concordance with the hit-count
   statistic is poor (13 of 22 shared), which is not what a sharpened real signal looks
   like.
3. **`cargo_null_control.py`** ran the identical statistic on **1,158 signal-free
   length-matched windows** — 3 random genomic and 3 TE-overlapping per real locus,
   everything downstream of sequence choice imported unchanged. Any p ≤ 0.05 there is a
   false positive by construction, so the per-length-bin rate *is* the null's calibration
   curve:

| statistic | <200 bp | 200–500 | 500 bp–2 kb | >2 kb |
|---|---|---|---|---|
| `hits_p` (original) | 11.1% | 13.6% | 30.8% | **64.9%** |
| `fam_sites_p` (family-merged) | 9.5% | 15.4% | 25.6% | **70.2%** |
| `any_sites_p` (all-motif-merged) | 2.1% | 3.7% | 3.4% | **3.5%** |

A dinucleotide shuffle preserves mono- and dinucleotide composition and nothing else, so
its bias is per-base and roughly constant while its sampling noise shrinks as sequence
grows — a fixed ~15% excess is unbeatable at 5 kb and invisible at 150 bp. **Only the
all-motif merge is calibrated at every length** (control FPR 2.8% overall, if anything
conservative): collapsing *any* overlapping hits is what reduces the unit to "distinct
pieces of DNA that match something", and a shuffle is a fair null for that count. Merging
within family is not enough, because a single GC- or AT-rich stretch recruits several
families at once. Note the TE arm is *less* inflated than the random-genomic arm, so this
is generic genomic sequence structure, not TE-ness.

**Head-to-head, there is no signal.** Real insertions sit at or *below* the false-positive
rate of signal-free sequence in **every** length bin of **every** statistic — pooled
`fam_sites_p` real 25/193 = 13.0% against a 16.5% control rate, and no bin reaches a
one-sided Fisher p < 0.43. On the one calibrated statistic the real loci give **5/193 =
2.6% against a 2.8% control rate**. The cargo README's "24 of 193 at p ≤ 0.05 vs 9.7
expected" excess, the 25 family-merged survivors, and the six-locus shortlist built on them
are all null mis-specification.

> The descriptive `cargo_*` columns stay — the sequence facts are still facts, and the
> merged-site columns are carried for auditing. What would revive this line is **a different
> control, not a different threshold**: the honest null for non-reference insertion sequence
> is the other panel SV insertions that are *not* GEA/GWAS hits — signal-free sequence of
> the same class — rather than reference windows. Until then no candidate gains or loses a
> line from its cargo, and `n_lines` tops out at 5.

### Which gene is a variant *about*? (`host_gene`, `prom_genes`)

Two columns look interchangeable and are not, and conflating them misnames genes.

| column | what it is |
|---|---|
| `nearest_gene` / `dist_to_gene` | the nearest GFF **`gene`** feature, by distance, ignoring strand |
| `gene` / `region` / `tier` | the gene the variant's **region call** is about — strand-aware |

`screen_sig_blocks.load_gff()` keeps only `typ == "gene"`, so
`transposable_element_gene` (3,903 features) and `pseudogene` are **absent from
`nearest_gene` by construction**. It therefore cannot name the TE gene a variant sits
inside, and it is never the right partner for a region call. Audited over the 193 cargo
insertion loci: of the 110 whose call claims a host gene, **36 (33%) would be given a
different gene** by `nearest_gene` — all 18 TE-gene calls and all 3 pseudogene calls (wrong
by construction, the named gene up to **193 kb** away), plus **15 of 57** promoters. The
worked case is the one the cargo dossier displays: Chr5:19,636,028 reads "promoter of
AT5G48460, 73 bp away", but AT5G48460 is on the minus strand and the insertion sits at its
3′ end — a terminator. The gene whose promoter it occupies is **AT5G48450 (SKS3)**, 364 bp
away.

So `build_sv_hit_dossier.py` publishes **`host_gene`** — the gene the call is about — with
`host_gene_basis` (`body` 31 / `te_pseudogene_body` 22 / `promoter` 57 / `none` 83) and
`host_dist_bp`, resolvable for 110/110. Body overlap wins over promoter adjacency because
one locus is both (Chr3_11718952 is inside TE gene AT3G29798 *and* in AT3G29800's
promoter). **An empty `host_gene` is meaningful**: the call claims no host at all (TE region
/ intergenic / gene desert, 83 loci), and `nearest_gene` there is a neighbour to be
reported as one — not as "*region* of *gene*". `cargo_line.py`'s independently derived
`target_gene` agrees with `host_gene` on 105/105 loci where both are populated.

**Promoters are not one-to-one.** Divergent (head-to-head) gene pairs share upstream DNA,
so a variant can sit in two promoters at once: **126 of 909** promoter variants in the pool
(13.9%) are in 2–4 promoters. `classify()` kept the first match *by coordinate* — where the
array happens to start, not a biological choice — and said nothing about the alternatives.
It now also emits `prom_genes` (all of them, nearest-TSS first), `n_prom_genes` and
`prom_gene_nearest_tss`. `gene` is unchanged, so nothing downstream moves; the ambiguity is
now visible instead of silently resolved. For reference, a nearest-TSS tie-break would
change the primary for **58 of 909** promoter variants across 43 gene pairs — including
shortlist genes (AT3G44010→AT3G44020, AT5G45740→AT5G45745, AT5G57270→AT5G57280,
AT1G13440→AT1G13448, AT5G48280→AT5G48290) — so it is **not** applied without a decision.

### ⚠ `size_inferred` — the multiallelic trap, in this pipeline

The per-garden GWAS arrays carry **only chrom/pos, no allele**, so `build_gwas_pool.py`
recovers the variant's length from the panel by taking the **largest record at the
position**. Where the position carries several records that is a *guess*, and
`size_inferred` flags it. For those rows `vclass`, `ftier`, the frameshift call, and any
join keyed on size are all unreliable — 21 rows in `functional_variants.csv`, 1 of them
carrying a frameshift call.

Worked example, because it hits the top candidate: the CRK18 frameshift indel at
Chr4:12,169,155 is **GWAS-derived**, and the panel has **8 records** at that position
with length differences 0,1,1,1,1,1,1,2. The pipeline picked the largest (2 bp). Seven
of the eight would still shift the reading frame, so "frameshift" survives — but "2 bp"
does not, and the tagging join on chrom+pos+size missed for the same reason. GEA-derived
variants are unaffected: `wza_in_clq09_tile` carries real ref_len/alt_len.

Run in the `kmate` env on a compute node, in that order. `build_gea_pool.py` reads
~2 GB of CSV (~35 s); `build_convergence.py` parses the full TAIR10 GFF and calls the
TAIR GO + UniProt annotator (needs outbound HTTPS).

## Result as built

7,283 GEA records → 2,938 unique variants; 150 GWAS hits → 82 unique markers. Union
3,019 variants, of which 1,762 land in a gene or promoter → **983 genes**.

| | genes |
|---|---|
| `L_gea` (gated GEA recurrence) | 59 |
| `L_gwas` (≥2 gardens) | 15 |
| `L_cross` (both scans) | 6 |
| **≥1 gated line** | **77 genes / 47 loci** |
| **≥2 independent lines** | **3** — CRK18, EPFL5, BT4 |

Sanity checks that the pool reproduces known work: **GPX6** (1,164 bp promoter deletion,
Chr4:7,011,705) comes out as the top GEA variant at 5/7 clusters, and **EMB1241** also
recovers — both were independently graded ROBUST in
`dissection/results/loci/CANDIDATE_VERDICTS.md`. The `L_cross` gene **CRK18**
(AT4G23260) sits at Chr4 ~12.1 Mb, the same region as the only block that survived
honest empirical-null recalibration in the retired WZA track (`wza-block-calibration`).

## Caveats that must travel with any candidate

- **Raw uncalibrated p on both sides.** This is a candidate net, not a calibrated hit
  list. Neither scan's multiple-testing correction accounts for the other, and the GEA
  Bonferroni is per class × axis, not across the 22 axes.
- **The GWAS cross-scan merge is keyed on chrom:pos**, because the per-garden result
  arrays carry no allele. 2.14% of arch3 positions are multiallelic
  (`panel-multiallelic-pos-key-trap`); `gwas_pos_multiallelic` flags every affected row.
  33/82 GWAS markers sit at such positions — in line with the 29.3% baseline for tested
  non-SNP markers, so this is not an enrichment, but the allele is ambiguous there.
- **61 of 82 GWAS markers are MAC < 12**, the worst-calibrated stratum.
- Nothing here has been through `dissection/` (founder-LD to the gene, haplotype panel).
  A shortlist entry is a hypothesis until it has.

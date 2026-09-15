# Result 1 — negative selection on SVs

**The question.** Comparing SVs and SNPs *at the same initial frequency*, is
there more negative selection on SVs, and does it track climate?

Everything here is per-variant and temporal: no haploblock unit, no GWAS. The
estimand is a selection coefficient `s` from the allele-frequency trajectory
across generations, with replicate plots within a site as the replication.

Outputs → `results/` (`sv_adaptive/`, `site_temporal/`); figures → `results/*/plots/`.

**Read the result in:** `../notebooks/temporal_s_consolidated.ipynb` — the
consolidated notebook (2026-07-15) that replaced seven separately-named ones.
Supporting: `sv_parallelism_climate.ipynb`, `sfs_shift_by_site.ipynb`,
`temporal_s_nofilter.ipynb`.

> **Streamlined 2026-08-25** to headline + direct support. Archived, each reduced to a
> one-liner in its archive README: bulk parallelism + PicMin
> (`archive/replicate_arm_bulk_2026-08-25/`) and `sfs_time_site4`
> (`archive/sfs_time_site4_2026-08-25/`).

Narrative + caveats: `SV_TEMPORAL_PURGING_SUMMARY.md` (this folder).

> ⚠ **Do not retire the `_temporal_s_*` / `_compute_s_*` family on name
> similarity.** They look like successive versions but were audited 2026-07-15
> and each tests a genuinely different statistic. The docs around them contain
> *retractions* of a superseded label, not declarations of one.

---

## A. The `s` estimator and the SNP-vs-SV comparison

| script | what it does | writes |
|---|---|---|
| `_temporal_selection_snp_vs_nonsnp.py` | the origin: pure per-variant, per-site temporal test | `temporal_selection_snp_vs_nonsnp_summary.csv` |
| `_temporal_s_snp_vs_nonsnp.py` | switches Δp → selection coefficient `s`, multi-generation | `temporal_s_snp_vs_nonsnp.csv` |
| `_temporal_s_plots_snp_vs_nonsnp.py` | `s` with **plots as replicates** — the power-gaining refinement | `temporal_s_plots_snp_vs_nonsnp.csv` |
| `_compute_s_dist_by_stratum.py` | **the headline computation**: distribution of `s` by class *within initial-frequency strata*, per site | `s_dist_by_stratum.csv`, `..._sitemeta.csv` |
| `_temporal_s_enrich_initqty.py` | SNP-null-free variant: matches on initial quantity instead of a SNP null | `temporal_s_enrich_initqty.csv` |
| `_nonsnp_temporal_category.py` | asks the same of non-SNP as a whole category, not SVs alone | — |
| `_temporal_sel_drift_maf.py` | robustness: explicit drift and MAF controls | `temporal_sel_drift_maf.csv` |
| `_build_temporal_s_consolidated_nb.py` | builds the consolidated notebook + its figures | `temporal_s_consolidated.ipynb` |

## B. Climate

| script | what it does | writes |
|---|---|---|
| `_compute_s_climate_slope.py` | the two climate-differential views, constructed to cancel the founder confound | `s_climate_slope*.npz`, `s_climate_slope_sign_by_site.csv` |

## C. Artifact controls — what makes section A defensible

Each answers a specific "could this be an artifact?" challenge. Keep them with
the result; they are the reason it survives.

| script | the challenge |
|---|---|
| `_audit_s_classes.py` | are the SNP/indel/SV `s`-distributions matching because of a bug, or really? |
| `noise_check_sv_s.py` | is the SV > SNP \|s\| excess a k-mer-support (noise) artifact? |
| `_founder_load_test.py` | founder-level test of the "selection on SVs" conclusion (audit 2026-07-06) |
| `_sv_temporal_direct.py` | do common SVs actually fall in frequency over generations? |
| `_sv_callqual_artifact.py` | calling-confidence proxy for the insertion-polarity artifact-vs-biology question |
| `_extract_vcf_callqual.sh` | pulls the `F_MISSING` / `CONFLICT` / `MA` fields the above needs |

## C2. Founder-level confound analysis — is the climate gradient real?

Added 2026-08-26. In global mode per-variant AF is a deterministic projection of the
founder mixture `h` through fixed carrier sets, so the per-variant climate slope carries
information only through *which founders carry the variant* and *how those founders
respond to climate*. These scripts test the founder-level claim that follows, and its
confounds. **Read the result in `../notebooks/sv_founder_confound.ipynb`.**

**Headline measure: `kb_ins` — kb of sequence present in a founder and absent from Col-0**
(median 121 kb, range 44–202 kb). Replaced `ins_frac` (count / total variant load) on
2026-08-26: raw bp is both the plainest to state and the best behaved on every test.

| script | what it establishes | writes |
|---|---|---|
| `_founder_load_test.py` | the founder-level fact (audit 2026-07-06): **insertion-rich founders decline in hot/arid gardens**. Its T2/T3 carrier-composition sections are a wiring sanity check, **not evidence** — see the note below | `founder_load_test.npz` |
| `_founder_climate_confound.py` | panel (cactus vs PG) and provenance controls; raw + partial correlations | `founder_climate_confound.{npz,csv}` |
| `_founder_sv_content.py` | per-founder SV content every way — count and **bp**, ins/del/both, raw and normalised — plus SNP-ALT divergence from Col-0 and **per-founder call rate** | `founder_sv_content.npz` |
| `_founder_sv_content_climate.py` | which measure tracks climate, and how entangled each is with Col-0 divergence | `founder_sv_content_climate.csv` |
| `_founder_sv_technical.py` | panel half on **bp** (never checked before), call-quality control, within-panel-half slopes | `founder_sv_technical.csv` |
| `_founder_col0_climate_distance.py` | climatic distance to the Col-0 backbone: direction vs distance, with a European-origin sensitivity anchor | `founder_col0_climate_distance.csv` |
| `_founder_graph_representation.py` | unfiltered divergence, kinship to the assembly set, within-PanGenie test, **assembly size** from the `.fai` indexes | `founder_graph_representation.{npz,csv}` |
| `_founder_assembly_technology.py` | **sequencing platform**, assembler, N50, gaps — from the assembly release metadata | `founder_assembly_technology.{npz,csv}` |
| `_build_sv_founder_confound_nb.py` | builds the notebook + its nine figures | `sv_founder_confound.ipynb` |

**Verdict — four confounds tested, none explains the effect:**
- **Panel asymmetry: ruled out for kb, REAL for counts.** kb vs cactus/PG ρ=−0.056 (p=0.40),
  but insertion **count** ρ=+0.378 (p<0.001) — long-read assemblies resolve more insertion
  *records*. Normalising by variant load does not fix it, it flips it (bp/load −0.224,
  p=0.001). This is what selects the measure.
- **Provenance: real, partial, not sufficient.** Origin bio1 is the strongest single predictor
  of the founder climate response (ρ=+0.632), and cold-origin ecotypes really do carry more
  inserted sequence (ρ=−0.542); but the attenuation **plateaus** — 2 → 19 origin bioclim axes
  barely moves it, and the chain's third link only falls −0.569 → −0.348 with origin fixed.
- **Col-0 reference bias: ruled out for kb.** Divergence from Col-0 is not climate-structured
  (vs origin bio1 ρ=−0.063, p=0.34). SV **counts** are near-pure divergence proxies (deletion
  count vs SNP divergence ρ=**+0.958**, SV count +0.947); **kb of insertions is not**
  (ρ=+0.088, p=0.18). Climatic distance to Col-0 gives a double dissociation: insertions are
  purely directional, deletions purely distance — and the deletion-distance link vanishes
  entirely once genetic divergence is held fixed (+0.312 → −0.071, p=0.29).
- **Panel construction: ruled out.** Call rate differs by half (cactus 0.9967 vs PG 0.9608,
  p=3e-33) and tracks origin climate (ρ=−0.249), but fixing it moves kb only −0.542 → −0.506.
  Decisively, the relationship holds **within the PanGenie half, which has no assemblies at
  all** (ρ=−0.576, p=1e-14, vs cactus −0.523) — an assembly-quality or long-read-technology
  effect cannot do that.

Deletions show nothing on any measure and serve as the internal null control. What survives is
founder-level, **not** evidence of per-SV selection.

- **Pangenome construction: ruled out.** Better instruments than the MAC-filtered `n_snp`:
  unfiltered divergence (`n_snp_all`, 524k vs 284k) behaves identically (vs origin bio1
  ρ=−0.073, p=0.27); **graph representation** (kinship to the 80 assembly founders) is
  unrelated to origin climate (ρ=−0.009, p=0.89) and controlling for it *strengthens* the
  result (−0.542 → −0.567); **assembly size** is flat against everything (vs origin bio1
  ρ=+0.069 p=0.54; vs inserted kb ρ=−0.043 p=0.70), so cold-origin accessions do not simply
  have bigger genomes. Within the PanGenie half alone, kb~origin goes −0.576 → **−0.633**
  once graph representation and unfiltered divergence are held fixed.
- **Sequencing technology: real effect, but orthogonal.** Metadata found at
  `pang/long_read_seq_ara/ASSEMBLIES_Best_version_of_dataset.csv` (`Primary_Sequencing_Technology`,
  `Assembler`, `N50_contigs`, `Gaps_Scaffolds`); 80/80 cactus founders matched. Platform **does**
  shift measured inserted kb (Kruskal-Wallis p=0.0076 — ONT 85 kb / ONT_R10.4 99 kb vs HiFi 125 kb
  / CLR 123 kb, a ~30% deficit for ONT-based assemblies) — a genuine technical effect not previously
  quantified. But platform is **not** confounded with ecotype origin (H=1.81, **p=0.77**), the
  relationship holds within every platform (CLR −0.52, HiFi −0.48, ONT −0.65), and technology +
  N50 + assembly size held together move kb~origin only −0.523 → **−0.518**. Contiguity runs the
  wrong way for an artifact: higher N50 gives *less* inserted sequence (ρ=−0.235, p=0.036).
>
> ⚠ **Col-0's climate record is wrong for this purpose.** The 1001G table places ecotype 6909
> in Missouri, USA (38.3, −92.3; bio1 13.1 °C) — where the Laibach/Rédei lineage was named, not
> collected. No absolute "climatic distance to Col-0" number is biologically meaningful; a
> European-origin anchor gives the same qualitative answer. Col-0 is also **not** among the 231
> founders or the 80 cactus assemblies, so there is no zero-divergence anchor.

> ⚠ `founder_load_test.npz` was stale until 2026-08-26 (built 2026-07-06 against inputs
> regenerated 2026-07-15 and 2026-08-25). Regenerate it whenever `s_climate_slope.npz` or
> `selection_s_matrix.npz` move.
>
> ⚠ **Do not cite T2/T3 of `_founder_load_test.py` as a result.** "Carrier composition predicts
> per-variant β" restates kMate's projection identity (`AF_v = h·V_pa/h·V_called`) and cannot
> fail; its ρ=0.67 measures only how lossy the ḡ scalar is, not a share of variance. The true
> share is 100% by construction. Keep it as a wiring check.
>
> ⚠ A K_snp mixed model was tried and its p-values are **not usable** — the REML variance
> component pins at the boundary (logL −7.91 at δ=1e-6 vs −30.7 at δ=1) and K_snp's condition
> number is 1.5e15.
>
> **Cut 2026-08-26 (user):** the structural/carrier-composition section, the genome-wide
> divergence control (superseded by the Col-0 analysis), and the population-structure controls
> (clade decomposition, within-clade permutation, jackknife). Those tests all passed; they were
> dropped for scope, not because they failed. `_founder_confound_structure.py` and
> `_founder_provenance_depth.py` still hold the structure and deep-provenance numbers if needed.

## D. Parallelism / PicMin arm

Does the same signal appear as parallel change across replicate plots and sites?

| script | what it does | writes |
|---|---|---|
| `_compute_parallelism.py` | per-variant parallelism across replicate plots (AF-vapeR / PicMin philosophy) | `parallelism.npz` |
| ~~`_picmin.py`~~ | PicMin (Booker et al.) — **archived 2026-08-25**: ~1.07-1.12x near parity, and direction-agnostic (`|z|` input) so it could not corroborate a purging result regardless. `archive/replicate_arm_bulk_2026-08-25/` | — |
| `_site_parallelism.py` | per-site parallelism of *founder* frequency change | `site_parallelism.csv` |
| `_plot_site_parallelism.py` | per-site purging vs per-site parallelism — the confound scatter | `site_parallelism_vs_purging.png` |
| `_build_sv_parallelism_climate_nb.py` | builds the surviving replicate-arm section (per-site parallelism excess vs climate) | `sv_parallelism_climate.ipynb` |

## E. Site-frequency-spectrum arm

A different signature from the central-tendency test: the SFS *mean* shift is
null, but SV **extinction rate** is elevated.

| script | what it does | writes |
|---|---|---|
| `_compute_sfs_shift_by_site.py` | folded-SFS shift per site, per class, per p0 decile | `sfs_shift_by_site.csv`, `..._sitemeta.csv` |
| `_build_sfs_shift_nb.py` | builds the SFS-shift notebook | `sfs_shift_by_site.ipynb` |
| `_compute_sfs_time_site4.py` | AF spectrum over generations 0→3 at site 4, by class | `sfs_time_site4.npz`, `..._summary.csv` |
| `_build_sfs_time_site4_nb.py` | builds the site-4 SFS-over-time notebook | `sfs_time_site4.ipynb` |

## F. Per-site block enrichment (site-4 pilot → cross-site)

The one block-unit arm in r1: are temporally-selected blocks SV-enriched?

| script | what it does | writes |
|---|---|---|
| `site_variant_temporal_scoef.py` | per-variant `s` at one site, by class | `site{N}_scoef_*.npz` |
| `site_sv_enrichment.py` | SV enrichment in selected clq0.9 blocks at one site | `site{N}_sv_enrichment.json`, `site{N}_clq90_blocks.csv.gz` |
| `run_site_enrichment.sbatch` | SLURM driver for the above across sites | — |
| `aggregate_sites_enrichment.py` | cross-site synthesis | `cross_site_enrichment.csv/.png` |
| `_enrich_threshold_sweep.py` | is the enrichment an artifact of the top-fraction cut? size-matched permutation null | stdout only |
| `_render_site4_enrichment_fig.py` | renders the site-4 figure from precomputed plot data | `site4_sv_enrichment.png` |

---

## G. What is in the insertions? (site context + cargo)

Added 2026-08-26, downstream of C2: if cold-origin founders carry more inserted sequence
and decline as gardens warm, what *is* that sequence? **Read the result in
`../notebooks/sv_insertion_content.ipynb`.** All 172,220 SV insertions, **no MAC floor**
(every other section of this arm filters MAC>=12; that is wrong here, because purifying
selection lives in the rare tail that filter removes).

Three layers, because an insertion is by definition absent from TAIR10 and so has no
reference coordinates of its own:

| script | what it establishes | writes |
|---|---|---|
| `_insertion_genomic_context.py` | **where it landed** — CDS/UTR/intron/non-coding-exon/intergenic + TE-overlap flag, from the TAIR10 `genes_transposons` GFF, plus carrier-mean origin bio1/bio18 | `insertion_context.{csv,npz}` |
| `_extract_insertion_seqs.py` | the inserted sequences themselves, per chromosome | `seq/insertions_Chr{N}.fa` |
| `_liftback_to_assemblies.py` | locates each insertion in a **carrier's assembly**, where it is an exact substring | `liftback_assignment.csv` |
| `_liftback_recover_short.py` | recovers the short insertions the first pass length-biased away | — |
| `_annotate_liftback.py` | **what it carries** — intersects the lifted interval with the Helixer / Liftoff / TRASH / repeat-compartment tracks precomputed on that assembly | `liftback_annotation.csv` |
| `_insertion_vs_tair10.py` | **where else the sequence occurs** — dc-megablast vs TAIR10: duplication-vs-novel + **TE family** from the `Alias` attribute | `insertion_tair10_class.csv`, `blast/Chr{N}.tsv` |
| `_build_insertion_content_nb.py` | builds the notebook + its four figures | `sv_insertion_content.ipynb` |

**Result.** Insertions are **5.2x depleted in CDS** (log2 obs/exp -2.37) with a clean
purifying gradient through UTR (-0.65) -> intron (-0.28) -> exon_noncoding (-0.24) ->
intergenic (+0.65). 89.7% lift back. **52.9% of the cargo is TE-derived** across 314
families, led by the non-autonomous ATREP/Helitron group; 28.3% is duplicated Col-0
sequence; the cold-origin excess of section C2 is specifically a **TE excess** (cold > warm
in 5/5 frequency-matched strata). Organellar contamination screen clean (0.14%).

Three things the per-family and per-class breakdowns add (notebook sections 4-5):
- The TE excess is **broad-based, not family-driven** — 16/33 families with n>=300 are
  cold-skewed, a coin flip. It is a compositional shift across the TE landscape.
- **ATCOPIA78 is the one real outlier**: cold:warm 2.3x (log2 +1.22), the only family of 33
  to survive Bonferroni, and the **youngest cargo in the dataset** (median 99.4% identity to
  TAIR10 vs 94.6% for TE-derived overall, median 4,961 bp = near full length). ATCOPIA78 is
  ONSEN, the heat-activated retrotransposon. **Carry this as a lead, not a result** — the
  family name is TAIR10's `Alias`, nothing here tests activation, and the direction
  (cold-origin founders carrying more of a heat-responsive element) is not the naive
  prediction.
- **Cargo class barely predicts frequency**, against the sharp CDS/intergenic gradient of
  section 1. Purifying selection here is about *where the insertion landed*, not what it
  carried.

> **Two design choices worth not re-litigating.**
> **Lift back, do not annotate the fragments.** Helixer's minimum record length is 25 kbp
> and its land-plant window is 21-107 kbp; the median insertion is 754 bp. Running an ab
> initio gene finder on the fragments is out-of-domain, not merely less accurate.
> **dc-megablast, not minimap2.** minimap2's `asm5/asm10` presets use k=19,w=19 minimizers
> and chain scores tuned for assembly-scale contigs, so at this query length they silently
> miss short and diverged copies — the same preset trap that length-biased the first
> lift-back pass.

> ⚠ **`novel` measures alignment power, not biology.** A short query has little statistical
> power in a homology search, and the class share shows it directly: `novel` is 64.1% of
> insertions <100 bp and **1.0%** of those >10 kb. The pooled 18.7% figure is an artifact
> ceiling. Always read the cargo composition **within a size bin**.

> ⚠ **`carrier_bio1` is mathematically coupled to allele frequency** — a singleton's
> carrier-mean *is* that one founder's origin, a common insertion regresses to the panel
> mean (SD falls 4.83 -> 1.42 across carrier-count strata). A naive quartile split pulls
> singletons into both tails, so any frequency comparison between the groups is confounded
> *by construction* (unmatched: 56.5% vs 78.9% private). Section 2 splits cold/warm
> **within** carrier-count strata, and section 3/4 inherit that same matching. The unmatched
> numbers are printed for contrast only and must not be quoted.

> ⚠ The track named `02_annotation_RepeatMasker` on the assemblies is **not** a TE-family
> annotation despite the name — it holds only centromere, telomere, 45S/5S rDNA,
> chloroplast, mitochondria and N_stretch. It is useful as a free **organellar-contamination
> screen**. TE families come from `_insertion_vs_tair10.py` instead.

> **Resolved 2026-09-15: the ab initio TE-ORF worry does not bite, it runs the other way.**
> The concern was that Helixer reads retrotransposon *gag*/*pol* as coding, inflating the
> 28.0% de novo gene rate specifically in the TE fraction. Measured against the section-4 TE
> calls, TE-derived cargo has the **lowest** Helixer rate of any class (22.4%, vs 66.2% for
> `gene_dup`), and dropping it *raises* the overall rate to 34.3%. The separate caveat that
> "overlaps a Helixer gene" is not "carries a gene" still bounds the number from above.

> ⚠ **Still open.** Tandem repeats beyond TRASH (via ULTRA;
> TRF mis-annotates >30% on AT-rich genomes); the per-variant climate-slope beta split
> (exists only for MAC>=12 and is indexed by AF-store column order, so joining it needs an
> ordinal two-pointer walk, **not** a `chrom:pos` join — 2.14% of arch3 positions are
> multiallelic and a position join silently matches the wrong allele); and the long tail
> above ~20 kb (p99 = 19.8 kb, max 291 kb), more plausibly segmental duplication or
> mis-assembly than insertion.

---

## Filed here but arguably belonging elsewhere

Both are **builders whose notebooks are about a different question** — flagged,
not moved, pending a decision:

- `_build_sv_selection_audit_nb.py` → `sv_selection_haplotype_audit.ipynb`. This
  is the haplotype-unit SV-enrichment audit; its compute scripts live in
  `../extras/`. Only the builder stayed behind.
- `_build_sv_polarity_manhattan_nb.py` → `sv_polarity_enrichment_clq90.ipynb`.
  Content is a WZA/clq0.9 GEA-block enrichment (Kendall + LFMM K=16 + binomial
  vs bio1), i.e. r2 material.

## Also here

- `BAYPASS_TEMPORAL_PLAN.md` — a **plan, not completed**; its own banner marks
  its haploblock-frequency inputs as stale (the retired Pipeline-B chain).

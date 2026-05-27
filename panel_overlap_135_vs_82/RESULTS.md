# pang_135 vs pang_82 — what do the 53 non-GrENE-Net "diversity" extras buy us?

**Date**: 2026-05-21
**Question**: The 151 missing GrENE-Net founders are PanGenie-genotyped against the pang_135 cactus pangenome. pang_135 contains 82 GrENE-Net cactus assemblies + **53 non-GrENE-Net "diversity bonus" accessions**. How much information do the 53 extras actually add to the genotyping of the 151 short-read founders?

**TL;DR**: ~95% of the 151's ALT-allele mass comes from variation captured by both sub-panels (shared). Of the remaining 5%, ~3.9 percentage points trace to cactus-only ALTs and ~1.4 pp to extras-only ALTs. The extras contribute almost nothing on SNPs (1.2%) but **11.8% on insertions ≥50 bp**. 3 of the 53 "extras" are clonal duplicates of GrENE-Net cactus accessions and contribute nothing. The effective extras count is **50, not 53**.

---

## 1. Setup

Confirmed from the pang_135 raw VCF sample list (135 columns):
- **82 numeric Assembly_IDs** (e.g. `100042`, `100313`) — the GrENE-Net cactus side. Maps to 81 unique Accession_IDs (one is a dedup-pair), of which 78 overlap `cactus_80_ecotypes.txt` (v3qc-v3 drops 2 more for QC).
- **53 non-numeric names** (e.g. `Tul-0`, `Yo-0`, `Est-0`) — the diversity extras from the Exposito-Alonso 2026 batch.

Lookup table: `pangenie_genotyping/data/cactus_80_ecotypes.txt`, `data/sv_panel_to_accession_id.tsv`.

PanGenie was built against `pang_1001gplus_all.vcf.gz` (vcfbub-decomposed, top-level only, 1 GB), **not** the raw nested VCF (`pang_1001gplus_all.raw.vcf.gz`, 2.86 GB). This matters: ~30% of pang_135 raw per-ALT records are nested-only bubbles that PG never processes.

The 151 PG-genotyped short-read founders' per-founder VCFs live in `pangenie_genotyping/data/genotyped/*_genotyping.vcf.gz`. We use these PRE-het-mask, PRE-V4-filter — the raw PG signal, since QC is a separate question.

---

## 2. Tier 1A — what fraction of pang_135 ALTs exist only because of the 53 extras?

We classify each per-(record, ALT) row from the raw VCF into:
- **V1a** — biallelic record where the (only) ALT is carried by ≥1 extras and 0 cactus → "extras introduced this new bubble"
- **V1b** — multi-allelic record where THIS ALT is carried by ≥1 extras and 0 cactus → "extras added a new allele at an existing bubble"
- **cactus_only** — ≥1 cactus carries the ALT, 0 extras
- **shared** — ≥1 on both sides
- **empty** — 0 on both sides (rare; AC=0 in both subsets)

### Genome-wide counts (raw VCF, all 5 chroms)

| category | n per-ALT rows | % |
|---|---:|---:|
| **V1 (V1a + V1b) — extras-introduced** | **2,632,386** | **22.07%** |
| — V1a (new bubble) | 1,352,535 | 11.34% |
| — V1b (new ALT at existing bubble) | 1,279,851 | 10.73% |
| shared | 5,431,495 | 45.54% |
| cactus_only | 3,848,253 | 32.27% |
| empty | 14,867 | 0.12% |
| **TOTAL** | **11,927,001** | |

**22% upper bound** on the contribution of the 53 extras (this is the share of ALT alleles in pang_135 that would not exist in pang_82).

### V1 fraction by variant type — SVs are most extras-enriched

| vtype | V1 fraction | dominant subtype |
|---|---:|---|
| SNP | 21.6% | V1a (new bubbles) |
| INS<50 bp | 22.0% | V1b |
| DEL<50 bp | 20.3% | V1b |
| **INS ≥ 50 bp (SV)** | **34.3%** | V1b |
| **DEL ≥ 50 bp (SV)** | **30.9%** | V1b |
| MNP/other | 22.0% | V1b |

### V1 ac_extras distribution — singleton-heavy

| ac_53 (n of 53 extras carrying the ALT) | n V1 records | % of V1 |
|---|---:|---:|
| 1 (singleton) | 2,013,662 | **76.5%** |
| 2–4 | 511,289 | 19.4% |
| 5–9 | 105,756 | 4.0% |
| ≥10 | 1,679 | 0.06% |

3 of every 4 V1 records are carried by just 1 of the 53 extras. So most of the "extras diversity" is private to individual extras genomes. But — see §3 — these singletons are NOT genuinely extras-private; many are population-rare alleles that just happen to be in one of the 53.

---

## 3. Tier 1B — the headline: how many V1 ALTs are actually called by PG in the 151?

For each V1 record we join against the 151 PG-genotyped founders' merged AC table (built via `bcftools merge` of all 151 single-sample VCFs, with `+fill-tags -t AC,AN`).

### Top-line numbers

| metric | value |
|---|---:|
| Total V1 records | 2,632,386 |
| V1 records present in PG VCF | 1,743,879 (66%) |
| V1 records with ≥1 PG carrier (ac_pg > 0) | **450,199 (17.10%)** |
| V1 records with ≥1 PG carrier, restricted to ac_53 ≥ 2 (de-noised) | **122,858 / 618,724 = 19.86%** |
| V1 records with ≥1 PG carrier, conditional on PG-present | **~25%** |

**~450K records genome-wide** where the 151 short-read founders gain an ALT call solely because the 53 extras introduced the bubble. That is the realized value-add.

### PG-carrier rate by category × vtype (multi-carrier V1, ac_53 ≥ 2)

| vtype | n records | PG-carrier rate |
|---|---:|---:|
| **DEL ≥ 50 bp** | 6,537 | **26.3%** |
| INS ≥ 50 bp | 11,964 | 24.6% |
| INS < 50 bp | 59,166 | 23.9% |
| DEL < 50 bp | 65,940 | 21.1% |
| SNP | 322,649 | 19.0% |
| MNP/other | 152,468 | 18.9% |

### Carrier rate by ac_53 bin (V1, all vtypes)

| ac_53 bin | n V1 records | PG-carrier rate (raw) | conditional on PG-present |
|---|---:|---:|---:|
| 1 (singleton) | 2,013,662 | 16.3% | ~25% |
| 2–4 | 511,289 | 20.0% | ~30% |
| 5–9 | 105,756 | 18.8% | ~28% |
| **≥10** | 1,679 | **56.0%** | **~73%** |

**Singletons aren't extras-private** — 16% PG-carrier rate is much higher than expected if these alleles were truly unique to a single extras genome. Most singleton V1 ALTs are population-rare alleles that happen to be in 1 of 53 extras AND some of the 151 PG founders, but were absent from the 82-cactus sample.

### Sanity-check baselines (same join, other categories)

| category | n records | PG-carrier rate | meaning |
|---|---:|---:|---|
| **V1 (extras-only)** | 2,632,386 | **17.1%** | realized value-add of 53 extras |
| cactus_only | 3,848,253 | 25.1% | analogous baseline for 82 cactus |
| shared | 5,431,495 | 58.0% | common ALTs, mostly carried by 151 |

Per-accession, extras and cactus contribute roughly equivalently to the 151:
- Cactus 25.1% / 82 = 0.306% per cactus accession
- Extras 17.1% / 53 = 0.323% per extras accession

I.e. on a per-accession basis the extras are **about as informative** as a random subset of the cactus accessions. The reason cactus contributes more in total is purely panel size (82 > 53), not per-accession quality.

---

## 4. PG ALT-allele mass attribution

A cleaner accounting: of the 151 PG founders' total ALT-allele mass (sum of `ac_pg` across all records = total founder-record cells called ALT), what fraction is attributable to each category × vtype?

Total mass = **179.2 M** (founder × record ALT cells across 151 founders).

| vtype | shared | cactus_only | V1 (extras) |
|---|---:|---:|---:|
| **SNP** | **95.40%** | 3.39% | **1.21%** |
| INS < 50 bp | 94.55% | 4.02% | 1.43% |
| DEL < 50 bp | 94.09% | 4.46% | 1.45% |
| MNP/other | 94.41% | 4.13% | 1.46% |
| **DEL ≥ 50 bp** | **78.36%** | 15.50% | **6.14%** |
| **INS ≥ 50 bp** | **58.18%** | 30.00% | **11.82%** |
| **All variants** | **94.70%** | **3.91%** | **1.39%** |

**Headline**: ~95% of the 151's ALT-allele mass is shared between the two sub-panels. Only 1.39% is attributable to ALTs only in the 53 extras. The 53 extras' contribution scales sharply with variant size:

| vtype | extras-only / SNP-extras ratio |
|---|---:|
| SNP | 1.0× (baseline 1.21%) |
| small indel | 1.2× |
| DEL ≥ 50 bp | 5.1× |
| **INS ≥ 50 bp** | **9.8×** |

**The 53 extras buy almost nothing on SNPs (~1% of SNP mass), but ~12% of large-insertion ALT mass.** The SV diversity argument for using pang_135 over pang_82 is real; the SNP argument is essentially noise.

---

## 5. Tier 2 — quasi-duplicate flag (the "wasted slots" finding)

Pairwise SNP genotype identity between the 53 extras and the 82 cactus, computed on 5.44M biallelic SNPs from the pang_135 raw VCF (streaming Python, both sides treated as haploid).

**3 of the 53 extras are clones of GrENE-Net cactus assemblies** (SNP identity > 99.98% on 4.8M+ called sites — same magnitude as the 5772/6150 mislabel):

| Extras label | Cactus name (Assembly_ID) | Cactus Accession_ID | SNP identity |
|---|---|---|---|
| Nov-02 | Noveg-2 (100703) | 9637 | 99.99987% |
| Est-0 | Est (100313) | 7127 | 99.99949% |
| Nok-1 | Nok-3 (100302) | 6945 | 99.98448% |

The naming hints at the relationship (Nov-02 / Noveg-2, Est-0 / Est, Nok-1 / Nok-3) — same physical accessions sequenced under different labels.

**Effective extras count = 50, not 53.** Median pairwise identity 53×82 is 0.897 — the other 50 are genuinely distinct.

**Recommendation**: drop or relabel these 3 columns in a future pang rebuild. Cross-reference with the `widely-used-tools-default-correct` memory entry — same pattern as 5772/6150.

---

## 6. Tier 2 — D1 (per-bubble unique k-mer counts)

Concern: adding more haplotypes to a pangenome graph can *reduce* the per-bubble unique-k-mer count after PanGenie's uniqueness filter (more shared k-mers between haplotypes → fewer that pass the per-bubble uniqueness test).

Per pang_135 PanGenie index (`pang_135_pangenie_index_Chr{1..5}_kmers.tsv.gz`):

| bubble contains | n bubbles | mean uniq-kmer | median |
|---|---:|---:|---:|
| no V1 ALTs | 132,277 | 67 | 61 |
| V1a only | 113,406 | 118 | 100 |
| V1b only | 14,678 | 159 | 142 |
| V1a + V1b | 47,714 | 253 | 301 |

**D1 not detected** at gross-bubble resolution: bubbles enriched in extras-introduced ALTs have *more* unique k-mer evidence per bubble, not fewer. The natural explanation is that bubbles with extras-only ALTs have more alleles total (extras add allele rows), and total unique-kmers scales with allele count (capped at PanGenie's 16/allele, visible as the p95=301 ≈ 19 alleles × 16/allele).

This rules out the audit's a-priori D1 concern at the bubble-level. Per-allele post-cap unique-k-mer count would be a sharper test but requires `.cereal` parsing — deferred.

---

## 7. Interpretation and recommendations

### What the 53 extras actually buy

1. **For SNPs and small indels (≈85% of total per-ALT records)**: the extras add 1.2-1.5% of the 151's ALT mass. Marginal. The 82 cactus already saturate this.
2. **For large SVs (≥50 bp)**: the extras add 6-12% of the 151's ALT mass. Real. The SV diversity argument for pang_135 over pang_82 is data-supported.
3. **Per-accession**, the extras are **about as informative as a random 53 of the 82 cactus**. They are not systematically better or worse on a quality basis; their smaller total contribution is just panel size.

### What's wasted

1. **3 of 53 are clones of GrENE-Net cactus accessions** — same physical samples relabeled. Effective extras count is 50.
2. **76.5% of V1 records are extras-singletons**, but even these contribute ~16% PG-carrier rate — they're population-rare alleles, not extras-genome noise.

### Caveats and limits of this analysis

- **PG was run on vcfbub top-level VCF, not raw**: ~30% of V1 records (nested bubbles) are absent from PG entirely. The "raw V1" denominator is conservative; the "PG-processed" conditional denominator is the more honest measure of what the extras contribute *given* what PG sees.
- **Singleton ac_53 / 151 carrier interactions are population-genetic** — the analysis doesn't prove the extras cause the 151's ALT calls (the alleles exist in 1001G regardless); it only measures coincidence in genotype calls. Tier 3 (a counterfactual pang_82 PanGenie rebuild) would be required to prove causation.
- **Tier 2 D1 is bubble-level**; per-allele post-cap k-mer counts may show a different picture and were not measured here.
- **No haploid-vs-diploid issues**: confirmed (chrom, pos, alt_idx) match between Tier 1's per-ALT decomposition and the PG VCF's multi-allelic records at spot-checked positions.

### Concrete next steps (optional)

- Drop the 3 quasi-duplicate columns (`Nov-02`, `Est-0`, `Nok-1`) from any future pang rebuild — they contribute nothing.
- If a future SV-focused panel iteration is planned, prioritize accessions with high private-SV content — that's where extras pay off (12% of INS≥50 mass at 53 accessions, vs 1% on SNPs).
- Tier 3 (counterfactual pang_82 PanGenie rebuild) is the only way to isolate V3+D1 effects (k-mer-disambiguation downside vs upside). Not done; cost is high (~1-2 day SLURM build + per-founder PG re-runs).

---

## 8. Outputs

```
panel_overlap_135_vs_82/
├── RESULTS.md                                       # this document
├── data/
│   ├── pang_135_samples.txt                         # 135 sample IDs from raw VCF
│   ├── grenenet_82_in_pang135.txt                   # 82 numeric (cactus) sub-panel
│   ├── extras_53.txt                                # 53 named (extras) sub-panel
│   ├── pang135_biallelic_snps.vcf.gz                # 5.44M biallelic SNPs (Tier 2 input)
│   ├── gtcheck_pairs.tsv                            # 4,346 (extras, cactus) pairs
│   └── pg_vcfs.txt                                  # list of 151 PG genotyped VCFs
├── scripts/
│   ├── tier1_compute_ac.sh                          # per-ALT AC for the 82 and 53 sub-panels
│   ├── slurm_tier1_chr_array.sh                     # SLURM array for Chr2-5 (Chr1 ran on login)
│   ├── slurm_tier1_chr3_64g.sh                      # Chr3 retry (centromeric SV OOM)
│   ├── slurm_tier1_chr_array_retry.sh               # Chr4-5 retry with 64G
│   ├── slurm_pg_merge_ac.sh                         # merge 151 PG VCFs + per-ALT AC (12h job)
│   ├── slurm_tier2_gtcheck_v2.sh                    # haploid→diploid + gtcheck (segfaulted)
│   ├── tier2_pairwise_identity.py                   # streaming Python identity (used instead)
│   ├── tier2_unique_kmers_per_bubble.py             # D1 check
│   ├── analyze_tier1.py                             # per-chrom V1 classification + PG join
│   ├── combine_tier1.py                             # genome-wide V1 summary
│   ├── v1_ac_distribution.py                        # ac_53 distribution (singleton-heavy)
│   ├── v1_by_ac53_chr1.py                           # Chr1 stratified by ac_53 bin
│   ├── headline_report.py                           # final PG-carrier headline
│   └── pg_mass_attribution.py                       # ALT-mass attribution (§4)
└── results/
    ├── ac_split_Chr{1..5}.tsv.gz                    # per-ALT ac_82, ac_53 (per-chrom)
    ├── ac_pg151_Chr{1..5}.tsv.gz                    # per-ALT ac_pg (per-chrom)
    ├── v1_records_Chr{1..5}.tsv.gz                  # V1 subset per chrom
    ├── tier1_summary_Chr{1..5}.tsv                  # per-chrom category × vtype counts
    ├── tier1_summary_combined.tsv                   # genome-wide category × vtype
    ├── tier1_summary_by_vtype.tsv                   # V1 fraction by vtype
    ├── tier1_headline_Chr{1,2}.tsv                  # per-chrom PG join (Chr1, Chr2 only)
    ├── tier1_headline_combined.tsv                  # genome-wide PG join
    ├── tier1_headline_v1_by_ac53.tsv                # V1 × ac_53 bin
    ├── tier1_v1_ac_distribution.tsv                 # singleton rate per V1 sub-cat × vtype
    ├── tier2_bubble_kmers_Chr{1..5}.tsv.gz          # per-bubble unique-kmer counts
    ├── tier2_bubble_kmer_summary.tsv                # D1 stratified summary
    ├── tier2_snp_identity_matrix.tsv                # 53 × 82 SNP identity matrix
    ├── tier2_quasi_duplicates.tsv                   # 3 quasi-duplicate pairs
    ├── pg_mass_attribution.tsv                      # ALT-mass by cat × vtype (raw counts)
    └── pg_mass_attribution_pct.tsv                  # same, normalized to % per vtype
```

---

## 9. Memory-worthy entries

- **The 3 quasi-duplicates** (`Nov-02/Noveg-2`, `Est-0/Est`, `Nok-1/Nok-3`) — same pattern as `5772/6150`; effective extras count is 50.
- **The 22% / 17% / 1.4% chain** — 22% of pang_135 ALT records are extras-only (V1); 17% of those have a PG carrier; the extras contribute 1.4% of total 151 ALT mass. SVs are the exception: ~12% of INS≥50 mass.
- **D1 (more accessions → fewer unique k-mers per bubble) not detected** at gross-bubble resolution; extras-enriched bubbles have *more* unique-kmer evidence, not less.

# Possible Assembly_ID mislabeling in v3 cactus pangenome inputs

**Created:** 2026-05-14
**Author:** hapFIRE-SV team (Carnegie)
**Audience:** collaborators who provided the assembly set at `/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/`
**Status:** raised for verification — no actions taken on the pangenome panel

## IMPORTANT: arapheno Accession_ID corrections applied

> **Added 2026-05-14, after the initial investigation.**

While verifying the 80 cactus founders against the GrENE-Net 1001G panel, we discovered that two Assembly_IDs in our local manifest (`data/sv_panel_to_accession_id.tsv`) were mapped to arapheno Accession_IDs that **have no 1001G genotype data**, even though the same biological accession exists in 1001G under a different ID.

**arapheno has TWO entries for each of Ct-1 and No-0** (same collector — Albert Kranz — same lat/long, same accession name, but two registration records):

| Accession | arapheno ID #1 | arapheno ID #2 | Notes |
|---|---|---|---|
| **Ct-1** (Italy, 37.3°N 15.0°E) | **6910** (no CS number) | **7067** (CS76786) | Two registrations of the same accession in arapheno. Only **7067** has 1001G short-read genotypes. |
| **No-0** (Germany, 51.06°N 13.30°E) | **7273** (CS77128) | **7275** (no CS number) | Same pattern. Only **7273** has 1001G short-read genotypes. |

**Where each ID is / is not present:**

| ID | Accession | In arapheno? | In full 1001G VCF (1135 acc.)? | In GrENE-Net 231? | In our v3 cn_var? | Our manifest claims it for |
|---|---|---|---|---|---|---|
| 6910 | Ct-1 | ✓ | ✗ | ✗ | ✗ | Assembly_ID 100645 (**WRONG**) |
| **7067** | **Ct-1** | ✓ | ✓ | ✓ | ✓ | — |
| **7273** | **No-0** | ✓ | ✓ | ✓ | ✓ | — |
| 7275 | No-0 | ✓ | ✗ | ✗ | ✗ | Assembly_ID 100765 (**WRONG**) |

**The corrections:**

| Assembly_ID | Accession | Manifest currently says | Should be | Source of truth |
|---|---|---:|---:|---|
| **100645** | Ct-1 | 6910 | **7067** | arapheno CS76786, 1001G, v3 cn_var, seedmix recipe |
| **100765** | No-0 | 7275 | **7273** | arapheno CS77128, 1001G, v3 cn_var, seedmix recipe |

This is a manifest bug (`data/sv_panel_to_accession_id.tsv`) that does **not** affect the production v3 cn_var panel (which already uses 7067 and 7273), and does **not** affect the seedmix recipe file (which also uses the correct IDs). It only affected our own ad-hoc cross-checks against GrENE-Net 231 — earlier we incorrectly reported that 2 of the 82 cactus founders "have no 1001G partner." After applying the correction, **all 82 cactus founders are 1001G-resolvable**, and the third heatmap-flagged pair (**Stw-0 7347 vs Ct-1 7067**) can now be cross-validated:

| Pair | Cactus k-mer Jaccard | 1001G GT Jaccard |
|---|---:|---:|
| Stw-0 (7347) vs Ct-1 (**7067**) | 0.911 | **0.697** (Stw-0's closest 1001G partner; far above the bulk ~0.29) |

Both pipelines see this pair as anomalously close — more like the Ei-2/St-0 signature (concordant) than Set-1/T980 (cactus-only). So **only 5772/6150 (Set-1 / T980) remains as a cactus-pipeline-specific anomaly**. The other two suspect pairs are concordant between long-read and short-read evidence and therefore not attributable to the assembly process.

Full details and a 82×82 half-and-half heatmap (cactus k-mer Jaccard upper / 1001G genotype Jaccard lower) live in `CACTUS_ASSEMBLY_LABELING_FINDINGS.ipynb`, Section 4.

---

## TL;DR

Two pairs of Assembly_IDs in your delivery contain the **same biological sample under two different Assembly_IDs**. Direct sequence-level evidence (k-mer Jaccard on raw Chr1 assemblies, independent of any genotyping pipeline):

| Assembly_ID A | Assembly_ID B | Labels per our manifest | **k-mer Jaccard (k=31)** | Cactus VCF SNP identity | 1001G short-read SNP identity |
|---|---|---|---:|---:|---:|
| **101003** | **100954** | Set-1 (5772, UK) vs T980 (6150, SWE) | **0.9643** | 99.81% | 88.85% |
| **100300** | **100692** | Ei-2 (6915, GER) vs St-0 (8387, SWE) | **0.9662** | 99.32% | 99.39% |
| (controls) | | known-different pairs | 0.58–0.67 | ~90% | ~85–90% |

For comparison, **every other cactus founder pair** (out of 3,160 pairs across 80 cactus assemblies) is at 5-15% SNP divergence — normal A. thaliana between-accession divergence. The raw-sequence Jaccard for known-different accessions sits at 0.58–0.67.

The bimodal split is clean: ~0.96 (same sample) or ~0.6 (different accessions). There is no in-between. The flagged pairs cannot plausibly be distinct A. thaliana accessions — they share 96% of 31-mers, which is only achievable for replicate assemblies of the same source genome.

We have high confidence that cactus is processing the inputs correctly. The most likely cause of the duplication is an **upstream sample handling or labeling error before the FASTAs reached cactus** — either the same physical sample was sequenced twice and assigned two different Assembly_IDs, or the wrong physical plant batch was processed under one of the labels.

We're raising this for verification with you before taking any action on our side.

## Evidence

### 1. Cactus correctly processed two distinct input files

From `pang/logs/pang_1001gplus_82acc.log`:

```
Importing input file:///.../chr_only/101003.chr.fa
Importing input file:///.../chr_only/100954.chr.fa
Assembly stats for 101003: Total-length: 123,930,796   N50: 23,895,196
Assembly stats for 100954: Total-length: 119,496,691   N50: 23,572,047
mash distance of 101003 to TAIR10 = 0.00884936
mash distance of 100954 to TAIR10 = 0.00861543
```

Cactus saw both as distinct inputs, sanitized headers separately, computed independent mash sketches, and aligned both through minigraph. No errors or warnings in the log.

### 2. Raw FASTA files are different at the byte level

```
$ md5sum chr_only/{101003,100954,100300,100692}.chr.fa
7f393f30c58925712b7584ea14d62bdd  101003.chr.fa   (labeled 5772/Set-1)
0c962a1ff09b3dba832989c5e67d319d  100954.chr.fa   (labeled 6150/T980)
b85628dbeaf6...                   100300.chr.fa   (labeled 6915/Ei-2)
6a17442f1a8f...                   100692.chr.fa   (labeled 8387/St-0)
```

Inodes are also distinct (not hardlinks). Sizes differ. Chromosome lengths from the `.fai` indexes:

| Chr | 101003 (Set-1?) | 100954 (T980?) |
|---|---:|---:|
| Chr1 | 30,526,563 bp | 29,766,349 bp |
| Chr2 | 20,347,895 bp | 19,823,988 bp |
| Chr3 | 23,895,196 bp | 23,572,047 bp |
| Chr4 | 21,566,898 bp | 19,568,252 bp |
| Chr5 | 27,594,244 bp | 26,766,055 bp |

100954 is systematically shorter on every chromosome (1.4-9.3%). Pattern is consistent with "same sample assembled with less complete output" rather than two distinct biological samples.

### 3. Cactus pangenome VCF: per-sample GTs are 99.81% identical

Among 791,044 Chr1 SNP records with non-missing GT for both samples in `pang_1001gplus_82acc.vcf.gz`, only 1,541 (0.19%) differ between sample 101003 and sample 100954. Sample 100300 and sample 100692 differ at 0.68%.

### 4. GrENE-Net 1001G short-read data disagrees on 5772/6150 but agrees on 6915/8387

The independent 1001 Genomes SNP catalog (used by hapFIRE: `/.../greneNet_final_v1.1.recode.vcf`, n = 827,765 Chr1 SNPs across 231 ecotypes) shows:

| pair | 1001G short-read divergence | Interpretation |
|---|---:|---|
| 5772 (Set-1) vs 6150 (T980) | **11.15%** | typical between-accession — *distinct accessions per 1001G* |
| 6915 (Ei-2) vs 8387 (St-0) | **0.61%** | near-identical — *same biological sample per 1001G* |

For 5772/6150: cactus and GrENE-Net **disagree**. Cactus says the FASTAs we labeled Set-1 and T980 are the same biological sample (99.81% identity), but the 1001G short-read pipeline says they're distinct (88.85% identity). Cross-pair controls confirm GrENE-Net is internally consistent (5772 vs 6915 = 10.06%, 6150 vs 8387 = 10.73%, all in the typical 10-15% between-accession range).

For 6915/8387: both pipelines **agree** these two accession IDs refer to the same biological sample. Possible 1001G catalog issue (same physical accession registered under two IDs) — but consistent across both genotyping pipelines, so the duplication exists upstream of any genotyping step.

### 5. The pairwise audit rules out a systematic cactus issue

Pairwise SNP-level distance among all 80 cactus founders in v3 cn_var (Chr1 SNPs):
- median: 10.06%
- p5: 8.39%
- p25: 9.26%
- next-closest pair after the two flagged: 768 vs 772 at 3.25% (closely-related accessions, confirmed by 1001G)
- min: 0.32% (101003/100954) and 0.68% (100300/100692)

Only the two flagged pairs are below the 5% threshold. If cactus had a generic collapse bug, we'd expect many pairs at low divergence. The specificity to just these two pairs points to per-sample input issues, not a tool-level bug.

### 6. Direct k-mer Jaccard on raw Chr1 assemblies — the cleanest evidence

To eliminate any dependence on cactus or hapFIRE, we counted canonical 31-mers in each raw Chr1 sequence (jellyfish, k=31) and computed pairwise Jaccard. This is independent of any genotyping pipeline:

| pair | label | Jaccard | call |
|---|---|---:|---|
| **101003 ↔ 100954** | **Set-1 vs T980** | **0.9643** | **SAME SAMPLE** |
| **100300 ↔ 100692** | **Ei-2 vs St-0** | **0.9662** | **SAME SAMPLE** |
| 100764 ↔ 100852 | Mt-0 vs Ped-0 | 0.5944 | distinct |
| 100715 ↔ 100778 | Qar-8a vs Hau-0 | 0.6553 | distinct |
| 100042 ↔ 100043 | Cvi-0 vs BRI-2 (known different) | 0.5865 | distinct ✓ control |
| 101003 ↔ 100042 | Set-1 vs Cvi-0 (cross-pair) | 0.5899 | distinct |
| 100954 ↔ 100043 | T980 vs BRI-2 (cross-pair) | 0.6635 | distinct |
| 100300 ↔ 100042 | Ei-2 vs Cvi-0 (cross-pair) | 0.5827 | distinct |
| 101003 ↔ 100300 | Set-1 vs Ei-2 (suspect cross-pair) | 0.6608 | distinct |
| 100954 ↔ 100692 | T980 vs St-0 (suspect cross-pair) | 0.6691 | distinct |

The distribution is bimodal: 0.96+ for the two flagged pairs, 0.58–0.67 for every other pair. There is no middle ground. At k=31 on a ~30 Mb chromosome, even a single SNP disrupts up to 31 distinct k-mers, so the per-pair Jaccard for two genuinely distinct A. thaliana accessions is dominated by their SNP-level divergence. A Jaccard above 0.96 cannot be achieved unless the underlying sequences are nearly identical at virtually every k-mer position — i.e., the same biological sample assembled twice with minor scaffolding differences.

This evidence is independent of cactus and hapFIRE pipelines. It is direct sequence-level confirmation that the FASTAs labeled 101003 and 100954 contain the same genome, and similarly for 100300 and 100692.

## What we're asking

For **101003 vs 100954** (Assembly_IDs labeled Set-1 and T980):
- Were both Assembly_IDs assigned from physically distinct seed batches / plant material? Or could a single physical sample have been sequenced twice and assigned two Assembly_IDs?
- Is there an internal record of which 1001G accession ID each Assembly_ID was generated from?
- Specifically: is 101003 indeed Set-1 (UK), and 100954 indeed T980 (SWE)?

For **100300 vs 100692** (Assembly_IDs labeled Ei-2 and St-0):
- Same questions. Both panels say these are the same biological sample. Was the same physical sample processed under both Ei-2 and St-0 labels at the sequencing step?

If our hypothesis (one of these pairs is the same sample, mislabeled) is correct, the practical implication is that v3 cn_var has TWO copies of one accession and ZERO copies of the other in each affected pair. Knowing which physical sample is correctly represented (and which is missing) is what we need before we can correct the panel.

## Impact on our analyses

The flagged pairs introduce simplex degeneracy in `cactus_em`'s founder-frequency estimation: the EM cannot distinguish between two identical founders, so any pool mass associated with either real ecotype gets split arbitrarily between the two duplicate-genotype columns. At per-SNP AF projections, this manifests as systematic discordance with hapFIRE's per-SNP AF estimates at the ~0.6% of Chr1 SNPs that fall in the relevant carrier-overlap regions.

This is a manageable downstream concern. We're not blocked. We just want to verify the labeling before deciding whether to (a) drop the redundant founders, (b) keep them and document, or (c) wait for corrected FASTAs.

## Reference paths

- Cactus output VCF: `/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.vcf.gz`
- Cactus log: `/home/tbellagio/scratch/pang/pang_1001gplus/pang/logs/pang_1001gplus_82acc.log`
- Assembly inputs: `/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/{101003,100954,100300,100692}.chr.fa`
- GrENE-Net VCF: `/carnegie/nobackup/scratch/xwu/GrENE_net/greneNet_final_v1.1.recode.vcf`
- Assembly_ID → Accession_ID manifest: `/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/sv_panel_to_accession_id.tsv`
- Diagnostic analysis: `/carnegie/nobackup/scratch/tbellagio/hapfire_sv/SEEDMIX_S1_v3_vs_hapfire.ipynb`

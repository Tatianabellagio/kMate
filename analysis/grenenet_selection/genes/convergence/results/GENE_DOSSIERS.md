# Gene-by-gene review — themed shortlist

23 genes, one dossier each.

⚠ **The `verdict` column is a 3-criterion count, not a judgement.** It counts whether the gene is recurrent (>=3 climate clusters or gardens, or found by both scans), common (>=20 of 231 founders carry the lead) and directly placed (CDS/UTR/promoter). Every threshold in that is arbitrary, and it visibly mis-ranks: CRK18 -- cross-scan, frameshift, r2=0.87 to its own gene, in the one locus that survived honest WZA recalibration -- is demoted to MODERATE purely by having 15 carriers instead of 20, while genes with r2~0 to their gene are promoted for being common. Use the component columns and the figures; treat the label as a sort key.

**LD is deliberately not in the verdict.** Every gene here is assigned at its variant's own position, so all leads sit physically inside their target gene and CARK-type detachment is impossible by construction. An earlier version of this script ranked on LD and graded GPX6 WEAK at r2=0.24, where `dissection/results/loci/CANDIDATE_VERDICTS.md` graded the same variant at the same r2 ROBUST -- correctly, because a promoter variant on its own haplotype background is a better-localized causal candidate, not a worse one. `ld_ambiguous` flags where a tightly-linked alternative could carry the signal instead: that is about which VARIANT, not which gene.

| gene | criteria met | why | r²(lead↔gene) | best local r² | carriers | themes | found by | locus (n genes) |
|---|---|---|---|---|---|---|---|---|
| CYP71B4 | **ROBUST** | recurrent, common (73 carriers), direct placement (F4_promoter) | 0.61 | 0.89 | 73 | climate_temp | GWAS | L0097 (36) |
| EPFL5 | **ROBUST** | recurrent, common (22 carriers), direct placement (F4_promoter) | 0.39 | 0.70 | 22 | stress | GEA+GWAS | L0096 (20) |
| ABHD11 | **ROBUST** | recurrent, common (20 carriers), direct placement (F4_promoter) | 0.25 ⚠ | 1.00 | 20 | climate_temp | GEA+GWAS | L0123 (7) |
| GPX6 | **ROBUST** | recurrent, common (27 carriers), direct placement (F4_promoter) | 0.24 ⚠ | 0.96 | 27 | stress | GEA | L0124 (65) |
| MJB24.8 | **ROBUST** | recurrent, common (69 carriers), direct placement (F4_promoter) | 0.22 | 0.22 | 69 | climate_temp | GEA | L0173 (35) |
| AT2G29210 | **ROBUST** | recurrent, common (70 carriers), direct placement (F2_CDS_inframe) | 0.16 ⚠ | 1.00 | 70 | stress | GEA+GWAS | L0067 (4) |
| GAPC2 | **MODERATE** | not recurrent, common (40 carriers), direct placement (F4_promoter) | 1.00 | 1.00 | 40 | climate_temp,stress | GEA | L0008 (28) |
| SD17 | **MODERATE** | recurrent, common (40 carriers), indirect placement (F6_intron) | 0.97 | 0.97 | 40 | stress | GEA | L0038 (7) |
| EMB1241 | **MODERATE** | not recurrent, common (27 carriers), direct placement (F3_UTR) | 0.92 | 0.92 | 27 | climate_temp | GEA | L0151 (20) |
| CRK18 | **MODERATE** | recurrent, rare (15 carriers), direct placement (F1_CDS_frameshift) | 0.87 | 0.93 | 15 | stress | GEA+GWAS | L0128 (19) |
| PLDALPHA2 | **MODERATE** | not recurrent, common (38 carriers), direct placement (F4_promoter) | 0.85 | 0.88 | 38 | stress | GEA | L0032 (31) |
| RALFL28 | **MODERATE** | recurrent, rare (12 carriers), direct placement (F1_CDS_frameshift) | 0.74 | 1.00 | 12 | stress | GEA | L0124 (65) |
| AT1G52180 | **MODERATE** | not recurrent, common (28 carriers), direct placement (F4_promoter) | 0.63 | 1.00 | 28 | stress | GEA | L0032 (31) |
| AIPP3 | **MODERATE** | recurrent, rare (11 carriers), direct placement (F4_promoter) | 0.43 | 0.78 | 11 | circadian_light,flowering | GEA | L0124 (65) |
| DOB1 | **MODERATE** | recurrent, rare (6 carriers), direct placement (F4_promoter) | 0.07 | 0.07 | 6 | stress | GWAS | L0130 (1) |
| BT4 | **MODERATE** | recurrent, rare (11 carriers), direct placement (F3_UTR) | 0.07 | 0.16 | 11 | stress | GEA+GWAS | L0178 (12) |
| GRXS2 | **MODERATE** | recurrent, rare (6 carriers), direct placement (F4_promoter) | 0.00 | 0.30 | 6 | stress | GWAS | L0151 (20) |
| NIP1-1 | **WEAK** | not recurrent, common (90 carriers), indirect placement (nan) | 0.98 | 0.98 | 90 | stress | GEA | L0127 (16) |
| GSH1 | **WEAK** | not recurrent, rare (11 carriers), direct placement (F3_UTR) | 0.45 | 0.46 | 11 | climate_temp,flowering,stress | GEA | L0128 (19) |
| GLR1.3 | **WEAK** | not recurrent, rare (5 carriers), direct placement (F3_UTR) | 0.25 | 0.32 | 5 | circadian_light | GWAS | L0169 (7) |
| APR2 | **WEAK** | not recurrent, rare (7 carriers), direct placement (F4_promoter) | 0.12 | 0.27 | 7 | stress | GEA | L0036 (10) |
| UPM1 | **FRAGILE** | not recurrent, rare (11 carriers), indirect placement (F6_intron) | 0.60 | 0.75 | 11 | circadian_light,climate_temp | GEA | L0166 (14) |
| PRT6 | **FRAGILE** | not recurrent, rare (5 carriers), indirect placement (F6_intron) | 0.12 | 0.25 | 5 | climate_temp,stress | GWAS | L0140 (9) |

---

## CYP71B4 (AT3G26280) — **ROBUST**

*recurrent, common (73 carriers), direct placement (F4_promoter).*

- **Protein** — Cytochrome P450 71B4 (EC 1.14.-.-)
- **Themes** — climate_temp
- **Location** — Chr3:9,630,200–9,632,009 (-) · locus `L0097` holding **36** candidate genes
- **Found by** — GWAS · GWAS 3 gardens (12,55,60), nlp 7.6, MAC 6
- **Evidence** — 1 independent line(s): GWAS gardens
- **Blocks** — Chr3_5089  ⚠ block attribution would name a different gene
- **Variants** — 1 (1 SV / 0 indel / 0 MNP), max 2421 bp, in promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr3:9632650 (sv, 2421 bp, promoter), 73 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.61** / med 0.08 · best local r² 0.89
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² not testable
- **Keywords** — Heme;Iron;Membrane;Metal-binding;Monooxygenase;Oxidoreductase;Reference proteome;Transmembrane;Transmembrane helix

![CYP71B4](plots/loci/CYP71B4_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr3:9,632,650 | sv | 2421 | promoter | F4_promoter | — | 0 | 3 |

---

## EPFL5 (AT3G22820) — **ROBUST**

*recurrent, common (22 carriers), direct placement (F4_promoter).*

- **Protein** — EPIDERMAL PATTERNING FACTOR-like protein 5 (EPF-like protein 5) [Cleaved into: CHALLAH-LIKE1]
- **Themes** — stress
- **Location** — Chr3:8,073,271–8,074,138 (-) · locus `L0096` holding **20** candidate genes
- **Found by** — GEA+GWAS · GEA 1 climate clusters / 1 axes (bio17), best bio15 nlp 9.9, λ_min 2.67 · GWAS 4 gardens (10,12,27,42), nlp 7.8, MAC 5
- **Evidence** — 2 independent line(s): GWAS gardens, both scans
- **Blocks** — Chr3_3955
- **Variants** — 1 (0 SV / 1 indel / 0 MNP), max 10 bp, in promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr3:8074451 (smallindel, 10 bp, promoter), 22 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.39** / med 0.02 · best local r² 0.70
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² 0.77
- **Function (UniProt)** — Controls stomatal patterning. Mediates differentiation of stomatal lineage cells to pavement cells and stomatal development inhibition (PubMed:23748792). TMM (AC Q9SSD1) functions to dampen or block CLL1 signaling. Acts as a growth-regulatory ligand for ERECTA family receptors. Promotes fruit growth and fertility (PubMed:22474391). {ECO:0000269|PubMed:19435754, ECO:0000269|PubMed:21862708, ECO:0000269|PubMed:22474391, ECO:0000269|PubMed:23748792}.
- **Keywords** — Developmental protein;Disulfide bond;Reference proteome;Secreted;Signal

![EPFL5](plots/loci/EPFL5_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr3:8,074,451 | smallindel | 10 | promoter | F4_promoter | 9.9 | 1 | 4 |

---

## ABHD11 (AT4G10030) — **ROBUST**

*recurrent, common (20 carriers), direct placement (F4_promoter).*

- **Protein** — Uncharacterized protein
- **Themes** — climate_temp
- **Location** — Chr4:6,270,241–6,273,241 (-) · locus `L0123` holding **7** candidate genes
- **Found by** — GEA+GWAS · GEA 1 climate clusters / 1 axes (bio17), best bio15 nlp 7.2, λ_min 2.67 · GWAS 1 gardens (55), nlp 7.9, MAC 5
- **Evidence** — 1 independent line(s): both scans
- **Blocks** — Chr4_2226
- **Variants** — 2 (1 SV / 1 indel / 0 MNP), max 7341 bp, in promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr4:6273895 (smallindel, 1 bp, promoter), 20 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.25** / med 0.01 · best local r² 1.00  ⚠ **LD-ambiguous**: a tightly-linked alternative in the window could carry the signal instead (about which VARIANT, not which gene)
- **SNPs** — a Bonferroni SNP within 2 kb · tagging r² not testable

![ABHD11](plots/loci/ABHD11_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr4:6,273,895 | smallindel | 1 | promoter | F4_promoter | 7.2 | 1 | 1 |
| Chr4:6,273,282 | sv | 7341 | promoter | F4_promoter | — | 1 | 1 |

---

## GPX6 (AT4G11600) — **ROBUST**

*recurrent, common (27 carriers), direct placement (F4_promoter).*

- **Protein** — Probable phospholipid hydroperoxide glutathione peroxidase 6, mitochondrial (AtGPX1) (PHGPx) (EC 1.11.1.12)
- **Themes** — stress
- **Location** — Chr4:7,009,769–7,011,375 (-) · locus `L0124` holding **65** candidate genes
- **Found by** — GEA · GEA 5 climate clusters / 6 axes (bio1,bio16,bio17,bio3,bio7), best pc1 nlp 10.2, λ_min 1.37
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr4_2749
- **Variants** — 1 (1 SV / 0 indel / 0 MNP), max 1164 bp, in promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr4:7011705 (sv, 1164 bp, promoter), 27 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.24** / med 0.10 · best local r² 0.96  ⚠ **LD-ambiguous**: a tightly-linked alternative in the window could carry the signal instead (about which VARIANT, not which gene)
- **SNPs** — a Bonferroni SNP within 2 kb · tagging r² not testable
- **Function (UniProt)** — Protects cells and enzymes from oxidative damage, by catalyzing the reduction of hydrogen peroxide, lipid peroxides and organic hydroperoxide, by glutathione. {ECO:0000250|UniProtKB:O70325}.
- **Keywords** — Mitochondrion;Oxidoreductase;Peroxidase;Reference proteome;Stress response;Transit peptide

![GPX6](plots/loci/GPX6_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr4:7,011,705 | sv | 1164 | promoter | F4_promoter | 10.2 | 5 | 0 |

---

## MJB24.8 (AT5G57270) — **ROBUST**

*recurrent, common (69 carriers), direct placement (F4_promoter).*

- **Protein** — Core-2/I-branching beta-1,6-N-acetylglucosaminyltransferase family protein
- **Themes** — climate_temp
- **Location** — Chr5:23,200,614–23,203,792 (-) · locus `L0173` holding **35** candidate genes
- **Found by** — GEA · GEA 3 climate clusters / 4 axes (bio1,bio17,bio3,pc3), best pc3 nlp 9.2, λ_min 1.99
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr5_12000,Chr5_12002
- **Variants** — 9 (0 SV / 5 indel / 4 MNP), max 4 bp, in intron,promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr5:23204156 (smallindel, 2 bp, promoter), 69 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.22** / med 0.18 · best local r² 0.22
- **SNPs** — a Bonferroni SNP within 2 kb · tagging r² 1.00
- **Keywords** — Reference proteome

![MJB24.8](plots/loci/MJB24.8_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr5:23,204,156 | smallindel | 2 | promoter | F4_promoter | 9.2 | 3 | 0 |
| Chr5:23,203,984 | smallindel | 2 | promoter | F4_promoter | 8.1 | 3 | 0 |
| Chr5:23,202,200 | smallindel | 4 | intron | F6_intron | 7.8 | 3 | 0 |
| Chr5:23,201,979 | smallindel | 1 | intron | F6_intron | 7.4 | 3 | 0 |
| Chr5:23,201,721 | smallindel | 1 | intron | F6_intron | 7.2 | 3 | 0 |

---

## AT2G29210 (AT2G29210) — **ROBUST**

*recurrent, common (70 carriers), direct placement (F2_CDS_inframe).*

- **Protein** — Splicing factor PWI domain-containing protein
- **Themes** — stress
- **Location** — Chr2:12,558,051–12,562,348 (+) · locus `L0067` holding **4** candidate genes
- **Found by** — GEA+GWAS · GEA 1 climate clusters / 1 axes (bio17), best bio15 nlp 9.0, λ_min 2.67 · GWAS 1 gardens (55), nlp 7.6, MAC 28
- **Evidence** — 1 independent line(s): both scans
- **Blocks** — Chr2_4974
- **Variants** — 2 (0 SV / 2 indel / 0 MNP), max 30 bp, in CDS,promoter; strongest mechanism **F2_CDS_inframe**
- **Lead** — Chr2:12557749 (smallindel, 9 bp, promoter), 70 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.16** / med 0.07 · best local r² 1.00  ⚠ **LD-ambiguous**: a tightly-linked alternative in the window could carry the signal instead (about which VARIANT, not which gene)
- **SNPs** — a Bonferroni SNP within 2 kb · tagging r² 1.00
- **Keywords** — mRNA processing;Proteomics identification;Reference proteome

![AT2G29210](plots/loci/AT2G29210_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr2:12,557,749 | smallindel | 9 | promoter | F4_promoter | 9.0 | 1 | 1 |
| Chr2:12,560,442 | smallindel | 30 | CDS | F2_CDS_inframe | — | 1 | 1 |

---

## GAPC2 (AT1G13440) — **MODERATE**

*not recurrent, common (40 carriers), direct placement (F4_promoter).*

- **Protein** — Glyceraldehyde-3-phosphate dehydrogenase GAPC2, cytosolic (EC 1.2.1.12) (NAD-dependent glyceraldehydephosphate dehydrogenase C subunit 2)
- **Themes** — climate_temp,stress
- **Location** — Chr1:4,608,196–4,610,647 (-) · locus `L0008` holding **28** candidate genes
- **Found by** — GEA · GEA 2 climate clusters / 3 axes (bio1,bio17), best bio15 nlp 9.3, λ_min 1.99
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr1_2343,Chr1_2345
- **Variants** — 9 (0 SV / 7 indel / 2 MNP), max 24 bp, in intron,promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr1:4609984 (smallindel, 3 bp, intron), 40 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **1.00** / med 0.52 · best local r² 1.00
- **SNPs** — a Bonferroni SNP within 2 kb · tagging r² 1.00
- **Function (UniProt)** — Key enzyme in glycolysis that catalyzes the first step of the pathway by converting D-glyceraldehyde 3-phosphate (G3P) into 3-phospho-D-glyceroyl phosphate. Essential for the maintenance of cellular ATP levels and carbohydrate metabolism (By similarity). Binds DNA in vitro (PubMed:22589465). Together with DNA polymerase II subunit B3-1 (DPB3-1) and GAPC1, enhances heat tolerance and promotes the expression of heat-inducible genes (PubMed:32651385). {ECO:0000250|UniProtKB:P25858, ECO:0000269|PubMed:22589465, ECO:0000269|PubMed:32651385}.
- **Keywords** — Alternative splicing;Cytoplasm;DNA-binding;Glutathionylation;Glycolysis;NAD;Nucleus;Oxidoreductase;Reference proteome;S-nitrosylation

![GAPC2](plots/loci/GAPC2_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr1:4,610,015 | smallindel | 24 | intron | F6_intron | 9.3 | 2 | 0 |
| Chr1:4,609,984 | smallindel | 3 | intron | F6_intron | 9.3 | 2 | 0 |
| Chr1:4,609,651 | smallindel | 15 | intron | F6_intron | 9.1 | 2 | 0 |
| Chr1:4,609,623 | smallindel | 2 | intron | F6_intron | 9.1 | 2 | 0 |
| Chr1:4,609,611 | smallindel | 3 | intron | F6_intron | 9.1 | 2 | 0 |
| Chr1:4,609,814 | smallindel | 1 | intron | F6_intron | 9.0 | 2 | 0 |
| Chr1:4,611,578 | smallindel | 1 | promoter | F4_promoter | 7.8 | 2 | 0 |

---

## SD17 (AT1G65790) — **MODERATE**

*recurrent, common (40 carriers), indirect placement (F6_intron).*

- **Protein** — Receptor-like serine/threonine-protein kinase SD1-7 (EC 2.7.11.1) (Arabidopsis thaliana receptor kinase 1) (S-domain-1 (SD1) receptor kinase 7) (SD1-7)
- **Themes** — stress
- **Location** — Chr1:24,468,932–24,472,329 (+) · locus `L0038` holding **7** candidate genes
- **Found by** — GEA · GEA 3 climate clusters / 3 axes (bio1,bio17,bio3), best bio3 nlp 9.4, λ_min 1.89
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr1_12318
- **Variants** — 1 (1 SV / 0 indel / 0 MNP), max 283 bp, in intron; strongest mechanism **F6_intron**
- **Lead** — Chr1:24470448 (sv, 283 bp, intron), 40 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.97** / med 0.22 · best local r² 0.97
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² not testable
- **Function (UniProt)** — Involved in the regulation of cellular expansion and differentiation. Mediates subcellular relocalization of PUB9 from nucleus to plasma membrane in a protein-phosphorylation-dependent manner. May be involved in the abscisic acid-mediated signaling pathway, at least during germination. {ECO:0000269|PubMed:18552232, ECO:0000269|PubMed:8811866}.
- **Keywords** — Abscisic acid signaling pathway;Alternative splicing;ATP-binding;Cell membrane;Disulfide bond;EGF-like domain;Glycoprotein;Kinase;Lectin;Membrane;Nucleotide-binding;Phosphoprotein;Receptor;Reference proteome;Serine/threonine-protein kinase;Signal;Transferase;Transmembrane;Transmembrane helix

![SD17](plots/loci/SD17_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr1:24,470,448 | sv | 283 | intron | F6_intron | 9.4 | 3 | 0 |

---

## EMB1241 (AT5G17710) — **MODERATE**

*not recurrent, common (27 carriers), direct placement (F3_UTR).*

- **Protein** — Co-chaperone GrpE family protein
- **Themes** — climate_temp
- **Location** — Chr5:5,839,343–5,841,728 (-) · locus `L0151` holding **20** candidate genes
- **Found by** — GEA · GEA 2 climate clusters / 7 axes (bio1,bio17), best bio15 nlp 11.9, λ_min 1.62
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr5_2856
- **Variants** — 1 (0 SV / 1 indel / 0 MNP), max 2 bp, in 3'UTR; strongest mechanism **F3_UTR**
- **Lead** — Chr5:5839507 (smallindel, 2 bp, 3'UTR), 27 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.92** / med 0.02 · best local r² 0.92
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² 1.00
- **Keywords** — Reference proteome

![EMB1241](plots/loci/EMB1241_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr5:5,839,507 | smallindel | 2 | 3'UTR | F3_UTR | 11.9 | 2 | 0 |

---

## CRK18 (AT4G23260) — **MODERATE**

*recurrent, rare (15 carriers), direct placement (F1_CDS_frameshift).*

- **Protein** — Cysteine-rich receptor-like protein kinase 18 (Cysteine-rich RLK18) (EC 2.7.11.-)
- **Themes** — stress
- **Location** — Chr4:12,167,354–12,170,086 (-) · locus `L0128` holding **19** candidate genes
- **Found by** — GEA+GWAS · GEA 2 climate clusters / 2 axes (bio1,bio7), best bio6 nlp 7.3, λ_min 1.80 · GWAS 1 gardens (12), nlp 7.8, MAC 7
- **Evidence** — 2 independent line(s): GEA recurrence, both scans
- **Blocks** — Chr4_6317,Chr4_6319
- **Variants** — 2 (0 SV / 2 indel / 0 MNP), max 2 bp, in CDS,promoter; strongest mechanism **F1_CDS_frameshift** · **frameshift**
- **Lead** — Chr4:12170418 (smallindel, 1 bp, promoter), 15 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.87** / med 0.06 · best local r² 0.93
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² 1.00
- **Keywords** — Alternative splicing;ATP-binding;Glycoprotein;Kinase;Membrane;Nucleotide-binding;Phosphoprotein;Receptor;Reference proteome;Repeat;Serine/threonine-protein kinase;Signal;Transferase;Transmembrane;Transmembrane helix

![CRK18](plots/loci/CRK18_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr4:12,170,418 | smallindel | 1 | promoter | F4_promoter | 7.3 | 2 | 1 |
| Chr4:12,169,155 | smallindel | 2 | CDS | F1_CDS_frameshift | — | 2 | 1 |

---

## PLDALPHA2 (AT1G52570) — **MODERATE**

*not recurrent, common (38 carriers), direct placement (F4_promoter).*

- **Protein** — Phospholipase D alpha 2 (AtPLDalpha2) (PLD alpha 2) (EC 3.1.4.4) (Choline phosphatase 2) (Phosphatidylcholine-hydrolyzing phospholipase D 2)
- **Themes** — stress
- **Location** — Chr1:19,583,940–19,587,050 (-) · locus `L0032` holding **31** candidate genes
- **Found by** — GEA · GEA 2 climate clusters / 3 axes (bio1,pc3), best pc3 nlp 11.2, λ_min 1.99
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr1_8914,Chr1_8915,Chr1_8916  ⚠ block attribution would name a different gene
- **Variants** — 6 (0 SV / 5 indel / 1 MNP), max 24 bp, in intron,promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr1:19587156 (smallindel, 1 bp, promoter), 38 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.85** / med 0.08 · best local r² 0.88
- **SNPs** — a Bonferroni SNP within 2 kb · tagging r² 1.00
- **Function (UniProt)** — Hydrolyzes glycerol-phospholipids at the terminal phosphodiesteric bond to generate phosphatidic acids (PA). Plays an important role in various cellular processes, including phytohormone action and response to stress, characterized by acidification of the cell. {ECO:0000250|UniProtKB:Q38882}.
- **Keywords** — Abscisic acid signaling pathway;Calcium;Cytoplasm;Cytoplasmic vesicle;Ethylene signaling pathway;Hydrolase;Lipid degradation;Lipid metabolism;Membrane;Metal-binding;Reference proteome;Repeat;Vacuole

![PLDALPHA2](plots/loci/PLDALPHA2_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr1:19,587,156 | smallindel | 1 | promoter | F4_promoter | 11.2 | 2 | 0 |
| Chr1:19,587,160 | smallindel | 1 | promoter | F4_promoter | 8.3 | 2 | 0 |
| Chr1:19,586,622 | smallindel | 24 | intron | F6_intron | 8.2 | 2 | 0 |
| Chr1:19,586,584 | smallindel | 8 | intron | F6_intron | 8.2 | 2 | 0 |
| Chr1:19,587,166 | smallindel | 2 | promoter | F4_promoter | 7.8 | 2 | 0 |

---

## RALFL28 (AT4G11510) — **MODERATE**

*recurrent, rare (12 carriers), direct placement (F1_CDS_frameshift).*

- **Protein** — Protein RALF-like 28
- **Themes** — stress
- **Location** — Chr4:6,984,051–6,984,308 (-) · locus `L0124` holding **65** candidate genes
- **Found by** — GEA · GEA 4 climate clusters / 6 axes (bio1,bio17,bio3,bio7), best bio6 nlp 11.0, λ_min 1.59
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr4_2725  ⚠ block attribution would name a different gene
- **Variants** — 13 (1 SV / 10 indel / 2 MNP), max 1092 bp, in CDS,intergenic,promoter; strongest mechanism **F1_CDS_frameshift** · **frameshift**
- **Lead** — Chr4:6984627 (smallindel, 17 bp, promoter), 12 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.74** / med 0.10 · best local r² 1.00
- **SNPs** — a Bonferroni SNP within 2 kb · tagging r² not testable
- **Function (UniProt)** — Cell signaling peptide that may regulate plant stress, growth, and development. Mediates a rapid alkalinization of extracellular space by mediating a transient increase in the cytoplasmic Ca(2+) concentration leading to a calcium-dependent signaling events through a cell surface receptor and a concomitant activation of some intracellular mitogen-activated protein kinases (By similarity). {ECO:0000250}.
- **Keywords** — Disulfide bond;Hormone;Reference proteome;Secreted;Signal

![RALFL28](plots/loci/RALFL28_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr4:6,984,627 | smallindel | 17 | promoter | F4_promoter | 11.0 | 3 | 0 |
| Chr4:6,984,722 | smallindel | 1 | promoter | F4_promoter | 11.0 | 3 | 0 |
| Chr4:6,984,655 | smallindel | 3 | promoter | F4_promoter | 11.0 | 3 | 0 |
| Chr4:6,982,683 | sv | 1092 | intergenic | F5_downstream | 10.9 | 3 | 0 |
| Chr4:6,984,038 | smallindel | 1 | intergenic | F5_downstream | 9.7 | 3 | 0 |
| Chr4:6,985,057 | smallindel | 1 | promoter | F4_promoter | 8.5 | 3 | 0 |
| Chr4:6,985,065 | smallindel | 2 | promoter | F4_promoter | 8.5 | 3 | 0 |
| Chr4:6,983,993 | smallindel | 1 | intergenic | F5_downstream | 8.0 | 3 | 0 |
| Chr4:6,984,382 | smallindel | 4 | promoter | F4_promoter | 7.5 | 3 | 0 |
| Chr4:6,984,361 | smallindel | 3 | promoter | F4_promoter | 7.5 | 3 | 0 |
| Chr4:6,984,052 | smallindel | 1 | CDS | F1_CDS_frameshift | 7.4 | 3 | 0 |

---

## AT1G52180 (AT1G52180) — **MODERATE**

*not recurrent, common (28 carriers), direct placement (F4_promoter).*

- **Protein** — Aquaporin-like superfamily protein
- **Themes** — stress
- **Location** — Chr1:19,424,944–19,425,928 (-) · locus `L0032` holding **31** candidate genes
- **Found by** — GEA · GEA 2 climate clusters / 6 axes (bio17,bio7), best pc1 nlp 7.7, λ_min 1.57
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr1_8825
- **Variants** — 1 (1 SV / 0 indel / 0 MNP), max 680 bp, in promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr1:19426362 (sv, 680 bp, promoter), 28 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.63** / med 0.10 · best local r² 1.00
- **SNPs** — a Bonferroni SNP within 2 kb · tagging r² 1.00
- **Keywords** — Membrane;Reference proteome;Repeat;Transmembrane;Transmembrane helix;Transport

![AT1G52180](plots/loci/AT1G52180_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr1:19,426,362 | sv | 680 | promoter | F4_promoter | 7.7 | 2 | 0 |

---

## AIPP3 (AT4G11560) — **MODERATE**

*recurrent, rare (11 carriers), direct placement (F4_promoter).*

- **Protein** — ASI1-immunoprecipitated protein 3 (Bromo-adjacent homology domain-containing protein 1) (BAH domain-containing transcriptional regulator 1) (Protein REPRESSOR OF VERNALIZATION 1) (AtRVR1)
- **Themes** — circadian_light,flowering
- **Location** — Chr4:6,999,896–7,003,493 (-) · locus `L0124` holding **65** candidate genes
- **Found by** — GEA · GEA 4 climate clusters / 5 axes (bio1,bio17,bio3,bio7,pc3), best bio6 nlp 10.9, λ_min 1.80
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr4_2730,Chr4_2733
- **Variants** — 4 (0 SV / 4 indel / 0 MNP), max 4 bp, in intron,promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr4:7003757 (smallindel, 4 bp, promoter), 11 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.43** / med 0.00 · best local r² 0.78
- **SNPs** — a Bonferroni SNP within 2 kb · tagging r² 0.51
- **Function (UniProt)** — Transcriptional repressor (PubMed:33433058). Together with PHD finger-containing proteins (e.g. PHD1, PAIPP2/PHD2, AIPP2/PHD3, PHD4, PHD5 and PHD6), cooperates to form a BAH-PHD bivalent histone reader complex able to read histone H3 lysine 27 trimethylation (H3K27me3) and low-methylated H3K4 histone marks in order to regulate transcription, especially to prevent early flowering; H3K27me3 reader of this complex (PubMed:33277495, PubMed:33433058). CPL2 is subsequently recruited to form a BAH-PHD-CPL2 complex (BPC) in order to silence several H3K27me3 and low-methylated H3K4 enriched loci, inclu
- **Keywords** — 3D-structure;Nucleus;Reference proteome;Repressor;Transcription;Transcription regulation

![AIPP3](plots/loci/AIPP3_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr4:7,003,757 | smallindel | 4 | promoter | F4_promoter | 10.9 | 4 | 0 |
| Chr4:7,003,747 | smallindel | 1 | promoter | F4_promoter | 8.4 | 4 | 0 |
| Chr4:7,003,703 | smallindel | 1 | promoter | F4_promoter | 7.7 | 4 | 0 |
| Chr4:7,002,464 | smallindel | 1 | intron | F6_intron | 7.2 | 4 | 0 |

---

## DOB1 (AT4G25670) — **MODERATE**

*recurrent, rare (6 carriers), direct placement (F4_promoter).*

- **Protein** — Stress response NST1-like protein
- **Themes** — stress
- **Location** — Chr4:13,085,198–13,087,090 (-) · locus `L0130` holding **1** candidate genes
- **Found by** — GWAS · GWAS 3 gardens (27,42,55), nlp 8.5, MAC 6
- **Evidence** — 1 independent line(s): GWAS gardens
- **Blocks** — Chr4_6746
- **Variants** — 1 (0 SV / 1 indel / 0 MNP), max 1 bp, in promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr4:13087216 (smallindel, 1 bp, promoter), 6 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.07** / med 0.01 · best local r² 0.07
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² 0.59
- **Keywords** — Reference proteome

![DOB1](plots/loci/DOB1_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr4:13,087,216 | smallindel | 1 | promoter | F4_promoter | — | 0 | 3 |

---

## BT4 (AT5G67480) — **MODERATE**

*recurrent, rare (11 carriers), direct placement (F3_UTR).*

- **Protein** — BTB/POZ and TAZ domain-containing protein 4 (BTB and TAZ domain protein 4)
- **Themes** — stress
- **Location** — Chr5:26,930,939–26,933,093 (-) · locus `L0178` holding **12** candidate genes
- **Found by** — GEA+GWAS · GEA 1 climate clusters / 1 axes (bio17), best bio15 nlp 7.9, λ_min 2.67 · GWAS 2 gardens (12,27), nlp 8.2, MAC 5
- **Evidence** — 2 independent line(s): GWAS gardens, both scans
- **Blocks** — Chr5_14048
- **Variants** — 3 (1 SV / 2 indel / 0 MNP), max 59 bp, in 5'UTR,promoter; strongest mechanism **F3_UTR**
- **Lead** — Chr5:26933736 (smallindel, 1 bp, promoter), 11 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.07** / med 0.02 · best local r² 0.16
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² 0.80
- **Function (UniProt)** — May act as a substrate-specific adapter of an E3 ubiquitin-protein ligase complex (CUL3-RBX1-BTB) which mediates the ubiquitination and subsequent proteasomal degradation of target proteins. {ECO:0000250}.
- **Keywords** — Alternative splicing;Cytoplasm;Metal-binding;Reference proteome;Ubl conjugation pathway;Zinc;Zinc-finger

![BT4](plots/loci/BT4_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr5:26,933,736 | smallindel | 1 | promoter | F4_promoter | 7.9 | 1 | 2 |
| Chr5:26,932,843 | sv | 59 | 5'UTR | F3_UTR | — | 1 | 2 |
| Chr5:26,932,825 | smallindel | 13 | 5'UTR | F3_UTR | — | 1 | 2 |

---

## GRXS2 (AT5G18600) — **MODERATE**

*recurrent, rare (6 carriers), direct placement (F4_promoter).*

- **Protein** — Monothiol glutaredoxin-S2 (AtGrxS2) (Protein ROXY 10)
- **Themes** — stress
- **Location** — Chr5:6,183,258–6,183,954 (-) · locus `L0151` holding **20** candidate genes
- **Found by** — GWAS · GWAS 6 gardens (11,12,42,55,60,9), nlp 9.1, MAC 6
- **Evidence** — 1 independent line(s): GWAS gardens
- **Blocks** — Chr5_3057
- **Variants** — 1 (0 SV / 1 indel / 0 MNP), max 1 bp, in promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr5:6184898 (smallindel, 1 bp, promoter), 6 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.00** / med 0.00 · best local r² 0.30
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² 0.83
- **Function (UniProt)** — May only reduce GSH-thiol disulfides, but not protein disulfides. {ECO:0000305}.
- **Keywords** — 2Fe-2S;Cytoplasm;Iron;Iron-sulfur;Metal-binding;Redox-active center;Reference proteome

![GRXS2](plots/loci/GRXS2_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr5:6,184,898 | smallindel | 1 | promoter | F4_promoter | — | 0 | 6 |

---

## NIP1-1 (AT4G19030) — **WEAK**

*not recurrent, common (90 carriers), indirect placement (nan).*

- **Protein** — Aquaporin NIP1-1 (NOD26-like intrinsic protein 1-1) (AtNIP1;1) (Nodulin-26-like major intrinsic protein 1) (NodLikeMip1) (Protein NLM1)
- **Themes** — stress
- **Location** — Chr4:10,421,525–10,423,487 (-) · locus `L0127` holding **16** candidate genes
- **Found by** — GEA · GEA 2 climate clusters / 6 axes (bio1,bio17), best bio1 nlp 12.6, λ_min 1.62
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr4_5223
- **Variants** — 1 (0 SV / 0 indel / 1 MNP), max 0 bp, in intron; strongest mechanism **nan**
- **Lead** — Chr4:10422792 (mnp, 0 bp, intron), 90 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.98** / med 0.57 · best local r² 0.98
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² not testable
- **Function (UniProt)** — Water channel probably required to promote glycerol permeability and water transport across cell membranes. {ECO:0000269|PubMed:11007982}.
- **Keywords** — Acetylation;Membrane;Phosphoprotein;Reference proteome;Repeat;Transmembrane;Transmembrane helix;Transport

![NIP1-1](plots/loci/NIP1-1_combined.png)

---

## GSH1 (AT4G23100) — **WEAK**

*not recurrent, rare (11 carriers), direct placement (F3_UTR).*

- **Protein** — Glutamate--cysteine ligase, chloroplastic (EC 6.3.2.2) (Gamma-ECS) (GCS) (Gamma-glutamylcysteine synthetase) (Protein ROOT MERISTEMLESS 1) (AtGCL) (Protein cadmium-sensitive 2) (Protein phytoalexin-deficient 2)
- **Themes** — climate_temp,flowering,stress
- **Location** — Chr4:12,102,135–12,107,460 (-) · locus `L0128` holding **19** candidate genes
- **Found by** — GEA · GEA 2 climate clusters / 3 axes (bio1,bio17), best pc1 nlp 9.4, λ_min 1.85
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr4_6271
- **Variants** — 1 (0 SV / 1 indel / 0 MNP), max 1 bp, in 3'UTR; strongest mechanism **F3_UTR**
- **Lead** — Chr4:12103235 (smallindel, 1 bp, 3'UTR), 11 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.45** / med 0.01 · best local r² 0.46
- **SNPs** — a Bonferroni SNP within 2 kb · tagging r² 0.67
- **Function (UniProt)** — Seems to play an important role in controlling the expression of resistance responses like the regulation of salicylic acid (SA) and phytoalexin (camalexin) production. Involved in resistance to fungal and bacterial pathogens. Required for the regulation of cell proliferation in root apical meristems through the GSH-dependent developmental pathway. Also participates in the detoxification process, the antioxidant response and is essential for embryo development and proper seed maturation. {ECO:0000269|PubMed:10634910, ECO:0000269|PubMed:11402187, ECO:0000269|PubMed:11722772, ECO:0000269|PubMed:
- **Keywords** — Alternative splicing;ATP-binding;Chloroplast;Disulfide bond;Glutathione biosynthesis;Ligase;Nucleotide-binding;Plant defense;Plastid;Reference proteome;Transit peptide

![GSH1](plots/loci/GSH1_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr4:12,103,235 | smallindel | 1 | 3'UTR | F3_UTR | 9.4 | 2 | 0 |

---

## GLR1.3 (AT5G48410) — **WEAK**

*not recurrent, rare (5 carriers), direct placement (F3_UTR).*

- **Protein** — Glutamate receptor 1.3 (Ligand-gated ion channel 1.3)
- **Themes** — circadian_light
- **Location** — Chr5:19,620,267–19,623,425 (+) · locus `L0169` holding **7** candidate genes
- **Found by** — GWAS · GWAS 2 gardens (55,60), nlp 9.0, MAC 5
- **Evidence** — 1 independent line(s): GWAS gardens
- **Blocks** — Chr5_9812,Chr5_9818  ⚠ block attribution would name a different gene
- **Variants** — 5 (0 SV / 3 indel / 2 MNP), max 16 bp, in 3'UTR,intergenic,promoter; strongest mechanism **F3_UTR**
- **Lead** — Chr5:19623282 (mnp, 0 bp, 3'UTR), 5 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.25** / med 0.01 · best local r² 0.32
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² 1.00
- **Function (UniProt)** — Glutamate-gated receptor that probably acts as a non-selective cation channel. May be involved in light-signal transduction and calcium homeostasis via the regulation of calcium influx into cells.
- **Keywords** — Glycoprotein;Ion channel;Ion transport;Ligand-gated ion channel;Membrane;Receptor;Reference proteome;Signal;Transmembrane;Transmembrane helix;Transport

![GLR1.3](plots/loci/GLR1.3_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr5:19,623,291 | smallindel | 1 | 3'UTR | F3_UTR | — | 0 | 2 |
| Chr5:19,619,649 | smallindel | 16 | promoter | F4_promoter | — | 0 | 2 |
| Chr5:19,623,444 | smallindel | 5 | intergenic | F5_downstream | — | 0 | 2 |

---

## APR2 (AT1G62180) — **WEAK**

*not recurrent, rare (7 carriers), direct placement (F4_promoter).*

- **Protein** — 5'-adenylylsulfate reductase 2, chloroplastic (EC 1.8.4.9) (3'-phosphoadenosine-5'-phosphosulfate reductase homolog 43) (PAPS reductase homolog 43) (Prh-43) (Adenosine 5'-phosphosulfate 5'-adenylylsulfate sulfotransferase 2) (APS sulfotransferase 2) (Thioredoxin-independent APS reductase 2)
- **Themes** — stress
- **Location** — Chr1:22,975,530–22,977,625 (-) · locus `L0036` holding **10** candidate genes
- **Found by** — GEA · GEA 2 climate clusters / 4 axes (bio1,pc3), best pc3 nlp 11.8, λ_min 1.99
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr1_10997
- **Variants** — 1 (0 SV / 1 indel / 0 MNP), max 1 bp, in promoter; strongest mechanism **F4_promoter**
- **Lead** — Chr1:22978065 (smallindel, 1 bp, promoter), 7 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.12** / med 0.01 · best local r² 0.27
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² 0.35
- **Function (UniProt)** — Reduces sulfate for Cys biosynthesis. Substrate preference is adenosine-5'-phosphosulfate (APS) >> 3'-phosphoadenosine-5'-phosphosulfate (PAPS). Uses glutathione or DTT as source of protons. {ECO:0000269|PubMed:11553635}.
- **Keywords** — Acetylation;Alternative splicing;Amino-acid biosynthesis;Chloroplast;Cysteine biosynthesis;Disulfide bond;Oxidoreductase;Plastid;Redox-active center;Reference proteome;Stress response;Transit peptide

![APR2](plots/loci/APR2_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr1:22,978,065 | smallindel | 1 | promoter | F4_promoter | 11.8 | 2 | 0 |

---

## UPM1 (AT5G40850) — **FRAGILE**

*not recurrent, rare (11 carriers), indirect placement (F6_intron).*

- **Protein** — S-adenosyl-L-methionine-dependent uroporphyrinogen III methyltransferase, chloroplastic (Urophorphyrin III methylase) (EC 2.1.1.107) (Urophorphyrin methylase 1) (AtUPM1)
- **Themes** — circadian_light,climate_temp
- **Location** — Chr5:16,366,353–16,369,099 (+) · locus `L0166` holding **14** candidate genes
- **Found by** — GEA · GEA 2 climate clusters / 2 axes (bio1,bio17), best pc1 nlp 8.2, λ_min 1.85
- **Evidence** — 1 independent line(s): GEA recurrence
- **Blocks** — Chr5_7243
- **Variants** — 1 (0 SV / 1 indel / 0 MNP), max 9 bp, in intron; strongest mechanism **F6_intron**
- **Lead** — Chr5:16367086 (smallindel, 9 bp, intron), 11 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.60** / med 0.01 · best local r² 0.75
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² 0.81
- **Function (UniProt)** — Essential protein required for siroheme biosynthesis (PubMed:20592802). Catalyzes the two successive C-2 and C-7 methylation reactions involved in the conversion of uroporphyrinogen III to precorrin-2 via the intermediate formation of precorrin-1 (PubMed:9006913). It is a step in the biosynthesis of siroheme (PubMed:9006913). Promotes nitrogen and sulfur assimilation as well as photosynthesis efficiency by triggering chlorophyll, nitrite reductase (NiR) and sulfite reductase (SiR) biosynthesis (PubMed:29472934). {ECO:0000269|PubMed:20592802, ECO:0000269|PubMed:29472934, ECO:0000269|PubMed:9006
- **Keywords** — Chloroplast;Methyltransferase;Plastid;Porphyrin biosynthesis;Reference proteome;S-adenosyl-L-methionine;Transferase;Transit peptide

![UPM1](plots/loci/UPM1_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr5:16,367,086 | smallindel | 9 | intron | F6_intron | 8.2 | 2 | 0 |

---

## PRT6 (AT5G02310) — **FRAGILE**

*not recurrent, rare (5 carriers), indirect placement (F6_intron).*

- **Protein** — E3 ubiquitin-protein ligase (EC 2.3.2.27)
- **Themes** — climate_temp,stress
- **Location** — Chr5:474,279–482,783 (+) · locus `L0140` holding **9** candidate genes
- **Found by** — GWAS · GWAS 2 gardens (12,27), nlp 10.1, MAC 5
- **Evidence** — 1 independent line(s): GWAS gardens
- **Blocks** — Chr5_308
- **Variants** — 1 (0 SV / 1 indel / 0 MNP), max 1 bp, in intron; strongest mechanism **F6_intron**
- **Lead** — Chr5:482054 (smallindel, 1 bp, intron), 5 founder carriers
- **LD (annotation, not the verdict)** — r²(lead↔gene) max **0.12** / med 0.01 · best local r² 0.25
- **SNPs** — no Bonferroni SNP within 2 kb (SNP-unique) · tagging r² 1.00
- **Function (UniProt)** — Ubiquitin ligase protein which is a component of the N-end rule pathway. Recognizes and binds to proteins bearing specific N-terminal residues that are destabilizing according to the N-end rule, leading to their ubiquitination and subsequent degradation. {ECO:0000256|RuleBase:RU366018}.
- **Keywords** — Metal-binding;Proteomics identification;Reference proteome;Transferase;Ubl conjugation pathway;Zinc;Zinc-finger

![PRT6](plots/loci/PRT6_combined.png)

| variant | class | size | region | mechanism | GEA nlp | clusters | gardens |
|---|---|---|---|---|---|---|---|
| Chr5:482,054 | smallindel | 1 | intron | F6_intron | — | 0 | 2 |

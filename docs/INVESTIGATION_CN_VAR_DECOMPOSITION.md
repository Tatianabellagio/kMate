# Investigation report: cn_var decomposition bug + open-loop validation framework

**Date:** 2026-05-19
**Scope:** Audit of cn_var construction across v3 → v3qc → v3qc_v2 → v3qc_v3 panels; survey of decomposition tools; review of validation methodology; consolidated path forward.
**Status:** Investigation complete. Implementation pending three open decisions (Section 11).

---

## 1. TL;DR

1. **There is a systematic bug in `cn_var` construction** that under-counts carriers at SNPs co-located with multi-allelic indels (the "cell [d]" failure mode). Caused by `bcftools norm -m -any` scattering carriers across decomposed biallelics at shifted positions.
2. **The bug went undetected through v3, v3qc, v3qc_v2, v3qc_v3** because the simulation validation is **closed-loop**: `compute_recomb_truth.py` computes truth-AF by indexing into the same cn_var matrix that cactus_em projects through. Any cn_var bias cancels out in the truth-vs-estimate comparison.
3. **Three "fixes" investigated and rejected in this session:**
   - Coord-aware coalescence patch: works for under-counting cases (10421645, 13843898) but introduces over-counting at indel-shifted bases (5870018-type).
   - `bcftools norm --atomize`: does positional byte-offset decomposition; produces ~35× false positives at indel-adjacent SNPs (verified against long-read assemblies).
   - `vcfwave`: silently deduplicates atomic records at converging-ALT sites and **drops per-sample GTs from all but one source ALT** (confirmed in `src/vcfwave.cpp`). 80 expected carriers → 1 observed at one PG bubble.
4. **The field solves this with graph-native or symbolic-ID architectures**, not pairwise alignment. HPRC, pantree (Aug 2025 preprint), and `pangenome/resolve-nested-genotypes` all bypass the alignment-based carrier-loss class.
5. **No published method does pool-seq + pangenome + multi-allelic AF estimation with open-loop validation.** Our cactus_em is a novel contribution; the field gap is named explicitly in Lehmann et al. 2025 *MBE* review.
6. **Recommended architecture:** `vg deconstruct -a -e` → `pangenome/resolve-nested-genotypes` (Arch 1). Field-validated components, graph-native both at catalog and propagation, escapes the bug class permanently.
7. **Recommended validation:** 3-tier (open-loop estimators + FASTA-truth + LOO concordance). Closed-loop comparison against `compute_recomb_truth.py` forbidden going forward.

---

## 2. The bug — mechanism and evidence

### 2.1 What's wrong with cn_var

Cactus pangenome multi-allelic snarls are decomposed in Phase A by `bcftools norm -f REF -m -any`. This:
1. Splits each multi-allelic record into one biallelic per ALT.
2. Per-pair, strips shared prefix/suffix between REF and ALT_j.
3. Anchors each biallelic at the position after prefix-stripping.

**Failure mode 1 — position scatter (under-counting).** When a 530bp multi-allelic snarl with 72 ALTs decomposes, biallelics scatter across many positions because each ALT's prefix-stripping produces a different anchor. At Chr1:10421232 (530bp × 72 ALTs), `bcftools norm` produces 72 biallelic records spanning positions 10421232–10421645+. A SNP at TAIR10:10421645 (T→C) only retains carriers of the one ALT path that happened to anchor there post-stripping. Other ALT paths with C at coord 10421645 are stored at different biallelic positions → invisible to a same-position carrier query.

Empirical: cn_var AC=2 of 231 vs xwu/hapFIRE truth AC=229. ~99% carrier loss at this site.

**Failure mode 2 — byte-offset misalignment (over-counting via patches/atomize).** When you try to "fix" by reading `ALT[byte_offset]` to attribute carriers across spanning records, the byte-offset only maps to TAIR10 coord (pos + offset) IF the ALT has no internal indels. For ALTs with internal indels, the bytes past the indel are shifted.

Empirical at Chr1:5869846 (260bp × 12 ALTs, with 8bp insertion in 9 of 12 ALTs):
- Byte offset 172 = `A` in 9 of 12 ALTs (naive read)
- But proper WFA pairwise alignment (vcfwave + edlib + parasail all agree): only 1–2 of 12 ALTs actually encode A at TAIR10:5870018
- Direct long-read assembly alignment (5 founders): only 100764 + 100861 carry A; matches WFA
- xwu truth AC≈5/231 (consistent)
- `bcftools norm --atomize` and naive coord-aware patches: claim 70–117 carriers (35–60× over-count)

### 2.2 Carrier loss IS in cn_var, not elsewhere

- **`cn_full` is fine.** Built from per-founder consensus FASTAs via `bcftools consensus`. Consensus correctly applies each biallelic substitution at its anchored position, so founder FASTAs end up with the correct base at each coord regardless of which scattered biallelic encoded the substitution. k-mer extraction from coord-correct FASTAs gives coord-correct cn_full.
- **EM math is fine.** Operating on coord-correct cn_full, the EM recovers founder weights `h` faithfully. Verified by single-founder pure-pool tests showing correct h-vector recovery.
- **AF projection is the broken step.** `h @ cn_var / h @ cn_var_called` projects through buggy cn_var → systematically biased per-record AF.

### 2.3 Where in the genome this hurts

Per the `coord_aware_full.tsv` analysis (Section 4): 2.81% of Chr1 SNPs receive extras when coord-aware coalescence is attempted. Stratified by F_MISSING × mixed-bubble:
- Cell [a] (low F_MISSING, pure SNP): 0.33% outlier rate — bug barely visible.
- Cell [b] (low F_MISSING, mixed bubble): 11.7% outlier rate — bug prominent.
- Cell [c] (high F_MISSING, pure SNP): 4.8%.
- Cell [d] (high F_MISSING, mixed bubble): **38.6%** — bug dominant.

Direction: 99% of cell [d] outliers are **under-prediction** — classic carrier-scatter signature. 83% of all outliers genome-wide are under-predictions.

---

## 3. Why we missed it — the closed-loop validation trap

### 3.1 What `compute_recomb_truth.py` actually does

The simulation truth computation (`sims/visor_freqk/scripts/compute_recomb_truth.py`) computes truth-AF by literally reading the same `cn_var` matrix that cactus_em uses for AF projection:

```python
cn = load_npz(args.cn_var).tocsr()  # SAME cn_var as cactus_em
...
for ind_idx, ind_id in enumerate(inds):
    founder_per_rec = assign_founder_per_record(rec_chrom, rec_pos, segs_by_chrom)
    for f in set(founder_per_rec):
        row_dense = np.asarray(cn[row_idx, :].todense()).flatten()
        ind_contrib[mask] = row_dense[mask]
    truth_af += w * ind_contrib
```

For non-recombinant pure-pool sims, this is exactly `truth_af = w @ cn_var`. For recomb sims, it's the same with per-position founder lookup. Both truth and cactus_em estimate are linear functions of cn_var. **Any cn_var bug is algebraically invisible to the comparison.**

### 3.2 hapFIRE on sims is also closed-loop

`FINAL_RESULTS_cov10_v3.ipynb` runs three "hapFIRE" variants for cross-comparison:
- `hapfire_v2panel_v3proj`: hapFIRE h-vectors from v2 panel projected through cn_var_v3
- `hapfire_v3panel`: hapFIRE on a SNP-only carrier panel built FROM cn_var_v3 carriers
- All also go through cn_var_v3 (either at projection step or input panel construction)

The 7–20 pp gap between cactus_em (~0.99 R²) and hapfire-v3panel (~0.89 R²) measures **h-vector quality difference**, not cn_var correctness. Both estimators share the cn_var bias.

### 3.3 Symptom was visible but mis-attributed

The notebook's own markdown observed:
- v3 projection costs ~10–23 R² pp vs v2 projection
- v3 estimates have slope ≈ 0.85 (systematic under-prediction)

Both are direct symptoms of carrier-scatter under-counting. Attributed in the notebook to "many atomized rows share predictions" and "Beagle was masking it" — partial truths that don't reach the actual mechanism.

### 3.4 First crack: real-data SEEDMIX_S1 vs xwu hapFIRE

The bug became visible only when comparing cactus_em against an **independent** estimator that doesn't go through our cn_var — xwu's hapFIRE on real GrENE-Net pool-seq data. That comparison surfaced the 0.88% outlier rate and the directional under-prediction at SNPs-in-bubbles that all closed-loop sim eval had hidden.

---

## 4. What we tried this session

| attempt | result | conclusion |
|---|---|---|
| Same-pos base-at-pos coalescence | Closed ~45% of cell [d] gap | Worked only where bcftools-norm didn't shift positions |
| Coord-aware coalescence patch (across all spanning records) | Closed ~68% of cell [d] gap, MAE 0.175 → 0.056 | Works for under-counting but introduces saturation/over-shoot at indel-adjacent sites |
| Verified the coord-aware patch via `coord_aware_full.py` on 516k Chr1 SNPs | Genome-wide MAE 0.0150 → 0.0134, outliers 0.88% → 0.45% | Real improvement but not the right fix |
| `bcftools norm --atomize` | At 10421645 and 13843898: appears to fix. At 5870018: introduces 35× false positives. | Byte-offset alignment, not pairwise. Wrong tool. |
| `vcfwave` on the over-counting case 5870018 | 117 → 2 carriers (matches long-read truth) | Correctly resolves over-counting |
| `vcfwave` on the under-counting cases 10421645, 13843898 | Recovers high AC (130/135, 109/109) | Correctly resolves under-counting |
| GT preservation audit on cactus haploid bubble | One-to-one edlib match for founder 100130 | vcfwave preserves haploid GTs correctly on this bubble |
| GT preservation audit on PG diploid bubble (K≥2 hom) | All K≥2 GTs cleanly atomized; missing GTs propagate | vcfwave handles diploid K≥2 correctly when no convergence |
| **GT preservation on PG bubble with converging ALTs (Chr1:87626)** | **80+ expected carriers → 1 observed in GT field** | **vcfwave silently drops carriers at converging atomic sites** |
| Source-code verification of vcfwave bug | `src/vcfwave.cpp` sums INFO/AC/AF on collision but doesn't re-union per-sample `RecGenotypes` | Bug is by design, not misconfiguration. Undocumented upstream. |

---

## 5. Field survey — what others do

### 5.1 Decomposition tools

| Tool | Behavior on convergence case | Use? |
|---|---|---|
| `bcftools norm -m -any` (no --atomize) | Scatters carriers across positions | NO |
| `bcftools norm --atomize` | Byte-offset, breaks on internal indels (bcftools #2239 OPEN) | NO |
| `vcfwave` (vcflib) | Silent per-sample GT loss on dedup (src confirmed) | Catalog-only OK; per-sample propagation FORBIDDEN |
| `vcfallelicprimitives` (vcflib, legacy) | Same as vcfwave | Superseded |
| `vt decompose_blocksub` | Biallelic-only; can't handle multi-allelic | N/A |
| `vg deconstruct -a -e` | Emits LV/PS/AT tags per snarl; graph-native | YES (with downstream) |
| `pangenome/resolve-nested-genotypes` | Path-based AT-substring propagation parent→nested | YES (paired with above) |
| `eblerjana convert-to-biallelic.py` | Symbolic-ID-based propagation; field-tested by HPRC | YES (alt) |
| `pantree` (Aug 2025) | Reference-tree variant definition; fixes vcfwave/graph inconsistency | EXPERIMENTAL |

### 5.2 HPRC's actual practice

`prepare-vcf-MC` → annotated multi-allelic VCF (per-ALT ID lists naming nested atomic variants) + biallelic catalog. PanGenie genotypes on the multi-allelic. `convert-to-biallelic.py` propagates per-sample GTs to biallelic catalog by **symbolic ID matching** — no alignment, correct convergence handling.

Used in HPRC v1.0 + v1.1, Liao et al. 2023 *Nature*. The atomic-variant catalog itself is still built via vcfwave, so it inherits a ~0.8% under-coverage gap (pantree paper quantifies as 239,375 SNPs missing from HPRC).

### 5.3 Pool-seq pangenome AF estimation

**No published method does what we're doing.** Every prior pool-seq tool (HARP, hapFIRE, HAF-pipe, PoPoolation2, PoolSNP, SNAPE-pool, grenedalf, Pool-hmm) operates on biallelic SNP VCF against a linear reference. PanGenie/vg-giraffe/Locityper/Varigraph are graph-native but per-individual genotypers.

GrENE-Net's published methods (Czech et al. 2022 bioRxiv; Wu, Bellagio et al. 2026 *Science*) use HARP via hapFIRE against the 1001G biallelic SNP catalog on linear TAIR10 — no SVs, no indels >1bp, no graph.

The gap our `cactus_em` fills is explicitly named in Lehmann et al. 2025 *MBE* review as an open methodological need.

### 5.4 Closed-loop validation is field-wide

Most pangenome genotyping validation is closed-loop (PanGenie vs DeepVariant on same individual reads; Long Lab AFs cross-checked against same-BAM REF/ALT counts). The HPRC + pangenome-aware DeepVariant (Cook et al. 2025) is one of the few that uses an independent truth (T2T-Q100). **No published pool-seq pangenome AF paper uses FASTA-based open-loop truth.**

---

## 6. Recommended architecture

### Arch 1 (PRIMARY): `vg deconstruct -a -e` + `pangenome/resolve-nested-genotypes`

Graph-native at both decomposition and per-sample propagation. No alignment-based carrier loss possible.

```
cactus pangenome graph (.gbz)
       │
       │  vg deconstruct -a -e
       │  (emit one record per snarl, all levels;
       │   tag with LV, PS, AT)
       ▼
multi-allelic VCF with snarl-tagged records
       │
       │  resolve-nested-genotypes
       │  (substring-match AT strings;
       │   propagate parent-snarl GTs → nested-snarl GTs)
       ▼
biallelic atomic VCF with correct per-sample GTs
       │
       │  build_cn_var.py (unchanged)
       ▼
cn_var, cn_var_called
```

**Pros:** graph-native; field-adjacent (used by HPRC for catalog generation); handles convergence by construction; preserves bubble metadata for downstream queries.

**Cons:** `resolve-nested-genotypes` documentation thin on convergence case (must verify on our 3 spot-checks); Rust toolchain dependency; vg-deconstruct -a output may be large (15–40M records genome-wide).

**Estimated cost:** 2–3 days.

### Arch 3 (FALLBACK): HPRC `prepare-vcf-MC` + `convert-to-biallelic.py`

Symbolic-ID-based propagation, field-tested by HPRC.

**Pros:** smallest diff to current pipeline; production-tested in Liao et al. 2023.
**Cons:** atomic-variant catalog still uses vcfwave internally → ~0.8% atomic variants missing.
**Estimated cost:** 1–2 days.

### Arch 2 (PARALLEL COMPARISON): `pantree` (Salehi Nowbandegani et al. bioRxiv 2025-08)

Most principled treatment (reference-tree formalization, Heng Li's group). Recovers 3.5M variants vcfwave misses on HPRC graph.

**Pros:** mathematically rigorous; explicitly fixes the vcfwave/graph inconsistency.
**Cons:** brand new (zero ecosystem adoption); variant-edge output requires our own coord-projection layer.
**Use:** parallel baseline for comparison, not primary path. +2 days.

### Arch 4 (NOT NOW): multi-allelic-native cn_var

Skip decomposition entirely, store cn_var as 3D ragged sparse with per-(record, alt) carriers. Architecturally cleanest endpoint. Effectively what `pantree` formalizes.

**Cons:** largest implementation lift; downstream rewrite (recipe, comparison, output schema).
**Use:** longer-horizon ideal; revisit after Arch 1 ships.

---

## 7. Validation framework — non-negotiable

The closed-loop trap took 6+ months to surface. The validation methodology adopted going forward IS itself a contribution to the field.

### Tier 1 — MANDATORY before declaring rebuild done

**A — Independent estimator on simulated pool BAMs.** Three estimators in parallel:
- `cactus_em` (our method)
- `grenedalf frequency` (Czech 2024, by the GrENE-Net consortium; model-free, BAM-direct)
- `bcftools mpileup -a FORMAT/AD` (rawest per-site allele depth)

Join on (chrom, pos, ref, alt). Pairwise MAE. Stratify by `is_in_multiallelic_bubble`. **Bug canary**: cactus_em should match grenedalf and mpileup at non-bubble sites AND at bubble-co-located sites if the rebuild worked.

**B — FASTA-lookup truth.** Per-coord truth as `Σ_i w_i × [base_at(founder_i.fa, coord) == ALT]`. Uses `pyfaidx`. Never touches cn_var. Algebraically independent. Canonical truth from here on.

### Tier 2 — Strongly recommended (+2 days)

**C — Stratified leave-one-out.** 10 founders × held-out-panel rebuild × `vcfeval` + `Truvari bench/phab`. HPRC-standard methodology. Gate: median wGC ≥ 0.95 (SNP), ≥ 0.85 (SV).

### Tier 3 — Quick polish (+1 day)

- `vcfdist` (Dunn 2023) on cactus_em vs hapFIRE: representation-independent variant comparison.
- `Truvari phab` on the cactus side specifically for SVs.

---

## 8. Gates for the rebuild

| gate | metric | target |
|---|---|---|
| G1a | Overall outlier rate \|d\|>0.10 vs hapFIRE | **<0.20%** (current 0.88%) |
| G1b | Cell [d] outlier rate | **<5%** (current 38.6%) |
| G1c | Cell [b] outlier rate | **<1.5%** (current 11.7%) |
| G1d | Cell [a] outlier rate | **≤0.3%** no regression |
| G2 | Per-cell pos/neg outlier split | **within 60/40 each cell** |
| G3 | Outlier concentration in known-issue regions | **≥60% in centromere + knob + flagged-founder positions** |
| G4 | cactus_em ↔ grenedalf MAE | **<0.02** at non-bubble sites; <0.05 at bubble sites |
| G5 | cactus_em ↔ FASTA-lookup truth MAE | **<0.02** (open-loop primary) |
| G6 | LOO median wGC (SNPs) | **≥0.95** |
| G7 | LOO median wGC (SVs) | **≥0.85** |

---

## 9. Sequenced implementation plan

### Week 1 — Pre-flight verification

- **Day 1:** Install Rust + build `resolve-nested-genotypes`. Run on a small Chr1 test region. Spot-check 5870018, 10421645, 13843898 vs known truth.
  - Expected: 5870018 → ~2 carriers; 10421645 → ~129; 13843898 → ~109.
  - If fails: fall back to Arch 3.
- **Day 2:** Run vg deconstruct -a on full Chr1. Measure runtime and record count. Extrapolate genome-wide.
- **Day 3:** Install pantree from source. Run on the same 3 spot-checks. Document comparison.
- **Day 4:** File vcfwave GT-loss bug upstream (vcflib issue) with our 80→1 example. Helps the community.

### Week 2 — Scaffold + Chr1 pilot

- **Day 5–6:** Replace `bcftools norm -m -any` with `vg deconstruct -a -e | resolve-nested-genotypes` in phase A scripts. Update build_cn_var.py if needed (likely no change). Haploidization, fill-tags, merge, AC=0 cleanup stay as-is.
- **Day 7:** Full Chr1 pilot. SEEDMIX_S1 end-to-end. Spot-check 3 known sites.
- **Day 8–9:** Tier 1 validation (grenedalf, mpileup, FASTA-truth) on Chr1.
- **Day 10:** Gate check (G1–G5). If pass → scale. If fail → diagnose.

### Week 3 — Genome scale + final validation

- **Day 11–13:** Genome-wide rebuild (all 5 chroms). cn_var + cn_full + FASTAs.
- **Day 14:** Tier 2 LOO validation (G6–G7).
- **Day 15:** Tier 3 + write up methodology for repo + paper.

---

## 10. Risk register

| severity | risk | mitigation |
|---|---|---|
| HIGH | `resolve-nested-genotypes` convergence case not explicitly documented | Day 1 spot-check is gating; fall back to Arch 3 if it fails |
| HIGH | Wrong architecture choice locks in another months-long bug cycle | Run Arch 1 + pantree in parallel on validation cases before committing |
| MEDIUM | vg deconstruct -a record count >50M, breaks cn_var memory | Bench on Chr1 first; chunk by region if needed |
| MEDIUM | PG side integration (PanGenie outputs not in AT-tagged form) | May need HPRC `convert-to-biallelic.py` for PG side; vg approach for cactus side |
| MEDIUM | FASTA-truth pyfaidx performance at 5M coords × 231 founders | Bench; likely sub-hour but verify |
| MEDIUM | cn_full rebuild required (founder FASTAs change with new VCF) | Document explicitly; estimated 6–8 h genome-wide |
| LOW | Methodological paper opportunity if framing right | Document validation framework carefully during rebuild |

---

## 11. Open decisions for the user

1. **Arch 1 commit vs Arch 1+Arch 3 parallel.** If Day 1 spot-check of `resolve-nested-genotypes` passes the convergence case, commit to Arch 1. If it fails, fall back to Arch 3 (+1 day). Want to also run Arch 3 in parallel as belt-and-suspenders?

2. **`pantree` parallel baseline (+2 days).** Most principled tool but brand new. Running alongside gives methodological-paper-grade comparison. Worth it or skip?

3. **Validation rigor.** Tier 1 mandatory. Tier 2 (LOO, +2 days) strongly recommended. Tier 3 (vcfdist + Truvari phab, +1 day) for full benchmark suite. How much?

---

## 12. Framing shift: novelty and contribution

The pool-seq pangenome literature search found **no published method does what cactus_em does**. The closest published methods (HARP/hapFIRE for pool-seq; PanGenie for graph genotyping) are each limited to half of our problem. The Lehmann et al. 2025 *MBE* review explicitly names the methodological gap.

Implications:
- `cactus_em` is a methodological contribution beyond the GrENE-Net Science paper, not just pipeline plumbing.
- The closed-loop validation pattern is field-wide; demonstrating an open-loop framework on a real bug class is publishable.
- The vcfwave per-sample GT loss is undocumented; filing benefits the community.
- A methods paper writes itself out of the rebuild work we're already doing.

---

## 13. Memory updates needed

After rebuild ships, update / create:
- [ ] `project_cn_var_decomposition_fix.md` — record what was fixed and how
- [ ] `project_validation_framework_open_loop.md` — document grenedalf + FASTA-truth + LOO methodology
- [ ] `feedback_compute_recomb_truth_is_closed_loop.md` — warn future-Claude not to trust truth_vs_estimate R² as panel correctness signal
- [ ] `feedback_vcfwave_dedup_loses_carriers.md` — record the vcfwave GT-loss bug for any future decomposition decisions
- [ ] `feedback_bcftools_atomize_byte_offset_bug.md` — document why --atomize isn't a fix
- [ ] Update `MEMORY.md` index

---

## 14. Sources

### Decomposition / architecture
- [vcfwave.cpp source confirming GT dedup bug](https://github.com/vcflib/vcflib/blob/master/src/vcfwave.cpp)
- [bcftools #2239 — atomize zygosity loss](https://github.com/samtools/bcftools/issues/2239)
- [bcftools #1689 — atomize unphases](https://github.com/samtools/bcftools/issues/1689)
- [bcftools #1668 — atomize REF check breaks](https://github.com/samtools/bcftools/issues/1668)
- [vcflib #360 — vcfwave sample-AC missing](https://github.com/vcflib/vcflib/issues/360)
- [vcflib #428 — vcfwave annotation scope](https://github.com/vcflib/vcflib/issues/428)
- [vg deconstruct wiki](https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct)
- [pangenome/resolve-nested-genotypes](https://github.com/pangenome/resolve-nested-genotypes)
- [resolve-nested-genotypes issue #2](https://github.com/pangenome/resolve-nested-genotypes/issues/2)
- [pantree paper, bioRxiv Aug 2025](https://www.biorxiv.org/content/10.1101/2025.08.04.668502v1)
- [pantree code](https://github.com/oclb/pantree)
- [eblerjana/genotyping-pipelines/prepare-vcf-MC](https://github.com/eblerjana/genotyping-pipelines/tree/main/prepare-vcf-MC)
- [convert-to-biallelic.py](https://github.com/eblerjana/pangenie/blob/master/pipelines/run-from-callset/scripts/convert-to-biallelic.py)

### Pool-seq + pangenome literature
- [GrENE-Net pilot paper, Czech et al. 2022 bioRxiv](https://www.biorxiv.org/content/10.1101/2022.02.02.477408v2)
- [Wu, Bellagio et al. 2026 Science](https://www.science.org/doi/10.1126/science.adz0777)
- [grenedalf, Czech et al. 2024 Bioinformatics](https://academic.oup.com/bioinformatics/article/40/8/btae508/7741639)
- [hapFIRE GitHub](https://github.com/moiexpositoalonsolab/hapfire)
- [HARP, Kessner et al. 2013 MBE](https://academic.oup.com/mbe/article/30/5/1145/993804)
- [HAF-pipe, Tilk et al. 2019 G3](https://pmc.ncbi.nlm.nih.gov/articles/PMC6893198/)
- [PanGenie, Ebler et al. 2022 Nat Genet](https://www.nature.com/articles/s41588-022-01043-w)
- [HPRC, Liao et al. 2023 Nature](https://www.nature.com/articles/s41586-023-05896-x)
- [Minigraph-Cactus, Hickey et al. 2023 Nat Biotech](https://pmc.ncbi.nlm.nih.gov/articles/PMC10638906/)
- [Bridging pangenomics + popgen review, Lehmann et al. 2025 MBE](https://academic.oup.com/mbe/article/42/3/msaf047/8052716)
- [Long Lab "Illusion of Polygenicity", Genetics 2026](https://academic.oup.com/genetics/advance-article/doi/10.1093/genetics/iyag068/8514529)
- [Pangenome-aware DeepVariant, Cook et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12157594/)

### Validation
- [Truvari, Genome Biology 2022](https://link.springer.com/article/10.1186/s13059-022-02840-6)
- [vcfdist, Nat Commun 2023](https://www.nature.com/articles/s41467-023-43876-x)
- [RTG vcfeval](https://github.com/RealTimeGenomics/rtg-tools)
- [hap.py](https://github.com/Illumina/hap.py)
- [pysamstats](https://github.com/alimanfoo/pysamstats)
- [precisionFDA Truth Challenge V2, Olson et al. Cell Genomics 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9205427/)

---

## 15. Local artifacts produced this session

- `scratch/f2_honest/coord_aware_full.tsv` — 516,712 Chr1 SNPs with original + coord-aware AFs (the patch test that revealed the bigger architectural issue)
- `scratch/f2_honest/coord_aware_full.rerun.log` — verification of deterministic reproduction
- `/tmp/vcfwave_uc/` — vcfwave vs --atomize comparisons on under-counting cases
- `/tmp/vcfwave_verify/` — direct vcfwave behavior verification (the 69-record cactus vs 1-record PG comparison that surfaced the dedup bug)
- `/tmp/vcfwave_12het/` — Chr1:87626 GT propagation audit (PG bubble where vcfwave dropped 80 carriers)
- `/tmp/vcfwave_gt_audit_*/` — agent-produced GT preservation audits
- `/tmp/bubble_*.tsv` — bubble structure analyses (offset-by-offset alignment of all ALTs)

All artifacts in `/tmp/` should be considered ephemeral; the persistent record is this document.

# Background — why this project exists, and why the approach is what it is

This file captures the reasoning that led to the current method. `HANDOFF.md` is the current-state pointer + next-steps queue. Read this if you need the *why*.

---

## The problem

Estimate **structural variant (SV) and SNP allele frequencies** in evolved Pool-seq populations from the **GrENE-Net** experiment.

- 2,415 evolved Pool-seq libraries (`MLFH*`), 8 SEEDMIX founder-pool replicates.
- Founder panel: 231 *A. thaliana* ecotypes, of which **80 have long-read assemblies** in the cactus pangenome; the remaining **151 are PanGenie-genotyped** from public 1001 Genomes short reads.
- Existing SV frequency method: **freqk** (k-mer-based per-sample AF). Suffers at low coverage and for SVs with non-unique flanking k-mers.

Goal: a more accurate per-sample SV allele-frequency table that slots into the downstream GEA pipeline.

---

## The method, in one line

**`kMate`**: per-sample k-mer Poisson EM on the 231-founder simplex, projected through a founder × variant matrix (`cn_var_231_arch3`):

```
c_k ~ Poisson( λ · h^T · cn_full[:, k] )   for each panel k-mer k
AF_r = h^T · cn_var[:, r]                  for each VCF record r
```

Joint EM across all ~80M panel k-mers buys identifiability for the 231-vector `h`, then a linear projection through `cn_var` gives per-record SNP + INS + DEL + SV alt frequencies in one pass. See `ALGORITHM.md` (prose) and `CACTUS_EM_MATH.md` (math).

The current kMate pipeline **replaced** an earlier hapFIRE-projection approach (run hapFIRE on SNPs, project its founder-frequency vector onto SVs via a founder × SV genotype matrix). The architectural pattern — recover `h` once, project to many variants — is the same; the inference is now k-mer-EM rather than HARP per-LD-block + CVXPY.

---

## Why a founder-mixture model (and not direct per-SV genotyping)

**hapFIRE** (Wu, Exposito-Alonso lab) and **HARP** (Kessner et al. 2013, PMID 23364324) demonstrated that recovering whole-genome founder frequencies from pool-seq is feasible and that per-SNP frequencies fall out as a linear projection. Conceptually similar to **HAF-pipe** but uses LD-defined haplotype blocks instead of fixed windows.

### Why we cannot just feed SVs into HARP

HARP's likelihood is **per-base**: each read contributes `P(observed base | true base, qual)` at SNP positions, summed over reads. The SNP-table format literally has A/C/G/T columns per haplotype per position. SVs don't fit that model — a deletion or insertion isn't a base substitution; the evidence is split reads / discordant pairs / coverage shifts. Making HARP SV-native would require replacing the entire per-base likelihood — a real rewrite.

### Why kMate instead of "hapFIRE then project"

Three reasons (see `old_docs/METHODS_TRIED.md` for the historical design-space map):

1. **Single pipeline handles SNPs + INS + DEL + SVs natively** — k-mer EM doesn't care what kind of variant the k-mer flanks.
2. **Graph-aware k-mers are richer evidence than per-SNP windows** — pangenome bubbles produce 16-32 k-mers per allele, capturing SV breakpoints directly.
3. **~3× faster end-to-end** than hapFIRE on real pool-seq samples (hapFIRE's Phase 2 HARP loop is serial across LD blocks).

kMate and hapFIRE agree at R²=0.996 on SNPs in real SEEDMIX_S1 — same inverse problem, different inference routes. kMate adds SVs; hapFIRE can't.

---

## Known limitations of this approach

In rough priority for downstream GrENE-Net analysis (current kMate era):

1. **Recombination violates the "every haplotype = a clean founder" assumption.** Evolved GrENE-Net samples have 1–3 generations of recombination → per-block ancestry. The window-mode + global-anchor + HMM-smooth (`★★`) recipe mitigates this; window_200kb still wins on the hardest n50_g3 regime by a narrow margin.
2. **Cactus-vs-PG k-mer asymmetry (+41% cactus h-bias).** Cactus founders have richer k-mer fingerprints (CV 2.6% across founders) than PanGenie-genotyped founders. The EM weighs k-mer evidence and over-credits cactus founders. Hurts the ~0.58% of records where carrier rates differ extremely (PG-specific variants, mostly Chr1q knob Mb 21-23). Mitigation strategies explored in `old_docs/BALANCING_KMERS.md`; production k-mer filter is **undecided** (see `docs/PIPELINE_STATE.md`).
3. **`cn_var` GT disagreement at multi-allelic atomization sites.** `bcftools norm -m -any` produced spurious per-founder carrier calls at SNPs adjacent to INDELs/SVs (`old_docs/OUTLIERS_SUMMARY.md`). **Fix in current production**: arch decomposition (annotate_vcf + convert-to-biallelic) avoids `norm -m -any` entirely. See `docs/INVESTIGATION_CN_VAR_DECOMPOSITION.md`.
4. **231-simplex identifiability.** With many near-identical founder haplotypes, multiple `h` solutions are near-optimal. Per-ecotype recovery is noisy; aggregates (sums over multiple carriers) are accurate.
5. **Beagle imputation hard-rejected for SVs.** LOO testing showed −25 to −30 pp concordance loss on small/medium SVs. Production VCF keeps PanGenie SV calls unimputed. See `panel/pangenie_genotyping/data/merged/README_GOLDEN_STANDARD.md`.
6. **Centromere alignment dead zone (Chr1 Mb 14–17 and pericentromeric regions on all chroms).** Structural; not fixable at EM or `cn_var` level. ~30% of per-SNP outliers live there.

---

## Data provenance

- **Production panel VCF**: `panel/arch3/chr1/merged_231_chr1_final.vcf.gz` — Arch 3 231-founder biallelic haploid panel (Chr1; graph-annotated, symbolic-ID decomposed). `panel/pangenie_genotyping/data/merged/founders_231_chr.vcf.gz` is the pre-Arch 3 mixed-ploidy catalog *(archive)*; see `panel/pangenie_genotyping/data/merged/README_GOLDEN_STANDARD.md`.
- **GrENE-Net SNP VCF** (`greneNet_final_v1.1.recode.vcf`): 231 founders, ~3.24M SNPs. Used by hapFIRE in the methods-comparison column and as a `cn_var` second-source in the Fix 2 hybrid `cn_var`. Note: the VCF header lists 232 names but the 232nd is blank.
- **Cactus pangenome**: `/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.vcf.gz` — minigraph-cactus on 82 long-read assemblies (80 unique GrENE-Net Accession_IDs after dedup; the two known-duplicate pairs are 5772/6150 and 6915/8387 — historical labeling query at `old_docs/CACTUS_ASSEMBLY_LABELING_QUERY.md`; resolved by dropping 5772 + 9947 in v3qc).
- **PanGenie genotypes**: 148 of the 151 missing founders genotyped from ENA fastqs (PRJNA273563 + PRJNA30811), 2 from xwu BAMs (100001, 100002), 1 absent.
- **Production cn matrices**: cn_var — `panel/arch3/chr1/cn_var_231_arch3_chr1.{cn_var,cn_var_called,meta}.npz` (Arch 3). **cn_full is the OPEN k-mer-filter decision** — many candidate builds under `data/cn_full_231_v3qc_v3*` (raw, `_filt2`, `_mixedloose` [deprecated], `_subsampMedian_refilt2`, `_subsampProtect{1,2}_refilt2`, …) under active test (2026-05-26); see `docs/PIPELINE_STATE.md` §1/§3. Earlier `cn_full_231_v3`, `cn_var_231_v3` etc. are *(archive)*.
- **Existing hapFIRE outputs from Xing Wu** for the SEEDMIX samples: `/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/seed_mix/s{1..8}_ecotype_frequency.txt`. For the 2,415 evolved samples: `…/hapFIRE_frequencies/samples/ecotype_frequency/MLFH*_ecotype_frequency.txt`.

---

## References

- **hapFIRE source**: github.com/xingwu2/hapFIRE (also github.com/moiexpositoalonsolab/HapFIRE). Local copy: `/carnegie/nobackup/scratch/xwu/haplotype_frequency_estimation/hapFIRE_sourcecode/`.
- **HARP**: Kessner et al. 2013, *Mol Biol Evol* — pubmed.ncbi.nlm.nih.gov/23364324/. Binary at `…/hapFIRE_sourcecode/bin/harp` (must be on PATH for hapFIRE to invoke it).
- **HAF-pipe**: Tilk et al. 2019 — conceptual sibling using fixed windows instead of LD blocks.
- **GrENE-Net pipeline / grenedalf** for downstream allele-frequency analysis once the SV table is produced.

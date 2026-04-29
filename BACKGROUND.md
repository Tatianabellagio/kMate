# Background — why this project exists, and why the approach is what it is

This file captures the reasoning that led to the current method. `SUMMARY.md` is the results writeup; `HANDOFF.md` is the next-steps queue. Read this if you need the *why*.

---

## The problem

Estimate **structural variant (SV) allele frequencies** in evolved Pool-seq populations from the **GrENE-Net** experiment.

- 2,415 evolved Pool-seq libraries (`MLFH*`), 8 SEEDMIX founder-pool replicates.
- Founder panel: 231 *A. thaliana* ecotypes, of which **81 have long-read assemblies** (others have only short-read SNP genotypes from 1001 Genomes / GrENE-Net VCF).
- Existing SV calls: minimap + syri on the 81 assemblies → SV VCF (50,446 SVs × 80 unique 1001G IDs after collapsing duplicate assemblies).
- Existing SV frequency method: **freqk** (k-mer-based per-sample AF). Suffers at low coverage and for SVs with non-unique flanking k-mers.

Goal: a more accurate per-sample SV allele-frequency table that slots into the downstream GEA pipeline.

---

## The method, in one line

Run **hapFIRE** unmodified on SNPs, then project its founder-frequency output onto SVs via the founder × SV genotype matrix:

```
f_SV(v) = Σ_e G_SV(e, v) · f_ecotype(e)
```

Mathematically identical to how hapFIRE itself computes per-SNP frequencies internally (`snp_frequency = haplotype_freq @ haplotype_to_SNP_matrix`). ~50 lines of post-processing Python on top of unmodified hapFIRE.

---

## Why hapFIRE / HARP, and why we don't modify them

**hapFIRE** (Wu, Exposito-Alonso lab — github.com/xingwu2/hapFIRE, also github.com/moiexpositoalonsolab/HapFIRE): infers SNP and founder-accession frequencies from Pool-seq BAMs given a phased founder panel.
- Uses **HARP** (Kessner et al. 2013, PMID 23364324) under the hood for windowed haplotype-frequency estimation.
- Two-stage internally:
  1. **HARP** estimates haplotype frequencies in each LD block from the BAM, using SNPs.
  2. **CVXPY** projects haplotype frequencies onto founder accession frequencies via constrained least-squares (`ecotype_frequency_estimation_selected`), using the founder × haplotype dosage matrix.
- Per-SNP frequencies are then `haplotype_freq @ haplotype_to_SNP_matrix`.
- Conceptually similar to **HAF-pipe** but uses LD-defined haplotype blocks instead of fixed windows.
- Documented in the GrENE-Net paper (Text S4); described as a tool that exploits whole-genome linkage so effective coverage exceeds true coverage.

### Why we cannot just feed SVs into HARP

HARP's likelihood is **per-base**: each read contributes `P(observed base | true base, qual)` at SNP positions, summed over reads. The SNP-table format (`snp_table_construction`) literally has A/C/G/T columns per haplotype per position. SVs don't fit that model — a deletion or insertion isn't a base substitution; the evidence is split reads / discordant pairs / coverage shifts. Making HARP SV-native would require replacing the entire per-base likelihood with a structural-evidence model — a real rewrite.

### Why we don't need to

hapFIRE *already* outputs the founder-frequency vector. SNP frequencies fall out by linear projection. Nothing about that projection is SNP-specific — swap the marker matrix from `founder × SNP` to `founder × SV` and the same machinery yields SV frequencies. **The estimator is mathematically identical to what hapFIRE does for SNPs.**

### Alternatives considered and ruled out

1. **Modify HARP's likelihood for split-read / breakpoint evidence.** Real engineering, brittle, unnecessary given option 3 works.
2. **Pangenome-graph approach** (vg map / giraffe; founders = paths through minigraph-cactus graph). More elegant, uses existing pangenome infrastructure, but it's a rewrite not an adaptation.
3. **Marker-substitution / projection** ← what we're doing. Smallest change, exact projection (no estimation noise added at the SV step), reuses unmodified hapFIRE.
4. **freqk-style per-SV estimates with founder SV genotypes as a design matrix.** Discards HARP's read-level information; loses the whole point of using LD.

---

## Known limitations of this approach

In rough priority for downstream GrENE-Net analysis:

1. **No-recombination assumption.** hapFIRE assumes every haplotype in the pool is one of the listed founders. The simulation does not test this (VISOR pools intact founder consensus FASTAs). Evolved GrENE-Net samples have 1–3 generations of recombination → expect degraded per-block accuracy. This is the assumption that *matters* and is *untested*.
2. **Panel ascertainment ceiling.** Only 80 of 231 GrENE-Net ecotypes have SV genotypes (≈33% of seedmix mass, 17.7% by recipe). SVs polymorphic only among the un-genotyped 151 are invisible to projection. Estimates for SVs that segregate in both groups will be biased low proportional to mass on un-assembled founders.
3. **231-simplex identifiability.** With many similar founder haplotypes, CVXPY can have multiple near-optimal solutions. Per-ecotype Pearson r between hapFIRE freqs and seedmix recipe is ~0; aggregates (sums over multiple carriers) are accurate.
4. **Per-window vs genome-wide.** hapFIRE currently saves only the genome-wide weighted average ecotype frequency (top-20% by haplotype diversity). Per-block freqs are computed but discarded (`hapFIRE.py:138-139` are commented out). Patching to use per-block freqs would mitigate (1) by isolating recombinants to their own blocks.

---

## Data provenance

- **GrENE-Net SNP VCF** (`greneNet_final_v1.1.recode.vcf`): 231 founders, used unmodified by hapFIRE. Chr1-only subset at `data/vcf/greneNet_Chr1_only.vcf`. Note: the VCF originally has a blank trailing column — header lists 232 names, only 231 are real.
- **Founder SV VCF** (50,446 SVs): minimap + syri on **81 long-read assemblies**, then bcftools merge. Pipeline at `/home/tbellagio/scratch/pang/sv_panel/`. Collapses to 80 unique 1001G IDs (one accession had two assemblies).
- **Founder × SV dosage matrix** (`data/founder_sv_matrix.parquet`): 50,446 × 82, 0/1 dosage. Built by `scripts/build_founder_sv_matrix.py`. *A. thaliana* inbreds are effectively haploid for the purposes of this matrix.
- **Assembly inventory** (`ASSEMBLIES_Best_version_of_dataset.csv`): 534 unique Accession_IDs total across all available long-read data; only **84** overlap GrENE-Net 231 (so a bigger panel from current data buys at most +4 over the 80 we have).
- **Existing hapFIRE outputs from Xing Wu** for the SEEDMIX samples: `/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/seed_mix/s{1..8}_ecotype_frequency.txt`. For the 2,415 evolved samples: `…/hapFIRE_frequencies/samples/ecotype_frequency/MLFH*_ecotype_frequency.txt`.

---

## References

- **hapFIRE source**: github.com/xingwu2/hapFIRE (also github.com/moiexpositoalonsolab/HapFIRE). Local copy: `/carnegie/nobackup/scratch/xwu/haplotype_frequency_estimation/hapFIRE_sourcecode/`.
- **HARP**: Kessner et al. 2013, *Mol Biol Evol* — pubmed.ncbi.nlm.nih.gov/23364324/. Binary at `…/hapFIRE_sourcecode/bin/harp` (must be on PATH for hapFIRE to invoke it).
- **HAF-pipe**: Tilk et al. 2019 — conceptual sibling using fixed windows instead of LD blocks.
- **GrENE-Net pipeline / grenedalf** for downstream allele-frequency analysis once the SV table is produced.

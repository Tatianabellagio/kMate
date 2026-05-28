# Simulation methods

Methods-ready description of the pool-seq simulation framework used to benchmark `kMate` (and comparators). This document describes the simulation pipeline (how pools and reads are generated, how truth is computed) and the regime matrix tested.

The code lives in two places:

- **Production sim driver and helpers** — `sims/visor_freqk/scripts/`
  - `make_recomb_mosaics.py` — mosaic founder FASTA generator
  - `compute_recomb_truth.py` — per-record AF truth from ancestry tracks
- **Active regime sweep** — `benchmarks/p80/scripts/`
  - `06_run_sim_p80.sh` — uniform-fraction pool driver
  - `06b_run_sim_p80_skewed.sh` — dominant-individual ("selection-like") variant
  - `make_recomb_mosaics_p80.py` — local clone of the mosaic builder with `--source-weights` support

## 1. Overview

Each simulated pool-seq dataset is generated in three stages:

1. **Pool construction** (Stage 1): draw `N` founder genomes from a source population (multinomial), optionally apply `G` generations of meiotic recombination at hotspot-aligned crossover positions, and stitch consensus FASTAs at recombination breakpoints to produce `N` per-individual haploid mosaic genomes.
2. **Short-read simulation** (Stage 2): pool the `N` mosaic genomes at specified clone fractions, simulate paired-end 150 bp Illumina-like reads with VISOR SHORtS (Bolognini et al. 2020) at target cumulative coverage `D`.
3. **Truth computation**: derive the per-record realized-pool allele frequency from the ancestry tracks via missing-at-random projection through the founder × variant matrix (`cn_var`).

Both binomial sampling layers of real pool-seq are explicit: the founder draw at Stage 1 and the read sampling within VISOR at Stage 2.

## 2. Pool construction

### 2.1 Source population

The source population is parameterised by a probability vector `p_f` over the `F` available founders. The default is uniform (`p_f = 1/F` for all f); a non-uniform vector may be supplied via `--source-weights <TSV>` to model differential fitness ("selection"). When non-uniform, the source-truth — the population-level allele frequency at infinite pool size — is `source_truth_af[r] = (Σ_f p_f · cn_var[f, r]) / (Σ_f p_f · cn_var_called[f, r])` (see `compute_recomb_truth.py`).

### 2.2 Gen-0 individuals

For each of `N` pool individuals, one founder is drawn from `p_f` by independent multinomial sampling. When `N > F`, draws are with replacement (each founder can be assigned to multiple individuals). The realised pool composition is recorded in `pool_weights.tsv` (founder, count, weight). Each individual at gen 0 carries a single founder's whole-chromosome ancestry.

### 2.3 Recombination (g ≥ 1)

For `n_generations ≥ 1`, each generation rebuilds the population: every new individual is formed by pairing two parents drawn independently and uniformly from the previous generation, then applying meiotic recombination chromosome-by-chromosome:

- Number of crossovers on chromosome `c` is `~Poisson(L_c × r)` where `L_c` is the chromosome length and `r` is the per-bp recombination rate. The default `r = 4×10⁻⁸/bp/meiosis` matches the A. thaliana literature consensus of ~4 cM/Mb (Salomé et al. 2012). On Chr1 (30.4 Mb) this is ~1.2 crossovers per meiosis in expectation.
- Crossover positions are sampled from a precomputed set of **LD-block boundaries** (BigLD on the GrENE-Net 231-panel; `hapfire_block_index_chr1.npz`), modelling the empirical concentration of A. thaliana recombination at hotspots. A uniform-position alternative is supported (`--no-hotspots`) but is not the default.
- The offspring's chromosome alternates parent identity at each crossover; the resulting ancestry track is a sequence of (start, end, founder) intervals.

At `n_generations = 0` the pool is the gen-0 set (no recombination; each individual is a single founder). At `n_generations = 3` individuals can derive ancestry from up to 8 distinct founders along the chromosome.

A. thaliana is treated as effectively haploid because the GrENE-Net populations are >97% selfing (Exposito-Alonso et al. 2018); each individual carries one mosaic haploid genome. The two diploid copies in real individuals are near-identical, and `cn_var` is stored per founder, not per haplotype.

### 2.4 Mosaic FASTAs

For each individual, founder-specific consensus FASTAs (one per cactus founder, generated previously via `bcftools consensus` against the production panel VCF) are concatenated at the ancestry-track breakpoints to produce a single mosaic haploid FASTA. The mosaic FASTAs are written to `haps/s_indNNN/h1.fa` for VISOR consumption. An accompanying `ancestry.tsv` records `(ind_id, chrom, start, end, founder)` for downstream truth computation.

## 3. Read simulation

VISOR SHORtS (Bolognini et al. 2020) simulates paired-end Illumina-like reads from the pool of mosaic FASTAs. Default parameters:

| Parameter | Value | Note |
|---|---|---|
| Read length | 150 bp | Standard for current Illumina chemistry |
| Sequencing error rate | 0.001 | VISOR default; lower than typical Illumina (~0.5%) — provides a low-noise read substrate so the inference floor is not dominated by base-call error |
| Insert size | VISOR default | Empirical Illumina-like distribution |
| Target coverage | 10× cumulative across the pool | Equivalent to per-individual ~0.2× at `N=50` |
| Reference | TAIR10 (`TAIR10.chr.iupacN.fa`) | IUPAC ambiguity codes replaced with `N` to match cactus VCF convention |

Each mosaic clone is assigned a fraction of the total pool reads via the `--clonefraction` argument. Two modes:

**Uniform** (default): every individual gets `100/N` percent of pool reads. Where `100/N` is not finite-precision representable (e.g. `N = 231`), the first individual absorbs the rounding residual via an iterative fixed-point adjustment to satisfy VISOR's strict-equality `sum == 100.0` check (see `06_run_sim_p80.sh`).

**Skewed ("selection-like")**: one designated individual receives `DOMINANT_FRAC` percent of pool reads (default 50%), the remaining `N-1` individuals split the rest equally. Breaks the ergodicity assumption of the uniform pool: a single ancestry pattern dominates the chromosome-wide signal, producing per-position founder mixtures that are spatially non-stationary. Tests robustness of AF estimators to over-represented individuals.

## 3.1 Read substrate: VCF-consensus vs raw-assembly reads (validation)

Reads in this framework are simulated from **VCF-consensus** founder FASTAs — `bcftools consensus` of the panel VCF applied to TAIR10 (§2.4) — not from the raw long-read assemblies. This is a deliberate choice with one obvious risk and one hard constraint, and it has been validated directly.

**The risk (read-side closed loop).** Consensus reads contain only the variants the panel VCF encoded, and `cn_full` (the EM's k-mer dictionary) is built from that same VCF. So the simulated read substrate is a guaranteed subset of the estimator's k-mers: no read carries a k-mer `cn_full` hasn't seen. Real pool-seq reads do — real genomes carry variation and repeat content the pangenome graph never captured. If that off-panel k-mer mass degraded the EM, a consensus-only simulation would hide it.

**The validation (g0).** We tested this on the two non-recombinant regimes by re-simulating reads **straight from the raw assemblies** (`…/pang_1001gplus/…/chr_only/<asm_id>.chr.fa`, via `pywgsim`, bypassing the VCF→consensus→VISOR path entirely; `benchmarks/p80/scripts/sim_from_raw_assemblies.py`) and re-running `kMate global` against the *same* `cn_full`, `cn_var`, and truth. Only the read source differs. Both substrates used the same 10× read budget (verified: 9.95× consensus / 10.00× raw of TAIR10 Chr1). Raw reads carry extra off-panel k-mer mass relative to consensus — the EM's λ̂ coverage estimate read ~7.2× for raw vs ~7.6× for consensus from the *identical* budget (λ̂ counts only k-mers present in `cn_full`; the ~7.x value is itself a uniform-h / partly-represented-panel artifact seen for both substrates, not missing depth).

Per-record AF error vs the same realized-pool truth (`benchmarks/p80/scripts/compare_af_vs_truth.py`):

| regime | class | MAE consensus | MAE raw | R² consensus | R² raw | \|d\|>0.1 raw |
|---|---|---|---|---|---|---|
| n50_g0  | ALL | 0.0048 | 0.0055 | 0.9985 | 0.9980 | 0.0006 |
| n50_g0  | SNP | 0.0049 | 0.0055 | 0.9985 | 0.9981 | 0.0006 |
| n50_g0  | SV  | 0.0028 | 0.0035 | 0.9983 | 0.9961 | 0.0012 |
| n231_g0 | ALL | 0.0048 | 0.0047 | 0.9986 | 0.9984 | 0.0008 |
| n231_g0 | SNP | 0.0048 | 0.0047 | 0.9986 | 0.9985 | 0.0008 |
| n231_g0 | SV  | 0.0030 | 0.0033 | 0.9980 | 0.9964 | 0.0014 |

**For g0, raw-assembly and VCF-consensus reads reach the same result** — ALL-class MAE within ±0.0007 and R² within 0.0005 at both pool sizes; at the production-relevant n231 "perfect mix" the raw substrate is in fact marginally *better* (MAE 0.0047 vs 0.0048, less bias). The only visible cost of the honest substrate is a tiny new outlier tail (≤0.08% of records, vs exactly 0 for consensus) and a proportionally larger but still small SV-RMSE penalty. The EM absorbs the extra off-panel k-mer mass real assemblies carry with no meaningful loss in g0 AF accuracy.

**Why we keep the VCF-consensus simulation (the recombination constraint).** For the recombinant regimes (g1, g3, dom500) the consensus path is not just convenient — it is required. A recombinant individual is a mosaic whose breakpoints are defined in **TAIR10 coordinates** (LD-block boundaries; §2.3). Consensus FASTAs are already in TAIR10 coordinates (TAIR10 + that founder's variants, near-identical length), so a mosaic is stitched by slicing each consensus FASTA at the TAIR10 breakpoint — exact and trivial. Raw assemblies are in their **own** coordinate frame (Chr1 ranges 29.5–34.6 Mb vs TAIR10's 30.4 Mb), so stitching a TAIR10-defined mosaic from raw sequence requires an assembly→TAIR10 alignment (e.g. minimap2) to project every breakpoint, plus a few-bp approximation at each crossover. Since g0 establishes that the read substrate is **not** the accuracy bottleneck, that added machinery would buy nothing for the recombinant regimes. **We therefore use VCF-consensus reads throughout the regime sweep, validated against raw-assembly reads at g0.**

The truth side is a separate matter: truth here is still computed through `cn_var` (§4), so this validation isolates the **read substrate**, not `cn_var` correctness. An independent, alignment-based truth would be required to close the truth-side loop and is out of scope for this framework.

## 4. Truth computation

Truth allele frequencies are computed deterministically from the ancestry tracks and the founder × variant matrix `cn_var` (and its called-mask `cn_var_called`). For each variant record `r`, the missing-at-random projection over the realised pool is:

$$
\text{truth\_af}[r] = \frac{\sum_i w_i \cdot \text{cn\_var}[\text{anc}_i(r),\ r]}{\sum_i w_i \cdot \text{cn\_var\_called}[\text{anc}_i(r),\ r]}
$$

where `w_i` is the pool-weight of individual `i`, `anc_i(r)` is the founder owning the ancestry segment containing record `r`'s position in individual `i`, and `cn_var[f, r]`, `cn_var_called[f, r]` are the carrier and called indicators (1 if founder `f` carries the ALT at record `r` / has a non-missing GT, 0 otherwise).

This projection is the same form `kMate` uses to map the EM-inferred founder mass vector to per-record AF, applied here with *exact* per-individual ancestry rather than an inferred mixture. Conceptually it is the optimal AF a perfect estimator could recover from the realised pool — independent of read sampling noise.

The truth file also reports `info[r] = Σ_i w_i · cn_var_called[anc_i(r), r]` — the h-weighted fraction of pool mass observable at record `r`. At records where `info → 0` (no individual's ancestry-founder is called), `truth_af` is emitted as `NaN`.

Optionally (when `source_weights.tsv` is present), a second truth column `source_truth_af` reports the source-population AF (`p_f · cn_var` / `p_f · cn_var_called`) — the population parameter the AF estimator is ultimately trying to recover. The gap `|truth_af − source_truth_af|` quantifies the Stage-1 sampling + drift noise floor — irreducible by any AF estimation method given the realised pool.

## 5. Regime matrix

The current canonical regime sweep (benchmarks/p80, Chr1, 10× coverage):

| Regime | N | G | Pool fractions | What it tests |
|---|---|---|---|---|
| `n50_g0` | 50 | 0 | uniform 2% | Small subset, no recombination — clean per-window founder identity, baseline for the `h chromosome` estimator |
| `n50_g1` | 50 | 1 | uniform 2% | Small subset + 1 generation recombination — moderate mosaicism |
| `n50_g3` | 50 | 3 | uniform 2% | Small subset + 3 generations recombination — fine mosaic, stresses `h chromosome` |
| `n231_g0` | 231 | 0 | uniform 100/231 % | "Perfect mix" approximation — N >> F, every founder represented ~2.9×, no recombination; cleanest possible pool |
| `n231_g1` | 231 | 1 | uniform 100/231 % | Perfect mix + 1 generation recombination |
| `n50_g3_dom500` | 50 | 3 | ind001 = 50%, others = ~1.02% each | Selection-like: one individual dominates pool reads. Hardest regime; the dominant individual's specific ancestry breaks ergodicity. |

Regimes are ordered (in result tables) easiest → hardest by empirical realised-pool truth MAE. The same naming scheme generalises: any (`n<N>_g<G>` ± `_dom<F>`) tuple identifies a simulation. The 231-panel production benchmark uses the same regimes against the production `cn_var` and `cn_full` matrices.

## 6. Output schema

Each simulation produces a self-contained directory `sims/cov<COV>_n<N>_g<G>_s<SEED>_hotspots[_dom<F>]_<panel>_chr<C>/`:

```
ancestry.tsv               (ind_id, chrom, start, end, founder)
pool_weights.tsv           realised pool composition (founder, count, weight)
source_weights.tsv         source-population p_f (defaults to uniform 1/F)
visor_pool_fractions.tsv   per-individual VISOR clonefractions (skewed-mode only)
haps/s_indNNN/h1.fa[.fai]  per-individual mosaic FASTA
reads/r1.fq, r2.fq         pooled paired-end reads
recomb_truth.tsv.gz        chrom, pos, ref_len, alt_len, truth_af, info[, source_truth_af, source_info]
region.bed                 VISOR genomic-region definition
```

The truth TSV is the join target for evaluating estimator outputs (`alt_freq`, `info`, `n_called`, `se` from `per_sample_per_chrom.py` — see `PIPELINE_STATE.md`).

## 7. Parameter defaults and configurable knobs

| Parameter | Default | Where to change |
|---|---|---|
| Recombination rate `r` | `4×10⁻⁸/bp/meiosis` | `--recomb-rate` on `make_recomb_mosaics_p80.py` |
| Crossover model | hotspot-aligned (LD-block boundaries) | `--crossovers-from-ld-blocks <npz>` (default) or `--no-hotspots` for uniform-position Poisson |
| Read length | 150 bp | `--length` in `VISOR SHORtS` |
| Sequencing error | 0.001 | `--error` in `VISOR SHORtS` |
| Coverage | 10× cumulative | `--coverage` in the sim driver (positional arg 3) |
| Founder probabilities | uniform 1/F | `--source-weights <TSV>` on `make_recomb_mosaics_p80.py` |
| Pool dominance fraction | 50% (skewed mode only) | `DOMINANT_FRAC` env var, fourth positional arg in `06b_run_sim_p80_skewed.sh` |
| Chromosomes | `Chr1` only currently | `CHROMS=…` env var in driver; controls VISOR region BED, mosaic builder, and downstream `--chroms` flag |
| Seed | 42 | Last positional arg in the sim driver |

## 8. Modelling caveats

1. **Recombination is forced outcrossing**, not Arabidopsis-realistic. Every individual at generation `g+1` is built from two random parents at generation `g`. Real A. thaliana selfs ≈99% of the time; outcrossing meioses are ~1% per generation. Our `G=3` therefore approximates several hundred generations of natural drift under empirical selfing rates. The simulation deliberately tests high-recomb regimes; results should be interpreted as upper-bound mosaicism, not realistic field dynamics.
2. **No drift between source and pool beyond multinomial sampling.** The source population is treated as a fixed `p_f` vector; we do not simulate within-generation Wright–Fisher drift on the source. To model this, replace the fixed `p_f` with a Wright–Fisher trajectory over `T` generations and condition Stage 1 sampling on the final state (would integrate naturally with MimicrEE2; not implemented here).
3. **No within-individual heterozygosity.** Each individual carries one mosaic haploid; diploid heterozygosity is collapsed under the high-selfing assumption. For non-selfing systems a per-individual paired-haplotype mode would be needed.
4. **Sequencing error rate is below typical Illumina.** `0.001` is VISOR's clean default. Real Illumina is closer to `0.005`. The clean substrate isolates the inference floor due to founder-panel + projection geometry rather than basecaller noise.

## 9. Tools and references

- **VISOR SHORtS**: Bolognini et al. 2020, *Bioinformatics* 36:1267. https://github.com/davidebolo1993/VISOR — paired-end pool-seq read simulation from haplotype FASTAs.
- **A. thaliana recombination rate**: Salomé et al. 2012, *Heredity* 108:447 — empirical genetic map yielding ~4 cM/Mb.
- **A. thaliana selfing rate**: Exposito-Alonso et al. 2018, *Mol Ecol* — >97% selfing in natural populations; basis for the haploid-individual assumption.
- **LD-block hotspots**: BigLD partitioning (Kim et al. 2018, *Bioinformatics* 34:359) on the GrENE-Net 231-panel SNP set (`hapfire_block_index_chr1.npz`).
- **Pool-seq sampling theory**: Futschik & Schlötterer 2010, *Genetics* 186:207 — variance decomposition `Var(p̂) ≈ p(1-p)[1/(2N) + 1/D]` characterising the Stage-1 + Stage-2 floor.
- **MimicrEE2** (referenced for comparison; not used here): Vlachos & Kofler 2018, *PLoS Comp Biol* — full forward-time pool-seq simulator with drift and selection.
- **kMate** projection math: `ALGORITHM.md` (`old_docs/CACTUS_EM_MATH.md` is the superseded version).

## 10. Reproducibility

- Python ≥3.10 in conda env `hapfm`; VISOR + `samtools`/`bcftools` in conda env `pang` (paths at top of driver scripts).
- All randomness is controlled by the seed argument. With the same seed + parameters, sims reproduce bit-identically.
- All sim output paths are deterministic from `(coverage, N_INDIV, N_GEN, SEED, panel, chrom, optional DOMINANT_FRAC)`.

---

For results, evaluation methodology, and the regime-by-method MAE tables, see `benchmarks/p80/results/FINAL_RESULTS_cov10_p80.ipynb` (notebook) and the `benchmarks/p80/results/plots/` figure set. The 231-panel regime sweep on the production cn_var/cn_full is the next planned step (see `PIPELINE_STATE.md` §3).

# Pool-seq frequency model — specification

**Status:** design draft, 2026-04-27. Not yet implemented.
**Reads:** `BACKGROUND.md` (why we're doing this), `CACTUS_PANGENOME_AUDIT.md` (the panel). Cite: PanGenie (Ebler 2022), HARP (Kessner 2013), HAF-pipe (Tilk 2019).

## Goal

Estimate per-SV (and per-SNP) allele frequencies in a pool-seq sample, given:
- a phased multi-sample panel VCF (cactus pangenome, top-level snarls)
- pool-seq short reads from N evolved samples
- a per-sample expected coverage `λ_i`

The output is a per-sample × per-variant frequency table that slots into the downstream GEA pipeline.

## Why a custom model

PanGenie is per-individual diploid only. HARP/hapFIRE handle pool-seq but only on SNPs and don't use graph-aware k-mers. HAF-pipe handles pool-seq on SNPs with LD windows. Nothing existing combines: pool-seq + graph-aware k-mer evidence + LD borrowing + multi-allelic SVs. We need to build it. The pieces are well-validated in adjacent contexts; the combination is what's novel.

## Inputs

| Symbol | Meaning | Shape | Source |
|---|---|---|---|
| `F` | founders in panel | scalar (82 → 231 after imputation) | cactus seqfile |
| `B` | LD blocks genome-wide | scalar (~30) | computed from bubble positions via Li–Stephens-style partition |
| `K_b` | unique k-mers in block b | varies (10² to 10⁴) | PanGenie-index output |
| `J_v` | alt alleles of bubble v | varies (1–80) | cactus VCF |
| `cn_k ∈ {0,1}^F` | whether each founder carries k-mer k | F-vector per k-mer | derived from PanGenie-index unique-kmer assignment to alts × cactus genotype (subject to PanGenie's per-allele cap — see below) |
| `c_{ik}` | observed count of k-mer k in sample i | scalar | jellyfish count of pool-seq BAM |
| `λ_i` | global k-mer coverage for sample i | scalar | total reads × read length / genome size |
| `enc[v, j, s] ∈ {0,1}` | does alt j of bubble v encode panel SV s | sparse matrix | computed via vcfwave or sequence-alignment of alts vs ref |

### PanGenie k-mer index caveats (relevant for our cn matrix)

We build `cn` from PanGenie's `kmers.tsv.gz` directly, so PanGenie's index-time
filtering is **inherited** by the cactus_em pipeline. Three filters apply
per pangenome bubble:

1. **Internal var-dedup:** a k-mer must be unique within the bubble's alleles
   (PanGenie source: `if (local_count > 1) continue`). Equivalent to our
   freqk var-dedup step.
2. **Internal ref-dedup:** a k-mer must not appear elsewhere in the graph or
   reference (`if (genomic_count - local_count != 0) continue`). Equivalent
   to freqk ref-dedup.
3. **Hard cap:** **16 k-mers per allele** for biallelic bubbles, **32 per allele**
   for multi-allelic. Biallelic bubbles thus carry up to 32 k-mers (ref + alt);
   N-allelic bubbles up to 16·N.

Distribution observed in `cn_full_Chr1.meta.npz` (206,631 bubbles, 20.8M total
k-mers, mean 100.6/median 84/max 301):

| n_kmers/bubble | n_bubbles | % | Note |
|---|---|---|---|
| 1-15 | 1,284 | 0.6 | Fewer than the cap — allele has limited unique k-mers |
| **= 16** | 2,291 | 1.1 | Biallelic-cap saturated |
| 17-31 | 1,971 | 1.0 | |
| **= 32** | 65,681 | 31.8 | Biallelic, both alleles at cap (16 + 16) |
| 33-100 | 60,517 | 29.3 | Multi-allelic, several alleles × 16 |
| 101-301 | 74,887 | 36.2 | Multi-allelic with many alleles |

**Implication for coverage analysis:** with ≥4 k-mers/bubble (99.7% of bubbles)
and ~80M unique k-mers genome-wide, the EM is massively over-determined for a
231-vector founder simplex even at 5× read coverage. This is consistent with
the observed cov 5x ≡ cov 30x saturation in `pool_sweep_82` (Tier 3 of the
simulation synthesis). The cap removes redundancy, not signal — every bubble
still has enough k-mers to discriminate its allele under the simplex prior.

**Implication for very low coverage (≪ 1×):** at ~0.5× whole-genome coverage,
expected hits per k-mer ≈ 0.5; per allele with cap=16 → expected ≈8 hits per
allele, which is the threshold where the cap can start to limit estimation.
Realistic GrENE-Net evolved data is 5–15× so this is well above the cap regime.

**Implication for downstream PanGenie genotyping (the `pangenie_genotyping/`
sister pipeline):** the cap is a built-in PanGenie behaviour. If a bubble's
most diagnostic k-mers are excluded by the cap, per-individual allele dosage
calls inherit that noise. Not a concern for our cactus_em pipeline (which
uses the same set), but worth flagging for the direct-genotyping route.

## Generative model

For sample i, k-mer k in block b, the count is Poisson with mean determined by founder frequencies:

```
E[c_{ik}] = λ_i · α_k · h_{b,i}^T · cn_k  +  ω_k · λ_i
            └────── true SV/SNP signal ──────┘   └─ contamination ─┘

c_{ik} | h_{b,i}, ω_k ~ Poisson(E[c_{ik}])
```

- `h_{b,i} ∈ Δ^F` is the **founder frequency vector for block b in sample i** (the latent we want)
- `α_k` is a per-k-mer "reach" (≈ 1 for clean breakpoint-spanning k-mers; can be set to 1 initially, learned later)
- `cn_k[f]` is 1 if founder f's haplotype carries k-mer k, 0 otherwise (precomputed from PanGenie-index)
- `ω_k` is per-k-mer **off-target contamination rate**, sample-independent at normalized coverage

The contamination term `ω_k` is what lets us avoid strict k-mer dedup. PanGenie keeps only unique k-mers and discards everything else; we keep all k-mers and let the bilinear model attribute counts to true vs off-target sources.

## Constraints

```
h_{b,i} ≥ 0
Σ_f h_{b,i}[f] = 1
ω_k ≥ 0
```

## Identifiability

Two sources combine to give identification:

1. **Across-k-mer variation within an SV / bubble.** K-mers within one bubble share the same `h_b` but have different `cn_k` and different `ω_k`. With ≥3 k-mers of varying contamination, we can separate `h_b` from the per-k-mer offsets.

2. **Across-block constraints from SNP-rich neighbors.** SNP-rich LD blocks have `K_b ≫ F`, so the SNP rows alone heavily constrain `h_b`. Adjacent SV-rich blocks inherit identification through LD smoothing.

For block b with K_b k-mers and F founders:
- K_b ≥ F + 1: over-determined, h_b uniquely identified up to noise
- K_b < F: under-determined (rare for SNP-rich blocks; possible for small SV-only blocks)
- Mitigation: solve adjacent blocks jointly with a smoothness penalty, or use HMM-style transitions

## Estimation: per-sample, per-block CVXPY

Poisson NLL is not DCP-compliant in CVXPY. Two clean workarounds, both DCP-valid:

### Option A: Anscombe-transformed weighted L2 (one-shot, faster)

```
y_k    = 2 · sqrt(c_{ik} + 3/8)
μ_k(h) = λ_i · α_k · h^T · cn_k + ω_k · λ_i
ŷ_k(h) = 2 · sqrt(μ_k(h) + 3/8)

minimize  Σ_k (y_k - ŷ_k(h))²
s.t.      h ≥ 0,  Σ h = 1
```

Approximates the Poisson likelihood. Variance is ~1 in the transformed space. CVXPY-DCP via `cp.sqrt(cp.pos(...))`.

### Option B: IRLS (iteratively reweighted least squares)

```
Initialize h^(0) = uniform 1/F
For t = 1, 2, ...
    μ_k = λ_i · α_k · cn_k^T · h^(t-1) + ω_k · λ_i
    w_k = 1 / sqrt(μ_k + ε)
    Solve: minimize  Σ_k w_k² · (c_{ik} - λ_i · α_k · cn_k^T · h - ω_k · λ_i)²
            s.t.     h ≥ 0,  Σ h = 1
    Update h^(t)
Until ||h^(t) - h^(t-1)|| < tol
```

Equivalent to Newton on Poisson NLL. Converges in 3–5 iterations typically.

**Pick A for first prototype.** Switch to B if accuracy needs improvement.

## Estimating ω: alternate solve

Bilinear in (h, ω). Decomposable:

```
Step 1 (given ω fixed): solve per-(sample, block) CVXPY for h
Step 2 (given all h's fixed): for each k-mer k,
    ω_k = max(0, mean over samples of (c_{ik} - λ_i · α_k · cn_k^T · h_{b,i}) / λ_i)

Iterate until convergence (typically 2–3 outer passes).
```

Identification of ω requires across-sample variation in `f`. For evolved GrENE-Net data with 2,415 samples spanning sites/years, plenty of variation.

## Decomposition: bubble alt freqs → per-SV freqs

Output of the per-block solve gives `h_{b,i}` (founder freqs). To get per-bubble alt-allele freqs:

```
For each bubble v in block b, alt j:
    p_alt(v, j, i) = Σ_f h_{b,i}[f] · I(founder f's haplotype carries alt j of bubble v)
```

The `I(founder f carries alt j of bubble v)` matrix is precomputed from the cactus VCF (each founder's GT picks one alt per bubble).

Then for each panel SV s contained in bubble v:

```
f_SV(s, i) = Σ_j p_alt(v, j, i) · enc[v, j, s]
```

The encoding matrix `enc` is built one-time from the cactus VCF: for each (REF, ALT_j) pair, run vcfwave or minimap2-asm to extract embedded sub-variants, match to syri-panel SVs by position+size+sequence.

This is the analog of PanGenie's `convert-to-biallelic.py` but works for plant data (`prepare-vcf-MC` is human-only and we can't use it).

## Block partitioning

Two options:

**A. Genomic distance.** Partition bubbles into blocks by gaps >X kb (X chosen so within-block LD is high). Simple, fast, transferable across samples.

**B. Founder-haplotype LD partition (hapFIRE-style).** Compute correlation matrix between bubble alt-allele genotypes across founders; cluster into blocks where intra-block correlation > cutoff. More accurate, slower, panel-dependent.

Start with A (simple). If accuracy suffers, switch to B by adapting hapFIRE's `CompleteLDPartition` to work on bubble genotypes instead of SNP genotypes.

## Algorithm pseudocode

```python
# ============== Preprocessing (one-time) ==============
# Build founder × k-mer copy-number matrix
def build_kmer_cn(pangenie_index, cactus_vcf):
    # Each unique k-mer k belongs to a specific bubble v and identifies one or more alt alleles
    # cn[k, f] = 1 iff founder f's GT at bubble v selects an alt that contains k-mer k
    # Output: sparse F × K matrix
    ...

# Build alt-encoding matrix
def build_alt_encoding(cactus_vcf, syri_panel_vcf):
    # For each bubble v, alt j: align ALT_j vs REF, identify embedded sub-variants
    # Match each sub-variant to syri-panel SVs by pos + size + sequence
    # Output: enc[v, j, s] sparse incidence
    ...

# Block partition
def partition_blocks(bubble_positions, max_gap_kb=100):
    # Simple genomic-distance partition
    ...

# ============== Per-sample k-mer counting ==============
def count_kmers(bam, kmer_set):
    # jellyfish count -m 31 -s ... <bam to fastq pipe>
    # query the resulting hash for each k-mer in kmer_set
    # Output: dict {kmer_str: count}
    ...

# ============== Per-sample, per-block solve ==============
def solve_block(c_block, cn_block, lambda_i, alpha=1.0, omega=None):
    # c_block: K_b vector of observed k-mer counts
    # cn_block: F × K_b matrix of copy numbers
    # Returns: h_b ∈ Δ^F
    if omega is None: omega = np.zeros(c_block.shape[0])
    h = cp.Variable(F)
    mu = lambda_i * alpha * (cn_block.T @ h) + omega * lambda_i
    y_obs = 2 * np.sqrt(c_block + 3/8)
    y_pred = 2 * cp.sqrt(cp.pos(mu) + 3/8)
    obj = cp.Minimize(cp.sum_squares(y_obs - y_pred))
    cons = [h >= 0, cp.sum(h) == 1]
    cp.Problem(obj, cons).solve(solver='SCS')
    return h.value

# ============== Outer loop: alternate h, ω ==============
def fit_pool_seq_freqs(samples, blocks, cn, lambda_, max_iter=3):
    omega = np.zeros(K_total)
    h_all = np.zeros((len(samples), B, F))
    for outer in range(max_iter):
        # Step 1: solve h per sample, per block (parallelizable)
        for i, sample in enumerate(samples):
            for b, block in enumerate(blocks):
                h_all[i, b] = solve_block(c[i, block.kmers], cn[block.kmers], lambda_[i], omega=omega[block.kmers])
        # Step 2: re-estimate omega per k-mer
        for k in range(K_total):
            block_b = which_block(k)
            residuals = [c[i, k] - lambda_[i] * cn[k] @ h_all[i, block_b] for i in range(N)]
            omega[k] = max(0, np.mean(residuals) / np.mean(lambda_))
    return h_all, omega

# ============== Decomposition ==============
def decompose_to_sv_freqs(h_all, founder_alt_assignment, enc):
    # founder_alt_assignment[v, f] = which alt index founder f takes at bubble v
    # h_all[i, b, f] = founder freq for sample i in block b
    # enc[v, j, s] = does alt j of bubble v encode SV s
    
    f_SV = np.zeros((N_samples, N_panel_SVs))
    for v in bubbles:
        b = bubble_to_block[v]
        for s in panel_SVs_in_bubble(v):
            for f in range(F):
                j = founder_alt_assignment[v, f]
                if enc[v, j, s]:
                    f_SV[:, s] += h_all[:, b, f]
    return f_SV
```

## Scale and runtime

For the GrENE-Net deliverable:
- N = 2,415 samples + 8 SEEDMIX
- F = 82 (panel) or 231 (after imputation)
- B ≈ 30 LD blocks
- Total k-mers: ~10M (PanGenie-index output)
- Per-block CVXPY solve: ~5,000 SNP rows + ~50 SV rows × F=231 vars → ~seconds with SCS
- Per sample: 30 blocks × 5s = 2.5 minutes
- 2,415 samples / 24 cores: ~5 hours wall

K-mer counting is the I/O-bound step:
- jellyfish on ~10× coverage Arabidopsis BAM (~100M reads): ~5 min/sample
- Embarrassingly parallel: ~1 day on 24 cores for 2,415

## Validation strategy

Before applying to real data:

1. **Simulation (visor_freqk).** Pool 231 ecotypes uniformly, inject SVs at known frequencies, simulate reads at varying coverage. Run new pipeline. Compare to truth and to freqk + hapFIRE-projection. Already-existing simulation framework — just plug in the new method.
2. **SEEDMIX (8 samples with recipe-truth).** Run new pipeline on the 8 SEEDMIX samples. Compare per-SV freq to recipe-based panel-truth. Should match the hapFIRE-projection r=0.99 we got, ideally without the panel-ascertainment bias (because we now have SV evidence directly from k-mers, not just projection).
3. **Recombinant-pool simulation.** Build mosaic ecotype FASTAs (e.g., 2 founders + 1 crossover per generation), simulate reads, run pipeline. Tests the no-recombination assumption embedded in the per-block independence. Spec'd in HANDOFF.md as next-step #2.

If (1) and (2) match expectations, move to applying on the 2,415 evolved samples.

## Open design questions

- **Block size.** Smaller blocks = better recombination handling, but each block's CVXPY is more under-determined. Pick block size adaptively (target K_b ≥ 5·F unique k-mers per block).
- **α_k.** Setting α_k = 1 implicitly assumes every k-mer at a junction has equal "reach". Likely good enough for SNPs; for SVs the per-k-mer reach varies (e.g., k-mers near the breakpoint vs deep inside). Consider learning α_k jointly with h, ω after the simple version works.
- **Imputation.** When using F=231 (impute G_SV onto the 151 unobserved GrENE founders), the cn matrix's 151 columns are imputed probabilities, not 0/1. The Poisson likelihood naturally handles fractional cn, but identifiability is harder for imputed columns. Worth a separate analysis.
- **Multiple SVs in one bubble.** A multi-allelic bubble may encode several syri-panel SVs at once. The decomposition step should not double-count: if alt j contains both SV s₁ and SV s₂, the sample's freq of alt j contributes to both s₁ and s₂. That's correct and is what the simple `enc[v, j, s]` linear formula captures.
- **Effective sample size for finite-pool sampling.** Pool-seq has additional variance from finite-genome sampling (a pool of 100 individuals at uniform freq has 1% sampling error per allele on top of read sampling). Currently the model treats `λ_i · h_{b,i}^T · cn_k` as the rate of a Poisson. If pool size is small (100s), should use Poisson-Binomial or Negative Binomial. For large pools (≥1000 individuals) the simple Poisson is fine.

## Implementation order

1. **Build cn matrix** from PanGenie-index output + cactus VCF GTs (~2 days)
2. **Build encoding matrix `enc`** via vcfwave or minimap2-asm (~3 days)
3. **K-mer counter wrapper** around jellyfish (~1 day)
4. **Per-sample CVXPY solver** with Anscombe loss (~1 day)
5. **Outer loop** for ω alternate-fit (~1 day)
6. **Decomposition** + per-SV output (~½ day)
7. **Run on simulation framework** for validation (~½ day after pipeline is stitched)
8. **Run on 8 SEEDMIX** samples (~1 day)
9. **Apply to 2,415 evolved samples** (~1 week wall)

Total to a working end-to-end pipeline: ~2–3 weeks of focused work, assuming PanGenie-index and the cactus pangenome are stable (which they are after 56173 lands).

---

## Future work — fine-block accuracy at ~10 kb resolution from k-mers alone

**Status:** design 2026-05-06. Production target: per-block AF accuracy at **10 kb resolution** (the user's required scale; coarser is not interesting).

### Motivation

Recomb sweep (`RECOMB_SWEEP_RESULTS.ipynb`, n50_g3 hardest regime, R² SNP):

| method | unit | k-mers/unit | R² SNP | evidence type |
|---|---|---|---|---|
| global | whole chrom | ~20.8 M | 0.917 | k-mer counts |
| window_200kb | 200 kb | ~140 k | **0.946** | k-mer counts |
| bigld_panel | BigLD ~7 kb | ~61 | 0.838 | k-mer counts |
| clean_smooth | BigLD ~7 kb | ~73 800 | ~0.90 (post-smoothing) | k-mer counts |
| hapfire_fine | BigLD ~7 kb | (read-level SNP) | 0.93 | HARP read-level |

Three independent bottlenecks visible in the table:

1. **Per-block k-mer evidence sparsity** (bigld_panel's 0.838) — fixed by clean_smooth's richer per-block k-mer source. Solved.
2. **Single-h misfit when ancestry is mosaic** (global 0.917 < window_200kb 0.946 on n50_g3) — partition into blocks helps. Solved at 200 kb scale.
3. **K-mer counting destroys read-level phase** (clean_smooth 0.90 < hapfire_fine 0.93 with same partition + much more evidence) — k-mers within a single read carry the linkage info HARP exploits, but treating each k-mer count as an independent Poisson observation throws it away.

The HARP-prior hybrid path was considered and dropped: adding HARP back as a prior gives up the speed advantage that's the whole reason to prefer k-mer methods.

### Direction — three composable changes, all k-mer-only

**Production target:** **R² ≥ 0.94 on n50_g3 at 10 kb resolution.** Anything coarser than 10 kb is out of scope — this is the resolution the downstream GEA pipeline needs. Beat the 200 kb result at finer scale, not match it at coarser scale.

#### 1. Global-anchor KL-prior per-block EM **— DONE 2026-05-06, production target hit**

Run a one-shot global EM first → `h_global`. Then per-block EM uses `h_global` as a Dirichlet-style prior:
```
h_b_new ∝ h_b · (cn_b @ cw)  +  λ · N · h_global    (then renormalize to simplex)
```
where `N` is total k-mer count in the window. Local k-mer evidence has to overcome the prior to move `h_b` away from `h_global`. In well-mixed regions `h_b ≈ h_global`; in real mosaic blocks it adapts. Cost: one extra EM solve per sample (≈ +1 min wall). Code change: ~15 lines in `em_solver.py` adding the prior term + `--global-anchor-weight λ` flag.

Replaces the dropped HARP-prior idea. Same regularization shape (KL toward a prior), but the prior comes from cactus_em's own global solution — much cheaper, no separate pipeline, and the prior is already trustworthy (R² 0.97-0.99 on easy regimes).

**Result on n50_g3 cov50** (Chr1 polymorphic, n=377,523 SNPs / 11,233 big SVs):

| λ | SNP R² | SV R² | resolution |
|---|---|---|---|
| 0.0 (`window_10kb` baseline) | 0.902 | 0.807 | 10 kb |
| 0.1 | 0.955 | 0.919 | 10 kb |
| **0.3 (production default)** | **0.963** | **0.934** | **10 kb** |
| 1.0 | 0.959 | 0.931 | 10 kb |
| (`window_100kb`, reference) | 0.963 | 0.932 | 100 kb |

**λ=0.3 ties `window_100kb`'s accuracy at 10× finer resolution.** +6.1 pp SNP, +12.7 pp big-SV vs unanchored `window_10kb`. Past the production bar of R² ≥ 0.94. See RESULTS_LOG.md 2026-05-06 for the full sweep and analysis.

**Production command**:
```
--block-mode window --window-bp 10000 --global-anchor-weight 0.3
```

Implementation lives in `poolfreq/src/em_solver.py` (`solve_em` `prior_h`+`prior_weight` args) and is threaded through `block_em.solve_em_per_block`, `per_sample_driver.py`, and `per_sample_per_chrom.py`.

#### 2. Overlapping windows (k-mer kept in 2-3 windows)

Replace 4,123 disjoint BigLD blocks (or fixed disjoint windows) with **~10 kb windows stepping every ~3 kb**. Each k-mer enters ~3 windows; per-record AF is the inverse-distance-weighted average across covering windows. Implicit smoothing via overlap — adjacent windows share most k-mers and can't drift independently. Cost: ~3× per-window solves but each is on a smaller k-mer subset (~1.5× wall). No new index needed.

Composes with #1.

#### 3. Kallisto-style read-level pseudoalignment EM (highest-promise)

The current EM treats each k-mer count as an independent Poisson observation. **K-mers from the same read carry phase information** — they came from the same founder haplotype — and we throw that away when we marginalize to per-k-mer counts.

The fix is the kallisto architecture (Bray et al. 2016, validated on 10⁸ reads in minutes):

1. **Index** = per-k-mer founder bitmask. Already built — `cn_full[f, k]` is exactly that. Each k-mer column is a 231-bit "which founders carry this k-mer" set.
2. **Pseudoalignment**: for each read, intersect the bitmasks of its k-mers within a 10 kb window → "founder compatibility set" S_r. Reads sharing identical S_r are grouped into **equivalence classes**.
3. **EM on equivalence-class counts** (the sufficient statistic; per-iter cost is O(unique_ec) not O(reads)):
   ```
   h_new[f] ∝ h[f] · Σ_c n_c · I[f ∈ S_c] / (Σ_{f' ∈ S_c} h[f'])
   ```
   Same multiplicative simplex update we already have, but with EC-level counts instead of per-k-mer counts.

Why it should work for our problem: a 150 bp PE read produces ~120 overlapping 31-mers, all from one founder. Intersecting their per-k-mer compatibility sets gives a much sharper read-level set than any single k-mer's. Reads from a within-block recombinant junction will produce empty intersections → drop them (correctly identifying recombinants as out-of-panel). Reads from intact founder backgrounds will have small high-confidence S_r → strong founder pull during EM.

This is exactly the read-phase information HARP uses, but recovered from k-mers directly with no alignment needed. The architecture is well-validated: kallisto/EMIRGE/Mora/PanTax all use it for transcript/strain abundance. **Applying it to founder pool-seq frequency estimation is a real gap in the literature** — the closest published work is Tang & Anderson 2023 (microhaplotype pool method), which is read-level EM but at base/alignment level not k-mer level.

References that compose into this design:
- kallisto pseudoalignment + EC-level EM (Bray et al. 2016, Nat Biotechnol)
- Mora's set-cover regularization for similar references (Zheng 2024, BMC Bioinformatics) — relevant because our 231 founders include many near-identical pairs
- KAGE 2's variant-correlation Bayesian model (Grytten 2023) — informs how to add cross-window coupling if EC-EM still leaves gaps

### Implementation order

1. **#1 alone** at 10 kb resolution (1 day). KL-prior add to `em_solver.py` + flag in driver. Benchmark on `cov50_n50_g3_s42_hotspots_p231_chr1` against `window_10kb` and `window_200kb`.
2. **#1 + #2** (1-2 days). Overlapping window definition in `block_em.py:define_windows` + projection averaging.
3. **#3 prototype** (3-5 days). New `per_sample_kallisto_em.py` driver: per-read 31-mer extraction → per-window founder-bitmask intersection → EC counts → per-window EM. Reuse `cn_full` and `cn_var` unchanged.
4. **#3 + #1 stack** if #3 works at all. Final config = read-level EC-EM with global-anchor prior at 10 kb resolution.

### Success criterion

R² ≥ 0.94 on n50_g3 at 10 kb resolution. Stretch goal: match or beat hapfire_fine's 0.93 on n50_g3 with k-mer-only evidence — that closes the read-phase gap without paying HARP's cost.

### Risks

- #3's per-read 31-mer extraction adds a step that jellyfish currently skips. If it can't be done as fast as jellyfish counting (~10⁸ reads in minutes), the speed advantage shrinks. Mitigation: kallisto demonstrates it at scale; the index lookup is O(1) per k-mer.
- Reads with empty intersections (recombinant junctions, sequencing errors) are dropped. Need to check this doesn't bias per-block estimates against high-recomb regions.
- At 10 kb scale the per-window k-mer count is smaller than at 200 kb; EC-EM trades quantity (k-mers) for quality (read-linkage). Whether the trade is favorable is the empirical question #3 answers.

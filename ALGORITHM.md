# kMate — algorithm, math, and production wiring (code-verified)

Single source of truth for the kMate algorithm. Supersedes the earlier prose
algorithm doc and `old_docs/CACTUS_EM_MATH.md` (out-of-date math). Every formula and pipeline step
below has been cross-checked against the production code; references are
`file:line` so they stay traceable.

Code-verified against:
- `src/em_solver.py` — EM core
- `src/build_kmer_pa.py` — kmer_pa builder
- `src/build_var_pa.py` — var_pa + var_called builder
- `src/per_sample_per_chrom.py` — production per-sample driver
- `src/block_em.py` — per-window / per-block EM
- `src/kmer_count.py` — jellyfish wrapper

Last verified: 2026-05-27.

> **Note (2026-05-26 cleanup):** the estimator was reduced to two modes,
> `global` and `window`; LD-block modes, overlapping windows, the older k-mer
> rebalancing (`--row-normalize-kmer_pa`, see §8.1), carrier-weighting (§8.2) and
> contamination-ω (§8.4) were archived to `src/archive/`. §7–§8 below
> are updated, but the inline `file:line` references elsewhere predate the
> cleanup and may have shifted. The authoritative file list and invocation
> recipes now live in `src/README.md`.
>
> **Note (2026-05-27 production decision):** the kmer_pa production filter
> `filt2inv` (drop k-mers with column-sum < 2 — the original `filt2` singleton
> rule — **and** column-sum = F invariants, the latter added 2026-05-29; §2.1) and the per-k-mer EM weight
> $\omega_k = 1/m_b$ (per-bubble de-replication; §4.2) were the 2026-05-27
> production defaults on the heterogeneous 231-founder panel. filt2inv still holds.
> $\omega_k = 1/m_b$ remains the **window-mode** weighting, but for **global mode**
> (the GrENE-Net production estimator) it was **superseded 2026-07-06** by
> per-founder M-step normalization + `--kmer-weight uniform` (§4.3;
> `docs/FOUNDER_NORMALIZATION_FIX.md`). Both weights are documented inline
> below; the panel-conditional caveat for $\omega_k=1/m_b$ is in §10 (M6) and
> `docs/METHODS_TRIED_AND_RESULTS.md` §3.

---

## 1. Problem statement

Given pool-seq short reads from a sample that is a mixture of $F$ known founder
genomes (the *panel*), recover the alternate-allele frequency at every variant
in the panel VCF.

Input per sample: paired FASTQs / one BAM.
Output: per-VCF-record alt-allele frequency TSV (chrom, pos, ref_len, alt_len,
alt_freq, info, n_called, se), one row per atomized record.

For the 231-founder *A. thaliana* panel used here: $F = 231$, $K \approx 8\times10^7$
distinct panel k-mers genome-wide, $R \approx 3.4\times10^6$ atomized VCF records
covering SNPs, short indels, and SVs $\geq 50$ bp.

---

## 2. Precomputed panel matrices

Three sparse $F \times \cdot$ binary **presence/absence (incidence)** matrices are
built once per panel. Each entry is a 0/1 indicator (the panel is haploid, so
these are membership flags, **not** integer copy numbers).

**Notation ↔ implementation map.** The math symbols below are the paper/algorithm
names; the code variables, on-disk file/dir names, and CLI flags now use the same
presence/absence nomenclature throughout (the legacy `cn`/`cn_var`/`cn_full` names
were fully renamed — see `scripts/rename_cn_to_kmer_pa.py`). Use this table to bridge
the math symbols to the implementation:

| symbol | shape | what it is | code var | on-disk / CLI |
|---|---|---|---|---|
| $K_{\mathrm{pa}}$ | $F\times K$ | founder × k-mer presence/absence | `kmer_pa` | dir `kmer_pa_*`, suffix `.kmer_pa`, flag `--kmer-pa-prefix` |
| $V_{\mathrm{pa}}$ | $F\times R$ | founder × variant alt-allele presence | `var_pa` | suffix `.var_pa.npz`, flag `--var-pa` |
| $V_{\mathrm{called}}$ | $F\times R$ | founder × variant genotype-call mask | `var_called` | suffix `.var_called.npz`, flag `--var-called` |

### 2.1 $K_{\mathrm{pa}}$ (`kmer_pa`) — founder × k-mer presence/absence (the EM evidence)

$$K_{\mathrm{pa}} \in \{0,1\}^{F \times K}, \quad (K_{\mathrm{pa}})_{f,k} = \mathbb{1}[\text{founder } f \text{ carries panel k-mer } k]$$

Built per-chromosome by `build_kmer_pa.py:96-230` from PanGenie's bubble-level
k-mer index (`kmers.tsv.gz`), the panel VCF, and the reference FASTA. For each
PanGenie bubble:

1. PanGenie's k-mer list for the bubble is taken as-is. PanGenie has already
   applied: (i) in-bubble dedup, (ii) genomic ref-dedup, (iii) a cap of 16
   k-mers/allele for biallelic bubbles and 32/allele for multi-allelic bubbles.
   We inherit those filters.
2. For each founder, we reconstruct their haplotype across the bubble by
   applying their genotype to the reference (`reconstruct_haplotype`,
   `build_kmer_pa.py:45-76`).
3. We compute the canonical k-mer set of the reconstructed haplotype
   (`canonical_kmer_set`, `build_kmer_pa.py:34-42`). A canonical k-mer is the
   lexicographically smaller of $(k, \text{revcomp}(k))$; k-mers containing N
   are skipped.
4. $(K_{\mathrm{pa}})_{f,k} = 1$ iff the bubble's k-mer $k$ is in founder $f$'s
   canonical k-mer set.

**Haploid panel input — the SOLE production VCF.** The panel VCF that feeds
**both** $K_{\mathrm{pa}}$ and $V_{\mathrm{pa}}$ is the **arch3 canonical biallelic panel
`panel/arch3/chr{N}/merged_231_chr{N}_final.vcf.gz`** (produced by the arch3
A1–A5 decomposition; `panel/arch3/README.md`). This is the only production panel
VCF — the older `founders_231_v3qc*.haploid.vcf.gz` naive-`norm` merges are
archived (`docs/PIPELINE_STATE.md` §0). It is pre-haploidized and biallelic: each
GT field is a single allele `0`, `1`, or `.` (PanGenie's diploid GTs for the 153
short-read founders are het-masked `het→.` and haploidized by arch3 A2 before the
cactus-78 merge, so the 231-founder panel is entirely haploid). Under pysam this
surfaces as single-element tuples (`(0,)`, `(1,)`, `(None,)`);
`build_kmer_pa.py:189-194` reads `gt[0]`.

> K_pa is *representation-invariant* — reconstructing a founder's haplotype from
> multi-allelic GTs or their decomposed biallelic equivalents yields the same
> sequence, hence identical k-mers (~100% carrier agreement,
> `benchmarks/p231/scripts/03c_compare_kmer_pa.py`). We build K_pa from the same
> merged_231 VCF as V_pa for single-source provenance.

**Missing-GT handling (`build_kmer_pa.py:189-204`):** controlled by
`--treat-missing-as-n`. The production `kmer_pa` (now `kmer_pa_231_arch3_filt2inv`)
is built with `treat_missing_as_n=True`: `./.` → N over the variant's REF span, so every k-mer
overlapping that span is dropped from that founder's reconstructed haplotype (the
founder gets neither REF nor ALT credit at the missing site). The `./.` → REF
alternative is **deliberately rejected** — imputing reference from a low-quality
no-call fabricates a confident genotype and biases the panel toward REF.
**Verified 2026-05-29:** *both* `v3qc_v2` and `v3qc_v3` were built N-on (build
scripts pass `--treat-missing-as-n`; build logs print
`--treat-missing-as-n ON: ./. → N`), and **no N-off production build exists on
disk** — the previous claim here that v3qc_v3 was N-off was incorrect. It was
hypothesized that N-on *causes* the cactus/PG private-k-mer imbalance (PG SV
missingness → dropped PG k-mers → cactus-skewed survivors). **Tested and rejected
2026-05-29** (`notebooks/MISSINGNESS_CAUSES_IMBALANCE.ipynb`,
`scripts/run_missingness_test.py`): in *fully-called* bubbles — where N-on ≡ N-off,
so missingness cannot contribute — the Chr1 private ratio is **16.5×**, *higher*
than the all-bubble 15.5×, and the ratio is flat across PG-missingness strata
(16.3× at 0% → 15.3× at >40%). The imbalance is **real biology** (long-read cactus
assemblies realize ~16× more private k-mers/founder than short-read PG-genotyped
founders, which snap onto common graph paths and cannot call truly-private
alleles), not a missingness artifact. See §10 M3.

**Production column filter (`build_kmer_pa.py --filter-production`, applied at
generation).** The production kmer_pa is written *already filtered*: a k-mer
column $k$ with allele-count $a_k = \sum_f (K_{\mathrm{pa}})_{f,k}$ is kept iff

$$K_{\mathrm{pa}}^{\text{prod}} = (K_{\mathrm{pa}})_{:,\,\{k\,:\,2 \,\leq\, a_k \,\leq\, F-1\}}.$$

This drops three useless classes in one cut:
- $a_k = 0$ — **dead**: tag k-mers for alleles no panel founder realizes (the
  pang_135-graph-vs-231-panel mismatch; ~22% of raw columns). $\mu_k=0$, never
  matched, zero contribution.
- $a_k = 1$ — **private singletons** (the original *filt2* rule): over-represent
  founder-private sequencing errors / private repeats surviving PanGenie's
  in-bubble dedup, zero cross-founder discrimination. The threshold sweep (drop
  $a_k<2,3,5,10$) found $a_k\!\geq\!2$ optimal; filt3/5/10 over-prune
  (`docs/METHODS_TRIED_AND_RESULTS.md` §1).
- $a_k = F$ — **invariant**: carried by every founder (monomorphic in panel).
  $\mu_k=\sum_f h_f\cdot 1=1$ for all $h$, so it adds an identical constant to
  every founder's M-step term — a useless mild dilution toward the current $h$.
  Upper bound added 2026-05-29.

On Chr1 (raw 22,673,541 columns) the filter drops 5,057,036 dead + 6,465,089
private + 216,980 invariant → **10,934,436 kept (48.2%)**. The keep rule is the
single source of truth `filter_kmer_pa_production.production_keep_mask`; it is applied
in-line by the production builder `scripts/build_kmer_pa_production_arch3.sh`
(`--filter-production`), and `filter_kmer_pa_production.py` applies the same rule to a
pre-built matrix. The production matrix is **`data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr{N}`**
(built from the arch3 `merged_231_chr{N}_final.vcf.gz` + the in-house index); all
$K_{\mathrm{pa}}$ references in §3–§4 use it.

### 2.2 $V_{\mathrm{pa}}$ (`var_pa`) — founder × variant alt-allele presence (the projection target)

$$V_{\mathrm{pa}} \in \{0,1\}^{F \times R}, \quad (V_{\mathrm{pa}})_{f,r} = \mathbb{1}[\text{any allele of founder } f \text{'s GT at } r \text{ is alt}]$$

Built by `build_var_pa.py:30-109` from the **same arch3 `merged_231_chr{N}_final.vcf.gz`**
used for $K_{\mathrm{pa}}$ — **NOT** `bcftools norm -m -any`, which scatters/drops
carriers at co-located multi-allelic sites; arch3's symbolic-ID decomposition is
exactly the fix (`panel/arch3/README.md`; `docs/INVESTIGATION_CN_VAR_DECOMPOSITION.md`).
A founder is marked as an alt-carrier iff any allele in its GT field is non-zero
(`build_var_pa.py:68`: `any(a is not None and a > 0 for a in gt)`). On a
haploid GT (`(0,)`, `(1,)`, ...) this is equivalent to checking `gt[0] > 0`,
so the alt-carrier definitions in $K_{\mathrm{pa}}$ and $V_{\mathrm{pa}}$ agree exactly under the
production panel.

### 2.3 $V_{\mathrm{called}}$ (`var_called`) — founder × variant call mask

$$V_{\mathrm{called}} \in \{0,1\}^{F \times R}, \quad (V_{\mathrm{called}})_{f,r} = \mathbb{1}[\text{founder } f \text{'s GT at } r \text{ is NOT } ./.]$$

Built alongside $V_{\mathrm{pa}}$ (`build_var_pa.py:62-66`). This call mask is **critical
for the AF projection** in §6 — without it, `./.` cells get silently treated as
REF, which under-counts AF at records with high `F_MISSING` (up to 0.66 on
cactus-only or PG-only records after the cactus-78 + PG-153 merge).

---

## 3. Observation model

Per pool-seq sample, we count occurrences of every panel k-mer with Jellyfish 2
(`kmer_count.py:32-99`), invoked in canonical mode (`-C`, k=31). Reads are
streamed from BAM via `samtools fastq -F 0x900` (drops secondary +
supplementary alignments; **does not** drop duplicates — preprocessing
must handle that for non-PCR-free libraries).

Let $c_k \in \mathbb{Z}_{\geq 0}$ be the observed canonical-k-mer count and
$\mathbf{h} \in \Delta^{F-1}$ the latent founder mixture. We model

$$\boxed{\; c_k \;\sim\; \mathrm{Poisson}\!\left(\lambda \cdot \mu_k(\mathbf{h})\right), \qquad \mu_k(\mathbf{h}) = \sum_{f=1}^F h_f \, (K_{\mathrm{pa}})_{f,k} \;}$$

as conditionally independent across $k$, with simplex constraint
$h_f \geq 0$, $\sum_f h_f = 1$.

The Poisson log-likelihood
$\ell(\mathbf{h}) = \sum_k \left[c_k \log \mu_k(\mathbf{h}) - \lambda \, \mu_k(\mathbf{h})\right] + \text{const}$
is concave on $\Delta^{F-1}$.

**Weighted likelihood (production).** Each k-mer carries a non-negative weight
$\omega_k \geq 0$ that scales its contribution to the log-likelihood:

$$\ell_\omega(\mathbf{h}) \;=\; \sum_k \omega_k \!\left[\,c_k \log \mu_k(\mathbf{h}) \;-\; \lambda \, \mu_k(\mathbf{h})\,\right] + \text{const}.$$

This is the standard composite-likelihood form (Lindsay 1988): $\omega_k$
adjusts the effective contribution of each pseudo-observation without
disturbing the simplex constraint or the Poisson concavity argument. The
unweighted MLE corresponds to $\omega_k \equiv 1$. The window-mode weighting
$\omega_k = 1/m_b(k)$ — per-bubble de-replication — is derived and motivated
in §4.2 (global-mode production uses $\omega_k\equiv 1$ + per-founder
normalization; §4.3).

**Coverage estimate $\lambda$.** Reported at load time as a sanity-check value
(`per_sample_per_chrom.py:104-106`):

$$\hat{\lambda} \;=\; \frac{F \cdot \sum_k c_k}{\sum_{f,k} (K_{\mathrm{pa}})_{f,k}}$$

(total counts $\times$ $F$ divided by total founder–k-mer carrier entries;
assumes uniform $h$). **$\lambda$ cancels in the M-step ratio (§4) and is not
consumed by the EM** — the value is logged for inspection only.

---

## 4. EM algorithm

Treat each unit of count at k-mer $k$ as a latent draw from one of the
founders carrying $k$. Let $z_{k,f}$ be the (unobserved) count attributed to
founder $f$, with $z_{k,f} = 0$ if $(K_{\mathrm{pa}})_{f,k} = 0$ and
$\sum_f z_{k,f} = c_k$.

**E-step.** Given current $\mathbf{h}^{(t)}$:

$$\mathbb{E}\!\left[z_{k,f} \mid c_k, \mathbf{h}^{(t)}\right] \;=\; c_k \cdot \frac{h_f^{(t)} \, (K_{\mathrm{pa}})_{f,k}}{\mu_k(\mathbf{h}^{(t)})}$$

**M-step (weighted form, production).** Constrained MLE on the simplex under
the weighted likelihood of §3:

$$\boxed{\; h_f^{(t+1)} \;\propto\; h_f^{(t)} \cdot \sum_{k=1}^K (K_{\mathrm{pa}})_{f,k} \cdot \frac{\omega_k\, c_k}{\mu_k(\mathbf{h}^{(t)})} \;}$$

followed by L1 renormalization $\sum_f h_f^{(t+1)} = 1$. The coverage $\lambda$
cancels. Setting $\omega_k \equiv 1$ recovers the unweighted Poisson MLE
update; $\omega_k = 1/m_b(k)$ is the per-bubble de-replication weight (§4.2) —
the **window-mode** setting. (For global mode, production is $\omega_k\equiv 1$
plus the per-founder M-step normalization of §4.3; superseded 2026-07-06 —
see `docs/FOUNDER_NORMALIZATION_FIX.md`.) This M-step is also renormalized
per-founder by default (§4.3).

> **Default normalization changed (2026-07-06): the M-step now divides each
> founder's numerator by its own effective k-mer content $K^w_f$ before the L1
> renormalization** (`normalize="per_founder"`, the new default), not by a single
> global scalar. The boxed update above is the legacy `global` (multinomial)
> form; it is retained via `--normalize global` for reproducing old runs but is
> deprecated. The per-founder form and its derivation are in §4.3, with the full
> "why" writeup in `docs/FOUNDER_NORMALIZATION_FIX.md`.

Vectorized implementation (`em_solver.py:80-110`):

```python
wc      = counts if omega is None else (omega * counts)  # ω_k · c_k
denom   = np.maximum(h @ kmer_pa, 1e-7)   # μ_k(h), with numerical floor
cw      = wc / denom                  # ω_k · c_k / μ_k
em_term = h * (kmer_pa @ cw)               # h_f · Σ_k kmer_pa[f,k] · ω_k · c_k / μ_k
h_new   = em_term / em_term.sum()     # L1 renormalize
```

**Initialization (`em_solver.py:70-71`):** $h_f^{(0)} = 1/F$ (uniform).
**Termination (`em_solver.py:108-111`):** halt when
$\lVert\mathbf{h}^{(t+1)} - \mathbf{h}^{(t)}\rVert_2 < \varepsilon$, with
$\varepsilon = 10^{-7}$ in production (`per_sample_per_chrom.py:225`,
`block_em.py:362`). Typical 30–80 iterations in float32.

**EM is run only on k-mers with $c_k > 0$** (`per_sample_per_chrom.py:207-211`).
The Poisson M-step sums $(K_{\mathrm{pa}})_{f,k} \, c_k / \mu_k$, so terms with $c_k = 0$
contribute zero — filtering is exact, not an approximation. Production
$n_{nz} / K \approx 0.05$–$0.30$ depending on coverage.

The update is the standard Lee–Seung multiplicative NMF update with simplex
renormalization; the same iteration arises in finite-mixture EM (Dempster,
Laird, Rubin 1977) and in RNA-seq abundance estimation (RSEM, kallisto,
salmon) under the substitution {transcripts → founders, reads → k-mers}.

### 4.1 Optional priors (per-window mode only, see §7)

Two regularizers are wired into `em_solver.solve_em` and used in `block_em.py`
when running per-window EM:

**Symmetric Dirichlet$(\alpha)$ prior on $\mathbf{h}$** (`em_solver.py:74-75, 92-94`):

$$h_f^{(t+1)} \;\propto\; h_f^{(t)} \sum_k (K_{\mathrm{pa}})_{f,k} \frac{c_k}{\mu_k} \;+\; (\alpha - 1)$$

$\alpha = 1$ recovers the MLE; $\alpha > 1$ pulls toward uniform.

**Anchor toward a prior $\mathbf{h}_{\text{prior}}$** (typically the
chromosome-wide global $\hat{\mathbf{h}}$, weighted by $\beta$;
`em_solver.py:78-94`):

$$h_f^{(t+1)} \;\propto\; h_f^{(t)} \sum_k (K_{\mathrm{pa}})_{f,k} \frac{c_k}{\mu_k} \;+\; \beta \, N \, h_{\text{prior},f}, \qquad N = \sum_k c_k$$

$\beta = 0$ is pure MLE; $\beta = 1$ weights the prior as much as the data;
$\beta \in [0.05, 0.5]$ is the working range.

### 4.2 Per-bubble de-replication weight $\omega_k = 1/m_b$ (window-mode weighting)

> **Scope (superseded 2026-07-06 for global mode).** $\omega_k = 1/m_b$ is the
> **window-mode** EM weighting. It was originally the global-mode production
> weighting too, but with the per-founder M-step normalization of §4.3 the
> cactus/PG completeness imbalance this weight patched is removed at its source,
> so global-mode production now uses `--kmer-weight uniform` + `--normalize
> per_founder` (ablation: uniform AF-MAE 0.0036 vs $1/m_b$ 0.0045). The
> derivation below stands as the rationale for the window-mode weight and as
> archaeology; see §4.3 / `docs/FOUNDER_NORMALIZATION_FIX.md`.

**Definition.** Every panel k-mer $k$ is contributed by exactly one PanGenie
bubble (the local pangenome graph object that produced it; `build_kmer_pa.py`
preserves the bubble assignment in `meta["bubble_id"]`). Let

$$m_b(k) \;=\; \#\{k' : \mathrm{bubble}(k') = \mathrm{bubble}(k)\}$$

be the number of canonical k-mers in $k$'s bubble (a per-bubble integer, in
[1, ~32] given PanGenie's per-allele caps). The production weighting is

$$\boxed{\;\omega_k \;=\; \frac{1}{m_b(k)}\;}$$

so a bubble contributing $m_b$ k-mers contributes total mass $m_b \cdot (1/m_b)
= 1$ to the M-step — one *locus* worth of evidence, independent of how many
k-mers it produced.

**Implementation** (`per_sample_per_chrom.py:121-128`, `:181-187`): computed at
load time from `meta["bubble_id"]` as
`m_b = np.bincount(bid)[bid]; omega = 1.0 / m_b` and passed into
`em_solver.solve_em(omega=...)`. Selected via the CLI flag
`--kmer-weight inv_mb` (default `uniform`, i.e. $\omega_k \equiv 1$). Wired
identically in the global and window modes.

**Why this helps the 231 panel.** The mixed panel (78 long-read cactus + 153
PG-genotyped short-read founders) has a structural imbalance: cactus founders
carry roughly 2× as many private k-mers *per founder* as PG founders (median
8.7K vs 4.5K; §8.1), and they disproportionately populate large bubbles
(multi-allelic / repetitive regions). *(This per-founder count ratio is a
different metric from the ~16× cactus/PG **modality-private** ratio measured in
fully-called Chr1 bubbles in §2.1/M3 — both quantify the same assembly-vs-short-read
asymmetry, at per-founder vs per-bubble-singleton granularity respectively.)* Under unweighted EM, "h proportional to k-mer count" then becomes "h
proportional to k-mer-rich-region count", over-crediting cactus founders by
$\sim 41\%$ in balanced ground truth (`docs/METHODS_TRIED_AND_RESULTS.md` §2).
$\omega_k = 1/m_b$ turns the estimand into "h proportional to *locus* count,"
which is panel-symmetric and cancels the over-credit at its source.

**Empirical result.** On the heterogeneous 231 panel, `filt2 + ω_k=1/m_b
global` wins per-record AF MAE in every cov/n/g regime tested (g0 sims,
replicate-validated; `docs/METHODS_TRIED_AND_RESULTS.md` §3). $\omega_k = 1/m_b$ is
therefore the **window-mode** EM weighting alongside the filt2 kmer_pa; for
global mode it is superseded by per_founder + uniform (§4.3).

**Equivalence to count rescaling.** Because the Poisson M-step is linear in
$c_k$, $\omega_k$ enters the update exactly as if each count had been rescaled
$c_k \leftarrow \omega_k c_k$ pre-EM. Implementation matches: `solve_em` does
`wc = omega * counts` once and runs the unweighted iteration on `wc`
(`em_solver.py:80-81`). This is *not* a coverage rescaling — $\omega_k$ varies
across k-mers and is decoupled from $\lambda$.

**Panel-conditional caveat (do not propagate into the paper).** On a
homogeneous all-long-read control panel (`benchmarks/p80`, 2026-05-27) without
the cactus/PG imbalance, $\omega_k = 1/m_b$ slightly *under-performs* plain
$\omega_k=1$ (per-record AF MAE +2% to +41% across regimes). $\omega_k = 1/m_b$
is therefore an imbalance-canceling correction rather than a universal
de-replication. For the GrENE-Net 231 panel this distinction is internal to
the project (see §10 M6); the paper treats $\omega_k = 1/m_b$ as the method's
default. The flag `--kmer-weight uniform` reproduces the unweighted MLE
byte-identically (`em_solver.py:80`, `:65`) for users on balanced panels.

### 4.3 Per-founder normalization $\omega$/$K^w_f$ (production default, `normalize="per_founder"`)

**The bug it fixes.** The M-step of §4 renormalizes each iteration by a single
**global** scalar $\mathrm{total\_c}=\sum_k \omega_k c_k$ — the same denominator
for every founder. At the fixed point the per-founder numerator is
$\mathrm{raw}_f = h^{\text{true}}_f\cdot K^w_f$, where

$$K^w_f \;=\; \sum_{k} \omega_k\,(K_{\mathrm{pa}})_{f,k}$$

is founder $f$'s $\omega$-weighted k-mer content **over the full estimation unit
(all k-mers, since $E[c_k/\mu_k]=1$ for every carried k-mer whether or not it
draws $c_k=0$)**. Normalizing by one global constant therefore converges to
$\hat h_f \propto h^{\text{true}}_f\cdot K^w_f$: k-mer-rich founders (the 78
cactus founders) are over-credited and k-mer-poor founders (many of the 153 PG
founders) are starved, and the multiplicative update compounds the starvation to
numerical zero. On the equimolar seed-mix this drove ~19–49 of 231 founders to
$\sim10^{-15}$–$10^{-31}$ and biased every $p_0$-anchored quantity. It is a
**completeness bias** (present even on noiseless counts), *not* read noise: the
per-founder error correlates $+0.84$ with $K^w_f$ under the global form.

**The per-founder M-step.** Divide each founder's numerator by its own $K^w_f$
before renormalizing to the simplex:

$$\boxed{\; B_f \;=\; \mathrm{total\_c}\cdot\frac{\mathrm{raw}_f / K^w_f}{\sum_{f'} \mathrm{raw}_{f'}/K^w_{f'}} \;}$$

then apply the Dirichlet / anchor pseudocounts of §4.1 and renormalize as before.
This makes $\mathbf{h}^{\text{true}}$ an **exact fixed point for any $K^w$
heterogeneity**. It is the RNA-seq **effective-length** correction: RSEM
(Li & Dewey 2011) forms $\tau_i=(\theta_i/\ell_i)/\sum_j(\theta_j/\ell_j)$ under
{transcripts $\to$ founders, reads $\to$ k-mers, effective length
$\ell\to K^w$}; kallisto/salmon carry the same $/\ell$ term the old kMate EM
omitted.

**Implementation notes.**
- The Dirichlet $(\alpha-1)$ and anchor $\texttt{prior\_weight}\cdot\mathrm{total\_c}\cdot h_{\text{prior}}$ pseudocounts stay **outside** the $/K^w$ division, so the anchor pull is governed by `prior_weight` independent of $K^w$.
- $K^w_f$ is computed over the **full** estimation unit (ALL k-mers, incl. $c_k=0$) with a $10^{-12}$ floor. This is the unbiased normalizer ($E[\mathrm{raw}_f]=h_f K^w_{f,\text{full}}$); summing over the **observed** ($c_k>0$) columns instead is a **survivorship bias** that shrinks $K^w$ for k-mer-poor/discriminative founders drawing bad-luck zeros and reintroduces the collapse (on a noisy $\sim$0.3× seed-mix: $\sim$24 absorbed vs 0). Callers that pre-slice `kmer_pa` to observed columns MUST pass `kfw` computed over the full unit (the global driver's `kfw_full`, `block_em`'s full-window `kfw`, the bootstrap's `kfw_full`); when the full matrix is passed directly the in-solver `kmer_pa @ w` fallback is already correct.
- When $K^w$ is constant across founders the division is a common factor that cancels in the L1 renormalization, so `per_founder` reduces **byte-for-byte** to the legacy `global` update.

**Wiring / status.** `normalize="per_founder"` is the **default** in
`em_solver.solve_em`, `block_em.solve_em_per_block`, and the
`per_sample_per_chrom.py` driver, selected via `--normalize
{per_founder,global}`; the legacy multinomial normalization (`--normalize global`)
is deprecated. **Global mode is the GrENE-Net production estimator** — the cohort
is heavily selfing, so there is little recombination mosaic and one founder
mixture per chromosome is appropriate for *both* the seed-mix $p_0$ and the
evolved samples. It emits **both** the per-sample founder $\mathbf h$ and the
per-record SNP/SV `alt_freq` (the driver always projects $\hat{\mathbf h}$ through
`var_pa`) in one pass; the production recipe is per_founder + filt2inv +
`--kmer-weight uniform` (**drop** $\omega=1/m_b$; §4.2), superseding the
$\omega=1/m_b$-for-global recommendation. **Window mode** (per-window EM + anchor +
HMM smoothing) is for *recombinant* pools, not the selfing production path;
per_founder is its single default too (it wins the local-only use case by ~18–21%
RMSE; the difference under the full `star2` recipe is a marginal bias-variance
effect). Full rationale, factorial ablation, and seed-mix / hapFIRE comparison:
`docs/FOUNDER_NORMALIZATION_FIX.md`.

### 4.4 Haploblock collapse — fit the distinct haplotypes, not all $F$ (default 2026-07-07)

Before the EM runs on a **unit** (a window in window mode, or the whole
chromosome in global mode), kMate **computes** how many distinct haplotypes the
panel actually resolves over that unit's k-mers, and fits the EM over those
$K_b \le F$ haplotypes rather than *assuming* all $F$ founders are separately
identifiable. Founders whose k-mer presence pattern over the unit is identical
(exact, `--haploblock-eps 0`, the default) are one haplotype; the EM fits a
$K_b$-vector $\hat{\mathbf h}^{\,c}$ and each class frequency is split **equally**
back to its member founders,
$\hat h_f = \hat h^{\,c}_{\text{class}(f)} / |\text{class}(f)|$.
This is the hapFIRE-like *block → haploblock → EM* design and removes the
flat-likelihood ridge among k-mer-indistinguishable founders (which is what
drove the multiplicative EM to spurious sparse vertices — cf. `HANDOFF`
non-identifiability). Implemented in `em_solver.haploblock_collapse_indices`
(packbits presence signature → grouping) and applied in both
`block_em.solve_em_per_block` and `per_sample_per_chrom.run_one_chrom_global`;
**on by default** (`haploblock_collapse=True`).

- **$K_b = F$ is an EXACT (byte-identical) no-op.** When every founder is
  distinct over the unit, the code fits the founders directly in native order (no
  permutation, no split), so it is bit-for-bit the pre-collapse result. On the
  231-founder arch3 panel at r²=0.1+ / whole-chromosome units, $K_b$ is
  essentially always 230–231 (the accessions are *not* k-mer-identical over
  regions large enough to carry many k-mers), so **for the GrENE-Net global
  production path the collapse is numerically inert** — it is the *framework*
  (derive haplotypes from data, don't assume) that matters, and it is what makes
  window mode well-posed where a 10 kb window may carry only a handful of
  haplotypes.
- **Uncertainty is computed on the same $K_b$ parametrization.** For `--emit-af-se`
  and the parametric bootstrap (`h_uncertainty.bootstrap_cov_h`), the covariance
  is estimated on the $K_b$ haplotype classes and mapped to founder space by the
  deterministic equal-split $\mathbf h = A\,\mathbf h^{\,c}$ (so
  $\Sigma_{\mathbf h} = A\,\Sigma_c\,A^{\top}$), never on the singular full-$F$
  Fisher information across identical class members. (Again a no-op at $K_b=F$.)
- **`--haploblock-eps` > 0** additionally merges founders whose presence differs
  in $\le \varepsilon\cdot(\text{unit k-mers})$ — an *approximate*, order-dependent,
  representative-based clustering for merging near-indistinguishable founders. The
  production default is $\varepsilon = 0$ (exact, lossless).

---

## 5. Per-chromosome execution (production)

The production driver `per_sample_per_chrom.py:68-121` runs **one EM per
chromosome**, not a single genome-wide EM. Rationale (`per_sample_per_chrom.py:1-13`):
loading the full genome-wide $K_{\mathrm{pa}}$ matrix would peak at ~74 GB float32; per-chrom
peaks ~5× lower, fitting in 32–64 GB SLURM allocations.

The per-chromosome $\hat{\mathbf{h}}_c$ vectors agree with the genome-wide
$\hat{\mathbf{h}}$ to within ~0.1% in practice (because each chromosome's
~16M k-mers already massively over-determines the 231-vector simplex).

For the methods section: state this as **one EM per chromosome** and treat the
genome-wide formulation as the conceptual model.

---

## 6. Allele-frequency projection (corrected formula)

Given the per-chromosome $\hat{\mathbf{h}}_c$, the alternate-allele frequency
at record $r$ on chromosome $c$ is the **missing-aware projection**
(`per_sample_per_chrom.py:129-150`):

$$\boxed{\; \widehat{\mathrm{AF}}_r \;=\; \frac{\hat{\mathbf{h}}_c^{\!\top} \, (V_{\mathrm{pa}})_{:,r}}{\hat{\mathbf{h}}_c^{\!\top} \, (V_{\mathrm{called}})_{:,r}} \;}$$

with a safe-divide floor of $10^{-12}$ on the denominator. This is the standard
"AF among called samples" (AC/AN). For uniform $\mathbf{h} = (1/F) \mathbf{1}$ it
collapses to AC/AN exactly.

**Why the call-mask matters.** Under the bare projection $\hat{\mathbf{h}}^\top
V_{\mathrm{pa}}$ (the formula in the old `ALGORITHM.md` and
`CACTUS_EM_MATH.md`), `./.` cells are treated as REF (because $V_{\mathrm{pa}}$ is 0
both for confirmed-REF and for missing). At records where a substantial
fraction of founders is `./.` (e.g. cactus-only or PG-only records after the
cactus-78 + PG-153 merge, where `F_MISSING` runs up to 0.66), this
systematically under-counts the true AF. The fix divides by the h-weighted
called mass at each record. See `build_var_pa.py:12-18` for the production
note.

### 6.1 Per-record uncertainty outputs (added 2026-05-21)

Alongside `alt_freq`, every record carries three uncertainty / QC fields, written
for both `global` and `window` modes (`per_sample_per_chrom.py:658-684`):

| field | definition | h-dependent? | code |
|---|---|---|---|
| `info` | $\hat{\mathbf{h}}_c^{\!\top}\,(V_{\mathrm{called}})_{:,r}\in[0,1]$ — the h-weighted called mass (the projection denominator) | yes | `:246-249` |
| `n_called` | $\sum_f (V_{\mathrm{called}})_{f,r}\in\{0,\dots,F\}$ — integer count of called founders at $r$ | no | `:594-600` |
| `se` | $\sqrt{\,p(1-p)/\max(n_{\text{called}},1)\,}$, $p=\widehat{\mathrm{AF}}_r$ clipped to $[0,1]$; NaN if AF is NaN or $n_{\text{called}}=0$ | only through $p$ | `:662-666` |

**What `se` is.** The Wald (binomial) standard error of a proportion
$p=\widehat{\mathrm{AF}}_r$ whose effective sample size is the number of *called
founders* at that record, $n_{\text{called}}$. For uniform $\mathbf{h}$,
$\widehat{\mathrm{AF}}_r=\text{AC}/\text{AN}$ and `se` is exactly the textbook SE
of that AC/AN proportion. So `se` answers **"how well is the per-record AF pinned
down by the panel genotypes available at that record?"** — a panel-completeness /
AC-AN sampling term. Records where most founders are `./.` (small
$n_{\text{called}}$) get a large `se`; fully-genotyped records get a small one.

**What `se` does NOT capture — carry into the paper and downstream.** `se` is a
*panel-side* term only. It deliberately omits the three variance sources that
actually dominate a pool-seq AF estimate:

1. **EM estimation error in $\hat{\mathbf{h}}$.** Finite reads make
   $\hat{\mathbf{h}}$ itself uncertain; that uncertainty propagates linearly
   through the projection but is invisible to `se`.
2. **Read-coverage (Poisson) sampling.** At low coverage the k-mer counts — and
   hence $\hat{\mathbf{h}}$ — are noisier, but `se` does not move with coverage.
3. **Pool finite-$N$ sampling.** The realized pool is a multinomial draw of $N$
   individuals from the source population (`docs/SIMULATIONS_METHODS.md` §2); the
   source-vs-realized gap is irreducible and not in `se`.

`se` also keys off the integer *count* $n_{\text{called}}$, not the h-weighted
mass `info`: a record whose called founders all carry negligible $\hat{h}$ can
still report a small `se` when $n_{\text{called}}$ is large. Treat `info` and
`n_called` as complementary — `info` is the *observed* confidence, `n_called` the
*panel* confidence.

**Downstream guidance.** Use `info` / `n_called` / `se` to **filter or flag**
low-panel-support records (e.g. drop records below an `n_called` threshold, as the
p80 notebook does at `missing_frac ≤ 10%`). Do **not** inverse-variance-weight on
`se` as if it were the full estimation error — that weights by panel completeness,
not AF precision. A precision estimate that folds in the EM/coverage terms would
require propagating $\hat{\mathbf{h}}$ covariance or a read/k-mer bootstrap, which
is **not implemented** in the production driver.

---

## 7. Per-window mode (recombinant pools)

For pools whose ancestry is a recombinant mosaic (evolved samples after
multiple generations of meiosis), partition the chromosome into windows
$\{W_w\}$ and run the EM independently within each window. Window-mode
projection becomes

$$\widehat{\mathrm{AF}}_r \;=\; \frac{\hat{\mathbf{h}}_{w(r)}^{\!\top} \, (V_{\mathrm{pa}})_{:,r}}{\hat{\mathbf{h}}_{w(r)}^{\!\top} \, (V_{\mathrm{called}})_{:,r}}$$

where $w(r)$ is the window containing record $r$.

**One estimator, selected by UNIT (`--unit`, 2026-07-07).** With haploblock
collapse + per-founder normalization + local-only, "global" and "window" are the
same algorithm at different unit sizes — `partition → per-unit (haploblock
collapse §4.4 → local EM) → project`. The unit:
- `--unit chrom` (**default**): one $\hat{\mathbf{h}}_c$ per chromosome — the
  production estimator for selfers / inbreds / F0 pools (SEEDMIX, GrENE-Net);
  robust on uniform and sparse panels; the only unit supporting `--h-only` and
  `--emit-af-se`.
- `--unit ld` (`--ld-r2 0.1`): r²-LD blocks from the panel's own
  `var_pa` via CompleteLDPartition (`ld_partition.py`; blocks are a panel property,
  computed once and cached; ~18 blocks/Chr1). A per-block option that **collapses
  in low-diversity blocks** (centromere) → not appropriate for selfing pools.
- `--unit bp` (`--window-bp`): fixed-bp windows (hard cuts, production 10 kb).
- `--unit tsv` (`--blocks-tsv`): an explicit block partition.

`--block-mode global|window` remain as **deprecated aliases**
(`global`→`--unit chrom`, `window`→`--unit bp`); `--unit chrom` is byte-identical
to the former `--block-mode global`. The earlier `ld_gabriel`, `ld_complete`,
`bigld_panel` modes and the overlapping-window variant were archived to
`src/archive/`.

**Window mode defaults to LOCAL-ONLY (2026-07-07):** `--local-only` is the
window-mode default. It is GLOBAL-FREE — **no** global anchor prior, **no**
fallback to the chromosome-wide $\hat{\mathbf h}_c$ (thin/empty windows and
unassigned records → NaN AF), and **no** cross-window smoothing. Rationale: a
recombinant window carries only a few haplotypes (§4.4), so the chromosome-wide
mixture is the wrong prior for it, and each window is fit purely on its own
k-mers over the $K_b$ haplotypes it actually contains. $\hat{\mathbf h}_c$ is
still computed and stored for diagnostics only. This supersedes the earlier
anchored+smoothed "star2" default.

**Legacy "star2" recipe** (`--no-local-only`): per-window EM anchored to the
chromosome-wide $\hat{\mathbf{h}}_c$ via the Dirichlet anchor of §4.1
(`--global-anchor-weight`, star2 value 0.3) and smoothed across windows (§8,
`--hmm-smooth-passes` 5). Available for reproducing old runs. Caveat: under a
genuine within-window haploblock collapse (§4.4) the anchor is summed per class
and split *equally* across k-mer-indistinguishable members, so it cannot preserve
intra-class prior asymmetry — one more reason local-only is the default.

**Window-mode caveat.** A 100 kb hard-window mode on the 827k-SNP panel
(hapFIRE-equivalent fine-scale) blew up haplotype-class count to 196/231 →
HARP-style class methods failed; the founder-simplex EM in kmer_pa-window mode
was unaffected (it doesn't form discrete classes). See
`PIPELINE_STATE_2026-05-19.md`.

---

## 8. EM variants — production status (updated 2026-07-08)

Current production:
- **M-step normalization `per_founder` (§4.3)** is the production default
  everywhere (`--normalize per_founder`), removing the founder-completeness bias
  at its source.
- **Per-bubble de-replication weight $\omega_k = 1/m_b$ (§4.2)** is the
  **window-mode** EM weighting (`--kmer-weight inv_mb`). For **global mode** (the
  GrENE-Net production estimator) it was superseded 2026-07-06 by
  `--kmer-weight uniform` + per_founder (§4.3). The `kmer_pa_231_arch3_filt2inv`
  filtered matrix (§2.1) is the production kmer_pa base in both modes.
- **Haploblock collapse (§4.4)** — fit the $K_b\le F$ distinct haplotypes per
  unit, not all $F$ — is **on by default** in both modes (`--haploblock-eps 0`,
  exact; an exact no-op at $K_b=F$, i.e. inert on the 231-founder global path).
- **Window mode defaults to `--local-only` (§7, 2026-07-07)**: no global anchor,
  no global fallback, no smoothing. Global anchor (`--global-anchor-weight`) and
  HMM smoothing (`--hmm-smooth-*`, §8.3) are now only reached via the legacy
  `--no-local-only` "star2" recipe (anchor 0.3, passes 5, α 0.5), not the default.

Removed / archived (do not use):
- **K-mer-budget balancing (§8.1)** — earlier row-normalization +
  $K_f^\alpha$ post-correction rebalancing approach. Superseded by $\omega_k =
  1/m_b$, which targets the same imbalance via the cleaner composite-likelihood
  route (§4.2). Archived to `src/archive/`.
- **Carrier-weighted counts (§8.2)** and **contamination-ω (§8.4)** —
  archived; never beat plain EM and are not part of the paper's algorithm.

The subsections below remain as archaeology for the archived variants; the
active production-EM weighting is documented in §4.2, not in §8.

### 8.1 K-mer-budget balancing (`--row-normalize-kmer_pa`, `--kf-correction-alpha`)

Empirically, the 78 cactus founders carry roughly 2× as many private k-mers
per founder as the 153 PanGenie-imputed founders (median 8.7K vs 4.5K),
producing a +41% h-bias toward cactus founders even in a balanced ground
truth. *(Per-founder count ratio; distinct from the ~16× modality-private
ratio in §2.1/M3, which counts cactus- vs PG-private singletons within
fully-called bubbles.)* The mitigation
(`per_sample_per_chrom.py:111-121, 153-183`):

1. Pre-EM: $(K_{\mathrm{pa}})_{f,k} \leftarrow (K_{\mathrm{pa}})_{f,k} / K_f$ where
   $K_f = \sum_k (K_{\mathrm{pa}})_{f,k}$. This makes each founder row sum to 1.
2. Post-EM: $h_f^{\text{proj}} \propto \hat{h}_f / K_f^{\alpha}$ (renormalize).
   $\alpha = 1$ is the mass-domain inverse exactly; $\alpha > 1$ over-corrects
   to suppress residual h-bias.

Tested in parallel SLURM arrays (k-mer-imbalance investigation,
2026-05-15); **none of the rebalancing variants beat the plain Poisson EM
(`filt2`) on per-record AF MAE**. The class h-bias and per-record AF MAE are
anti-correlated; balancing one degrades the other. Plain EM was the
production default at the time; both this approach and plain EM have since
been **superseded (2026-05-27) by the weighted EM of §4.2** ($\omega_k =
1/m_b$), which targets the same imbalance via a composite-likelihood form
that wins per-record AF MAE on the 231 panel.

### 8.2 Carrier-weighted counts (`--ac-weight-counts`)

Multiplies $c_k \leftarrow c_k \cdot a_k$ where $a_k = \sum_f (K_{\mathrm{pa}})_{f,k}$
before the EM (`per_sample_per_chrom.py:214-222`). Cancels the
"singleton voice" amplification: in plain EM, a k-mer carried by a single
founder gets a $1/a_k = 1$ weight in the M-step, vs. $1/231$ for a
universally-shared k-mer, which over-amplifies private k-mers. Off by default.

### 8.3 HMM smoothing across blocks (`--hmm-smooth-passes`)

Post-EM Li-Stephens-style smoothing of $\hat{\mathbf{h}}_w$ across adjacent
windows, with recombination-rate-weighted neighbor averaging
(`per_sample_per_chrom.py:375-391`; `block_haplotype_em.smooth_h_across_blocks`).
The **production window defaults** are 5 passes, $\alpha = 0.5$ (the `★★` recipe;
§7, §8, and `docs/PIPELINE_STATE.md` §0.1). An earlier `clean_smooth` sweep found
10 passes / $\alpha = 0.2$ / recomb rate $4 \times 10^{-8}$/bp best on that test;
the production defaults (5 / 0.5) supersede it.

### 8.4 Contamination ω (`em_solver.solve_em_with_omega`)

Augments the EM with a per-k-mer contamination rate $\omega_k$ that doesn't
depend on $\mathbf{h}$ (`em_solver.py:118-175`). Experimental — not invoked by
the production driver.

---

## 9. Per-sample output schema

`per_sample_per_chrom.py:672-684` writes a TSV with columns:

| column | meaning |
|---|---|
| `chrom` | chromosome (string, e.g. `Chr1`) |
| `pos` | 1-based position (matches `var_pa.meta.pos`; freqk users: see §10 L4) |
| `ref_len` | length of REF allele (bp) |
| `alt_len` | length of ALT allele (bp) |
| `alt_freq` | $\widehat{\mathrm{AF}}_r$, the missing-aware projection of §6 |
| `info` | h-weighted called mass per record, $\in [0,1]$ |
| `n_called` | integer panel-level called count, $\in [0,F]$ |
| `se` | Wald binomial SE assuming $n = n_{\text{called}}$ |

Also written: `*.h_per_chrom.npz` (global mode) or `*.h_blocks_per_chrom.npz`
(window mode) with the per-chrom or per-block h vectors.

---

## 10. Audit findings (carry into the paper / discussion)

The following are issues / assumptions worth knowing about. Severity tags:
`CRITICAL` = wrong scientific result, `HIGH` = wrong in realistic edge cases,
`MEDIUM` = fragile / undocumented, `LOW` = minor.

| ID | Sev | Issue |
|---|---|---|
| C1 | CRITICAL | Old AF projection formula in `ALGORITHM.md` and `CACTUS_EM_MATH.md` (`AF = h^T V_pa`) silently treats `./.` as REF and under-counts AF on records with high `F_MISSING`. **Fixed in production**: §6 above is the correct formula. The old MDs are why this single source of truth was needed. |
| H1 | HIGH | Production runs **one EM per chromosome**, not a single genome-wide EM. Empirically per-chrom $\hat{\mathbf{h}}_c$ agrees to ~0.1% with genome-wide $\hat{\mathbf{h}}$, but methods text must say "per-chromosome" not "genome-wide". |
| H2 | HIGH | The EM filters to $c_k > 0$ k-mers before solving. Mathematically equivalent under the M-step, but the matrix actually fed to EM is ~5–30% the size of $K$ depending on coverage. State this so reviewers don't get confused by "80M k-mer EM" claims. |
| H3 | ~~HIGH~~ → resolved | The kmer_pa-vs-var_pa heterozygote-handling asymmetry in code (`GT[0]`-only vs `any(a > 0)`) is **not active** for the production panel: the arch3 `merged_231_chr{N}_final.vcf.gz` is pre-haploidized (only `.`, `0`, `1` tokens; verified 2026-05-22). pysam returns 1-tuples, both definitions agree. **Latent hazard closed 2026-05-22**: both `build_kmer_pa.py:189-202` and `build_var_pa.py:62-73` now raise `ValueError` on any GT with `len(gt) != 1`, so feeding a diploid VCF fails fast instead of silently producing inconsistent matrices. |
| H4 | HIGH | The per-record `se` (`per_sample_per_chrom.py:662-666`) is a **panel-side** binomial SE: $\sqrt{p(1-p)/n_{\text{called}}}$ with effective N = number of genotyped founders. It captures panel completeness / AC-AN sampling only, and **omits** EM-$\hat{\mathbf{h}}$ estimation error, read-coverage Poisson noise, and pool finite-$N$ sampling. Safe to use for **filtering/flagging** low-support records; **not** safe to inverse-variance-weight on as if it were the full AF precision (that weights by panel completeness). Documented in §6.1; a coverage/EM-aware uncertainty is not implemented. |
| M1 | MEDIUM | The old MD coverage formula "$\lambda$ = total read length / genome size" doesn't match production code, which uses $\hat\lambda = F \cdot \sum c_k / \sum_{f,k} (K_{\mathrm{pa}})_{f,k}$ (uniform-h estimate). Cosmetic since $\lambda$ cancels; fixed in §3 above. |
| M2 | ~~MEDIUM~~ → resolved | Production tolerance is $\varepsilon = 10^{-7}$ throughout (`per_sample_per_chrom.py:225`, `block_em.py:362`); §4 above uses that single value. |
| M3 | MEDIUM | `--treat-missing-as-n` in `build_kmer_pa.py`. **Corrected 2026-05-29:** the then-production `v3qc_v3` (and `v3qc_v2`) builds — and the current arch3 `kmer_pa` — were/are built with this **ON** (`./.` → N, every k-mer over the span dropped), *not* OFF as this table previously stated. Verified from build scripts + build logs; no N-off build exists on disk. N-on is the deliberate choice — `./.`→REF would fabricate confident reference genotypes from low-quality no-calls. **Hypothesis tested and REJECTED 2026-05-29:** we suspected N-on *causes* the cactus/PG private-k-mer imbalance via PG SV missingness (dropped PG k-mers → cactus-skewed survivors). The missingness-causation test (`notebooks/MISSINGNESS_CAUSES_IMBALANCE.ipynb`, `scripts/run_missingness_test.py`) shows otherwise: in *fully-called* Chr1 bubbles (where N-on ≡ N-off, so missingness cannot contribute) the private ratio is **16.5×**, *higher* than the all-bubble 15.5×, and flat across PG-missingness strata (16.3× at 0% → 15.3× at >40%); cactus privates are no more enriched in missing bubbles (81.9%) than PG privates (82.7%). The imbalance is **real biology** — long-read cactus assemblies realize ~16× more private k-mers/founder than short-read PG-genotyped founders. **Consequences:** (1) the proposed per-founder `kmer_pa` call-mask / dosage fix is **shelved** — it would not reduce the imbalance, which is what `filt2`/ω=1/m_b (§2.1, §4.2) already target; (2) N-on remains the correct modeling choice on its own merits; (3) `kmer_pa` is still not missing-aware the way the §6 projection is (`var_called`) — a latent correctness nuance, but not the imbalance driver. For the paper: state the panel is N-on and that the cactus/PG private skew is a genotyping-modality (assembly vs short-read) effect, not a missingness artifact. |
| M4 | MEDIUM | `denom = max(h @ kmer_pa, 1e-7)` numerical floor in EM (`em_solver.py:87`). Kicks in only at simplex boundaries; mention if you want to be precise. |
| M5 | MEDIUM | `samtools fastq -F 0x900` filters secondary+supplementary but **not** PCR duplicates (`kmer_count.py:72`). SEEDMIX is PCR-free (memory: `seedmix_is_pcr_free`) so no dedup needed. **Evolved GrENE-Net samples are not PCR-free and require an upstream dedup step** (clumpify or equivalent) before kMate; otherwise PCR-duplicate reads inflate k-mer counts and bias the EM. State the preprocessing distinction in the paper's per-sample pipeline section. |
| M6 | MEDIUM | Weighting $\omega_k = 1/m_b$ (§4.2) is **panel-conditional**: it wins per-record AF MAE on the heterogeneous 231 panel (where it cancels the cactus/PG imbalance) but slightly under-performs $\omega_k=1$ on the homogeneous 80-cactus control panel (+2% to +41% MAE, `benchmarks/p80/results/filt2_mb_vs_uniform_summary.tsv`). **Superseded 2026-07-06 for global mode:** with the per-founder M-step normalization (§4.3) the completeness imbalance is removed at its source, so global-mode production uses `--kmer-weight uniform` + `--normalize per_founder` (uniform AF-MAE 0.0036 vs $1/m_b$ 0.0045; the old multinomial + $1/m_b$ was the worst factorial row); $\omega_k=1/m_b$ remains the **window-mode** weighting. Code preserves both via `--kmer-weight {uniform,inv_mb}`. See `docs/FOUNDER_NORMALIZATION_FIX.md` and `docs/METHODS_TRIED_AND_RESULTS.md` §3. |
| L1 | LOW | `kmer_pa` is conceptually genome-wide but physically per-chromosome (`build_kmer_pa.build_kmer_pa_for_chrom`). Cosmetic. |
| L2 | LOW | `solve_em_with_omega` (contamination) exists in code but is unused in production runs. |
| L3 | LOW | `freqk` reports `VCF_pos - 2`; add 2 before joining freqk output to any VCF or var_pa (memory: `freqk_pos_offset`). |
| L4 | LOW | `pos` in the kMate TSV is the original VCF position (1-based). Joining to other tools requires checking their offset convention. |

**Top 3 to carry into the paper:**

1. **C1** — state the AF projection with the missing-mask denominator, not the
   bare projection. This is the single biggest correctness item.
2. **H1** — say "per-chromosome EM" not "genome-wide EM" in the methods.
3. **M5** — state that the per-sample preprocessing for evolved GrENE-Net
   samples includes PCR-duplicate removal (SEEDMIX is PCR-free and skips this
   step).

---

## 11. Relationship to existing methods

| | kMate | HARP / hapFIRE | freqk |
|---|---|---|---|
| Latent | founder mixture $\mathbf{h}$ on simplex | per-LD-block haplotype mixture → CVXPY founder solve | none (no founder concept) |
| Evidence | canonical k-mer counts $c_k$ | per-base read pileup $P(\text{base}|\text{hap}, q)$ | per-bubble k-mer counts (independent bubbles) |
| Inference | weighted Poisson EM, multiplicative update (per-founder M-step norm §4.3; $\omega_k=1/m_b$ window / uniform global, §4.2) | per-base likelihood EM (HARP) | per-bubble alt/total ratio |
| Cross-bubble pooling | yes — bubbles with sparse k-mers borrow strength via h | yes within an LD block | no |
| Projection | $\hat{\mathbf{h}}^\top V_{\mathrm{pa}} / \hat{\mathbf{h}}^\top V_{\mathrm{called}}$ | $\hat{\mathbf{h}}_{\text{hap}}^\top \mathrm{H2S}$ + founder solve | per-bubble alt-count / total-count |
| Sees SVs | yes — k-mers span bubble alleles | no (per-base SNP-only) | yes (bubble-level) |
| Reference | pangenome graph | single linear reference | pangenome graph |

The M-step form
$h_f^{(t+1)} \propto h_f^{(t)} \sum_k (K_{\mathrm{pa}})_{f,k}\, \omega_k\, c_k / \mu_k$
also recurs in NMF (Lee & Seung 1999), the original EM derivation (Dempster
et al. 1977), and RNA-seq abundance estimation (RSEM/kallisto/salmon); the
$\omega_k$ weighting is the composite-likelihood generalization (Lindsay 1988).

---

## References

- Dempster A, Laird N, Rubin D (1977) *J. R. Stat. Soc. B* 39:1–38 — EM algorithm.
- Lee D, Seung S (1999) *Nature* 401:788–791 — NMF multiplicative update.
- Marçais & Kingsford (2011) — Jellyfish.
- Li & Dewey (2011) — RSEM. Bray et al. (2016) — kallisto. Patro et al. (2017) — salmon.
- Kessner et al. (2013) — HARP.
- Wu et al. (2024) — hapFIRE.
- Ebler et al. (2022) *Nat. Genet.* 54:518–525 — PanGenie.
- Danecek et al. (2021) — bcftools.
- Hickey et al. (2024) — minigraph-cactus.
- Lindsay B (1988) — Composite likelihoods. *Contemp. Math.* 80:221–239.

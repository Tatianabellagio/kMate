# kMate — algorithm, math, and production wiring (code-verified)

Single source of truth for the kMate algorithm. Replaces `ALGORITHM.md` (out-of-date
prose) and `CACTUS_EM_MATH.md` (out-of-date math). Every formula and pipeline step
below has been cross-checked against the production code; references are
`file:line` so they stay traceable.

Code-verified against:
- `src/em_solver.py` — EM core
- `src/build_kmer_cn.py` — cn_full builder
- `src/build_cn_var.py` — cn_var + cn_var_called builder
- `src/per_sample_per_chrom.py` — production per-sample driver
- `src/block_em.py` — per-window / per-block EM
- `src/kmer_count.py` — jellyfish wrapper

Last verified: 2026-05-27.

> **Note (2026-05-26 cleanup):** the estimator was reduced to two modes,
> `global` and `window`; LD-block modes, overlapping windows, the older k-mer
> rebalancing (`--row-normalize-cn`, see §8.1), carrier-weighting (§8.2) and
> contamination-ω (§8.4) were archived to `src/archive/`. §7–§8 below
> are updated, but the inline `file:line` references elsewhere predate the
> cleanup and may have shifted. The authoritative file list and invocation
> recipes now live in `src/README.md`.
>
> **Note (2026-05-27 production decision):** the cn_full singleton filter
> `filt2` (drop k-mers with column-sum < 2; §2.1) and the per-k-mer EM weight
> $\omega_k = 1/m_b$ (per-bubble de-replication; §4.2) are now production
> defaults on the heterogeneous 231-founder panel. Both are documented inline
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

Three sparse $F \times \cdot$ indicator matrices are built once per panel:

### 2.1 `cn` — founder × k-mer (the EM evidence)

$$\mathrm{cn} \in \{0,1\}^{F \times K}, \quad \mathrm{cn}_{f,k} = \mathbb{1}[\text{founder } f \text{ carries panel k-mer } k]$$

Built per-chromosome by `build_kmer_cn.py:96-230` from PanGenie's bubble-level
k-mer index (`kmers.tsv.gz`), the panel VCF, and the reference FASTA. For each
PanGenie bubble:

1. PanGenie's k-mer list for the bubble is taken as-is. PanGenie has already
   applied: (i) in-bubble dedup, (ii) genomic ref-dedup, (iii) a cap of 16
   k-mers/allele for biallelic bubbles and 32/allele for multi-allelic bubbles.
   We inherit those filters.
2. For each founder, we reconstruct their haplotype across the bubble by
   applying their genotype to the reference (`reconstruct_haplotype`,
   `build_kmer_cn.py:45-76`).
3. We compute the canonical k-mer set of the reconstructed haplotype
   (`canonical_kmer_set`, `build_kmer_cn.py:34-42`). A canonical k-mer is the
   lexicographically smaller of $(k, \text{revcomp}(k))$; k-mers containing N
   are skipped.
4. $\mathrm{cn}_{f,k} = 1$ iff the bubble's k-mer $k$ is in founder $f$'s
   canonical k-mer set.

**Haploid panel input.** The panel VCF that feeds both `cn` and `cn_var` is
*pre-haploidized* (e.g. `panel/pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz`):
each GT field contains a single allele, `0`, `1`, ..., or `.`. PanGenie outputs
diploid GTs for the 153 short-read–genotyped founders, which we haploidize
before merging with the 78 long-read–assembled cactus founders, so the resulting
231-founder panel is entirely haploid. Under pysam this surfaces as
single-element tuples (`(0,)`, `(1,)`, `(None,)`); `build_kmer_cn.py:189-194`
reads `gt[0]`, which is the only allele.

**Missing-GT handling (`build_kmer_cn.py:189-194`):** controlled by
`--treat-missing-as-n`. Production `cn_full_v3qc_v3` was built with
`treat_missing_as_n=False`, i.e. `./.` → REF in the reconstructed haplotype.
The earlier `v3qc_v2` build with `treat_missing_as_n=True` over-inflated
private-k-mer ratios via N-poisoning of short bubbles and is archived.

**Production singleton filter — `filt2` (`build_cn_full_filt2_v3qc_v3_chr1.sh`).**
After building $\mathrm{cn}$, the production matrix drops every k-mer column
$k$ with allele-count $a_k = \sum_f \mathrm{cn}_{f,k} = 1$ (so-called
*singletons* — k-mers carried by a single founder). The retained matrix is

$$\mathrm{cn}^{\text{filt2}} = \mathrm{cn}_{:,\,\{k\,:\,a_k \geq 2\}}.$$

Rationale: singletons over-represent founder-private sequencing errors and
private repeats that survive PanGenie's in-bubble dedup, and they contribute
zero discrimination across founders. The threshold sweep (drop $a_k < 2$, $3$,
$5$, $10$) found that **filt2 is the optimum**: filt3/5/10 over-prune real
signal (`docs/METHODS_TRIED_AND_RESULTS.md` §1). The production cn_full referenced
elsewhere in this document is `cn_full_231_v3qc_v3_filt2`; all `cn` references
in §3–§4 implicitly use the filt2-filtered matrix.

### 2.2 `cn_var` — founder × variant (the projection target)

$$\mathrm{cn}_{\text{var}} \in \{0,1\}^{F \times R}, \quad \mathrm{cn}_{\text{var},f,r} = \mathbb{1}[\text{any allele of founder } f \text{'s GT at } r \text{ is alt}]$$

Built by `build_cn_var.py:30-109` from the same haploid biallelic-decomposed
panel VCF used for `cn`, processed with `bcftools norm -m -`. A founder is
marked as an alt-carrier iff any allele in its GT field is non-zero
(`build_cn_var.py:68`: `any(a is not None and a > 0 for a in gt)`). On a
haploid GT (`(0,)`, `(1,)`, ...) this is equivalent to checking `gt[0] > 0`,
so the alt-carrier definitions in `cn` and `cn_var` agree exactly under the
production panel.

### 2.3 `cn_var_called` — founder × variant (the call mask)

$$\mathrm{cn}_{\text{var\_called}} \in \{0,1\}^{F \times R}, \quad \mathrm{cn}_{\text{var\_called},f,r} = \mathbb{1}[\text{founder } f \text{'s GT at } r \text{ is NOT } ./.]$$

Built alongside `cn_var` (`build_cn_var.py:62-66`). This call mask is **critical
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

$$\boxed{\; c_k \;\sim\; \mathrm{Poisson}\!\left(\lambda \cdot \mu_k(\mathbf{h})\right), \qquad \mu_k(\mathbf{h}) = \sum_{f=1}^F h_f \, \mathrm{cn}_{f,k} \;}$$

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
unweighted MLE corresponds to $\omega_k \equiv 1$. The production weighting
$\omega_k = 1/m_b(k)$ — per-bubble de-replication — is derived and motivated
in §4.2.

**Coverage estimate $\lambda$.** Reported at load time as a sanity-check value
(`per_sample_per_chrom.py:104-106`):

$$\hat{\lambda} \;=\; \frac{F \cdot \sum_k c_k}{\sum_{f,k} \mathrm{cn}_{f,k}}$$

(total counts $\times$ $F$ divided by total founder–k-mer carrier entries;
assumes uniform $h$). **$\lambda$ cancels in the M-step ratio (§4) and is not
consumed by the EM** — the value is logged for inspection only.

---

## 4. EM algorithm

Treat each unit of count at k-mer $k$ as a latent draw from one of the
founders carrying $k$. Let $z_{k,f}$ be the (unobserved) count attributed to
founder $f$, with $z_{k,f} = 0$ if $\mathrm{cn}_{f,k} = 0$ and
$\sum_f z_{k,f} = c_k$.

**E-step.** Given current $\mathbf{h}^{(t)}$:

$$\mathbb{E}\!\left[z_{k,f} \mid c_k, \mathbf{h}^{(t)}\right] \;=\; c_k \cdot \frac{h_f^{(t)} \, \mathrm{cn}_{f,k}}{\mu_k(\mathbf{h}^{(t)})}$$

**M-step (weighted form, production).** Constrained MLE on the simplex under
the weighted likelihood of §3:

$$\boxed{\; h_f^{(t+1)} \;\propto\; h_f^{(t)} \cdot \sum_{k=1}^K \mathrm{cn}_{f,k} \cdot \frac{\omega_k\, c_k}{\mu_k(\mathbf{h}^{(t)})} \;}$$

followed by L1 renormalization $\sum_f h_f^{(t+1)} = 1$. The coverage $\lambda$
cancels. Setting $\omega_k \equiv 1$ recovers the unweighted Poisson MLE
update; the production setting $\omega_k = 1/m_b(k)$ is the per-bubble
de-replication weight (§4.2).

Vectorized implementation (`em_solver.py:80-110`):

```python
wc      = counts if omega is None else (omega * counts)  # ω_k · c_k
denom   = np.maximum(h @ cn, 1e-7)   # μ_k(h), with numerical floor
cw      = wc / denom                  # ω_k · c_k / μ_k
em_term = h * (cn @ cw)               # h_f · Σ_k cn[f,k] · ω_k · c_k / μ_k
h_new   = em_term / em_term.sum()     # L1 renormalize
```

**Initialization (`em_solver.py:70-71`):** $h_f^{(0)} = 1/F$ (uniform).
**Termination (`em_solver.py:108-111`):** halt when
$\lVert\mathbf{h}^{(t+1)} - \mathbf{h}^{(t)}\rVert_2 < \varepsilon$, with
$\varepsilon = 10^{-7}$ in production (`per_sample_per_chrom.py:225`,
`block_em.py:362`). Typical 30–80 iterations in float32.

**EM is run only on k-mers with $c_k > 0$** (`per_sample_per_chrom.py:207-211`).
The Poisson M-step sums $\mathrm{cn}_{f,k} \, c_k / \mu_k$, so terms with $c_k = 0$
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

$$h_f^{(t+1)} \;\propto\; h_f^{(t)} \sum_k \mathrm{cn}_{f,k} \frac{c_k}{\mu_k} \;+\; (\alpha - 1)$$

$\alpha = 1$ recovers the MLE; $\alpha > 1$ pulls toward uniform.

**Anchor toward a prior $\mathbf{h}_{\text{prior}}$** (typically the
chromosome-wide global $\hat{\mathbf{h}}$, weighted by $\beta$;
`em_solver.py:78-94`):

$$h_f^{(t+1)} \;\propto\; h_f^{(t)} \sum_k \mathrm{cn}_{f,k} \frac{c_k}{\mu_k} \;+\; \beta \, N \, h_{\text{prior},f}, \qquad N = \sum_k c_k$$

$\beta = 0$ is pure MLE; $\beta = 1$ weights the prior as much as the data;
$\beta \in [0.05, 0.5]$ is the working range.

### 4.2 Per-bubble de-replication weight $\omega_k = 1/m_b$ (production)

**Definition.** Every panel k-mer $k$ is contributed by exactly one PanGenie
bubble (the local pangenome graph object that produced it; `build_kmer_cn.py`
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
carry roughly 2× as many private k-mers per founder as PG founders, and they
disproportionately populate large bubbles (multi-allelic / repetitive
regions). Under unweighted EM, "h proportional to k-mer count" then becomes "h
proportional to k-mer-rich-region count", over-crediting cactus founders by
$\sim 41\%$ in balanced ground truth (`docs/METHODS_TRIED_AND_RESULTS.md` §2).
$\omega_k = 1/m_b$ turns the estimand into "h proportional to *locus* count,"
which is panel-symmetric and cancels the over-credit at its source.

**Empirical result.** On the heterogeneous 231 panel, `filt2 + ω_k=1/m_b
global` wins per-record AF MAE in every cov/n/g regime tested (g0 sims,
replicate-validated; `docs/METHODS_TRIED_AND_RESULTS.md` §3). $\omega_k = 1/m_b$ is
therefore the production EM weighting alongside the filt2 cn_full.

**Equivalence to count rescaling.** Because the Poisson M-step is linear in
$c_k$, $\omega_k$ enters the update exactly as if each count had been rescaled
$c_k \leftarrow \omega_k c_k$ pre-EM. Implementation matches: `solve_em` does
`wc = omega * counts` once and runs the unweighted iteration on `wc`
(`em_solver.py:80-81`). This is *not* a coverage rescaling — $\omega_k$ varies
across k-mers and is decoupled from $\lambda$.

**Panel-conditional caveat (do not propagate into the paper).** On a
homogeneous all-long-read control panel (`control_p80`, 2026-05-27) without
the cactus/PG imbalance, $\omega_k = 1/m_b$ slightly *under-performs* plain
$\omega_k=1$ (per-record AF MAE +2% to +41% across regimes). $\omega_k = 1/m_b$
is therefore an imbalance-canceling correction rather than a universal
de-replication. For the GrENE-Net 231 panel this distinction is internal to
the project (see §10 M6); the paper treats $\omega_k = 1/m_b$ as the method's
default. The flag `--kmer-weight uniform` reproduces the unweighted MLE
byte-identically (`em_solver.py:80`, `:65`) for users on balanced panels.

---

## 5. Per-chromosome execution (production)

The production driver `per_sample_per_chrom.py:68-121` runs **one EM per
chromosome**, not a single genome-wide EM. Rationale (`per_sample_per_chrom.py:1-13`):
loading the full genome-wide `cn` matrix would peak at ~74 GB float32; per-chrom
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

$$\boxed{\; \widehat{\mathrm{AF}}_r \;=\; \frac{\hat{\mathbf{h}}_c^{\!\top} \, \mathrm{cn}_{\text{var},:,r}}{\hat{\mathbf{h}}_c^{\!\top} \, \mathrm{cn}_{\text{var\_called},:,r}} \;}$$

with a safe-divide floor of $10^{-12}$ on the denominator. This is the standard
"AF among called samples" (AC/AN). For uniform $\mathbf{h} = (1/F) \mathbf{1}$ it
collapses to AC/AN exactly.

**Why the call-mask matters.** Under the bare projection $\hat{\mathbf{h}}^\top
\mathrm{cn}_{\text{var}}$ (the formula in the old `ALGORITHM.md` and
`CACTUS_EM_MATH.md`), `./.` cells are treated as REF (because `cn_var` is 0
both for confirmed-REF and for missing). At records where a substantial
fraction of founders is `./.` (e.g. cactus-only or PG-only records after the
cactus-78 + PG-153 merge, where `F_MISSING` runs up to 0.66), this
systematically under-counts the true AF. The fix divides by the h-weighted
called mass at each record. See `build_cn_var.py:12-18` for the production
note.

### 6.1 Per-record uncertainty outputs (added 2026-05-21)

Alongside `alt_freq`, every record carries three uncertainty / QC fields, written
for both `global` and `window` modes (`per_sample_per_chrom.py:658-684`):

| field | definition | h-dependent? | code |
|---|---|---|---|
| `info` | $\hat{\mathbf{h}}_c^{\!\top}\,\mathrm{cn}_{\text{var\_called},:,r}\in[0,1]$ — the h-weighted called mass (the projection denominator) | yes | `:246-249` |
| `n_called` | $\sum_f \mathrm{cn}_{\text{var\_called},f,r}\in\{0,\dots,F\}$ — integer count of called founders at $r$ | no | `:594-600` |
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

$$\widehat{\mathrm{AF}}_r \;=\; \frac{\hat{\mathbf{h}}_{w(r)}^{\!\top} \, \mathrm{cn}_{\text{var},:,r}}{\hat{\mathbf{h}}_{w(r)}^{\!\top} \, \mathrm{cn}_{\text{var\_called},:,r}}$$

where $w(r)$ is the window containing record $r$.

**Partitioning schemes** (selected via `--block-mode`; two modes only after the
2026-05-26 cleanup):
- `global` (default): one $\hat{\mathbf{h}}_c$ per chromosome — for selfers /
  inbreds / F0 pools (SEEDMIX).
- `window` (fixed-bp): hard cuts every `--window-bp` (production default 10 kb).
  Per-window EM is anchored to the chromosome-wide $\hat{\mathbf{h}}_c$ and
  smoothed across windows (§8); each record takes the $\hat{\mathbf{h}}$ of its
  window (hard assignment).

The earlier `ld_gabriel`, `ld_complete`, `bigld_panel` modes and the
overlapping-window variant were archived to `src/archive/`.

**Anchoring** (`--global-anchor-weight`): per-window EM can be anchored to
the chromosome-wide $\hat{\mathbf{h}}_c$ via the Dirichlet anchor of §4.1 with
$\beta$ = `--global-anchor-weight`. Used when low-evidence windows would
otherwise collapse onto a single founder.

**Window-mode caveat.** A 100 kb hard-window mode on the 827k-SNP panel
(hapFIRE-equivalent fine-scale) blew up haplotype-class count to 196/231 →
HARP-style class methods failed; the founder-simplex EM in cn-window mode
was unaffected (it doesn't form discrete classes). See
`PIPELINE_STATE_2026-05-19.md`.

---

## 8. EM variants — production status (updated 2026-05-27)

Current production:
- **Per-bubble de-replication weight $\omega_k = 1/m_b$ (§4.2)** is the
  production EM weighting on the 231 heterogeneous panel, selected via
  `--kmer-weight inv_mb`. The `cn_full_231_v3qc_v3_filt2` singleton-filtered
  matrix (§2.1) is the production cn_full base.
- **Global anchor (`--global-anchor-weight`) and HMM smoothing
  (`--hmm-smooth-*`, §8.3)** are part of the production *window* recipe — they
  are the window-mode defaults (anchor 0.3, passes 5, α 0.5; see §7), not
  experimental.

Removed / archived (do not use):
- **K-mer-budget balancing (§8.1)** — earlier row-normalization +
  $K_f^\alpha$ post-correction rebalancing approach. Superseded by $\omega_k =
  1/m_b$, which targets the same imbalance via the cleaner composite-likelihood
  route (§4.2). Archived to `src/archive/`.
- **Carrier-weighted counts (§8.2)** and **contamination-ω (§8.4)** —
  archived; never beat plain EM and are not part of the paper's algorithm.

The subsections below remain as archaeology for the archived variants; the
active production-EM weighting is documented in §4.2, not in §8.

### 8.1 K-mer-budget balancing (`--row-normalize-cn`, `--kf-correction-alpha`)

Empirically, the 78 cactus founders carry roughly 2× as many private k-mers
per founder as the 153 PanGenie-imputed founders (median 8.7K vs 4.5K),
producing a +41% h-bias toward cactus founders even in a balanced ground
truth. The mitigation
(`per_sample_per_chrom.py:111-121, 153-183`):

1. Pre-EM: $\mathrm{cn}_{f,k} \leftarrow \mathrm{cn}_{f,k} / K_f$ where
   $K_f = \sum_k \mathrm{cn}_{f,k}$. This makes each founder row sum to 1.
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

Multiplies $c_k \leftarrow c_k \cdot a_k$ where $a_k = \sum_f \mathrm{cn}_{f,k}$
before the EM (`per_sample_per_chrom.py:214-222`). Cancels the
"singleton voice" amplification: in plain EM, a k-mer carried by a single
founder gets a $1/a_k = 1$ weight in the M-step, vs. $1/231$ for a
universally-shared k-mer, which over-amplifies private k-mers. Off by default.

### 8.3 HMM smoothing across blocks (`--hmm-smooth-passes`)

Post-EM Li-Stephens-style smoothing of $\hat{\mathbf{h}}_w$ across adjacent
windows, with recombination-rate-weighted neighbor averaging
(`per_sample_per_chrom.py:375-391`; `block_haplotype_em.smooth_h_across_blocks`).
Best params from `clean_smooth`: 10 passes, $\alpha = 0.2$, recomb rate
$4 \times 10^{-8}$/bp.

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
| `pos` | 1-based position (matches `cn_var.meta.pos`; freqk users: see §10 L4) |
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
| C1 | CRITICAL | Old AF projection formula in `ALGORITHM.md` and `CACTUS_EM_MATH.md` (`AF = h^T cn_var`) silently treats `./.` as REF and under-counts AF on records with high `F_MISSING`. **Fixed in production**: §6 above is the correct formula. The old MDs are why this single source of truth was needed. |
| H1 | HIGH | Production runs **one EM per chromosome**, not a single genome-wide EM. Empirically per-chrom $\hat{\mathbf{h}}_c$ agrees to ~0.1% with genome-wide $\hat{\mathbf{h}}$, but methods text must say "per-chromosome" not "genome-wide". |
| H2 | HIGH | The EM filters to $c_k > 0$ k-mers before solving. Mathematically equivalent under the M-step, but the matrix actually fed to EM is ~5–30% the size of $K$ depending on coverage. State this so reviewers don't get confused by "80M k-mer EM" claims. |
| H3 | ~~HIGH~~ → resolved | The cn-vs-cn_var heterozygote-handling asymmetry in code (`GT[0]`-only vs `any(a > 0)`) is **not active** for the production panel: `founders_231_v3qc_v3.haploid.vcf.gz` is pre-haploidized (only `.`, `0`, `1` tokens; verified 2026-05-22). pysam returns 1-tuples, both definitions agree. **Latent hazard closed 2026-05-22**: both `build_kmer_cn.py:189-202` and `build_cn_var.py:62-73` now raise `ValueError` on any GT with `len(gt) != 1`, so feeding a diploid VCF fails fast instead of silently producing inconsistent matrices. |
| H4 | HIGH | The per-record `se` (`per_sample_per_chrom.py:662-666`) is a **panel-side** binomial SE: $\sqrt{p(1-p)/n_{\text{called}}}$ with effective N = number of genotyped founders. It captures panel completeness / AC-AN sampling only, and **omits** EM-$\hat{\mathbf{h}}$ estimation error, read-coverage Poisson noise, and pool finite-$N$ sampling. Safe to use for **filtering/flagging** low-support records; **not** safe to inverse-variance-weight on as if it were the full AF precision (that weights by panel completeness). Documented in §6.1; a coverage/EM-aware uncertainty is not implemented. |
| M1 | MEDIUM | The old MD coverage formula "$\lambda$ = total read length / genome size" doesn't match production code, which uses $\hat\lambda = F \cdot \sum c_k / \sum_{f,k} \mathrm{cn}_{f,k}$ (uniform-h estimate). Cosmetic since $\lambda$ cancels; fixed in §3 above. |
| M2 | ~~MEDIUM~~ → resolved | Production tolerance is $\varepsilon = 10^{-7}$ throughout (`per_sample_per_chrom.py:225`, `block_em.py:362`); §4 above uses that single value. |
| M3 | MEDIUM | `--treat-missing-as-n` flag in `build_kmer_cn.py`. Production v3qc-v3 was built with this **OFF** (./. → REF). The v3qc-v2 build with the flag ON inflated K_f ratios and is archived. State which build the paper uses. |
| M4 | MEDIUM | `denom = max(h @ cn, 1e-7)` numerical floor in EM (`em_solver.py:87`). Kicks in only at simplex boundaries; mention if you want to be precise. |
| M5 | MEDIUM | `samtools fastq -F 0x900` filters secondary+supplementary but **not** PCR duplicates (`kmer_count.py:72`). SEEDMIX is PCR-free (memory: `seedmix_is_pcr_free`) so no dedup needed. **Evolved GrENE-Net samples are not PCR-free and require an upstream dedup step** (clumpify or equivalent) before kMate; otherwise PCR-duplicate reads inflate k-mer counts and bias the EM. State the preprocessing distinction in the paper's per-sample pipeline section. |
| M6 | MEDIUM | Production weighting $\omega_k = 1/m_b$ (§4.2) is **panel-conditional**: it wins per-record AF MAE on the heterogeneous 231 panel (where it cancels the cactus/PG imbalance) but slightly under-performs $\omega_k=1$ on the homogeneous 80-cactus control panel (+2% to +41% MAE, `control_p80/results/filt2_mb_vs_uniform_summary.tsv`). For the GrENE-Net 231 application the win is decisive; the paper treats $\omega_k=1/m_b$ as the method default. Code preserves both via `--kmer-weight {uniform,inv_mb}` so balanced-panel users can opt out. See `docs/METHODS_TRIED_AND_RESULTS.md` §3. |
| L1 | LOW | `cn_full` is conceptually genome-wide but physically per-chromosome (`build_kmer_cn.build_cn_for_chrom`). Cosmetic. |
| L2 | LOW | `solve_em_with_omega` (contamination) exists in code but is unused in production runs. |
| L3 | LOW | `freqk` reports `VCF_pos - 2`; add 2 before joining freqk output to any VCF or cn_var (memory: `freqk_pos_offset`). |
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
| Inference | weighted Poisson EM, multiplicative update ($\omega_k=1/m_b$, §4.2) | per-base likelihood EM (HARP) | per-bubble alt/total ratio |
| Cross-bubble pooling | yes — bubbles with sparse k-mers borrow strength via h | yes within an LD block | no |
| Projection | $\hat{\mathbf{h}}^\top \mathrm{cn}_{\text{var}} / \hat{\mathbf{h}}^\top \mathrm{cn}_{\text{var\_called}}$ | $\hat{\mathbf{h}}_{\text{hap}}^\top \mathrm{H2S}$ + founder solve | per-bubble alt-count / total-count |
| Sees SVs | yes — k-mers span bubble alleles | no (per-base SNP-only) | yes (bubble-level) |
| Reference | pangenome graph | single linear reference | pangenome graph |

The M-step form
$h_f^{(t+1)} \propto h_f^{(t)} \sum_k \mathrm{cn}_{f,k}\, \omega_k\, c_k / \mu_k$
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

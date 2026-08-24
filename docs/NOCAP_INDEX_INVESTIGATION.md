# Do we need PanGenie's per-bubble k-mer caps? (no-caps investigation)

**Date:** 2026-05-31 · **Status:** index-level result in; EM-accuracy test pending
**Verdict (current):** **Probably NOT necessary in `inv_mb` (window) mode.** The
caps' statistical role is already played — more precisely — by our `inv_mb` EM
weight. Removing them is a compute-vs-marginal-variance tradeoff, not a
correctness fix. Default stays capped pending an empirical accuracy check.

> **Update 2026-07-06 (global-mode production weighting changed — see
> `docs/FOUNDER_NORMALIZATION_FIX.md`, ALGORITHM.md §4.3).** Global mode (the
> GrENE-Net production estimator) now uses `--kmer-weight uniform` +
> `--normalize per_founder`, superseding the `inv_mb`-for-global default this
> doc assumed. That flips the premise below: under `uniform` there is **no**
> per-bubble normalization, so by §3's own logic the caps **do** matter for
> global-mode production. Net effect: the "keep caps" recommendation is
> **reinforced** for global mode; the "`inv_mb` subsumes the caps" argument now
> applies to **window mode** only.

---

## 1. What the caps are

In-house index builder `panel/pangenie_index/scripts/build_kmers_tsv.py` reproduces
PanGenie-index's per-bubble unique-k-mer selection, including three caps (PanGenie
defaults, faithfully copied):

| Cap | Value | What it limits |
|---|---|---|
| `--cap-biallelic` | 16 | unique k-mers kept **per allele** in a biallelic bubble |
| `--cap-multiallelic` | 32 | unique k-mers kept **per allele** in a multi-allelic bubble |
| per-bubble total | `max(nr_paths, 301)` | total k-mers kept across the **whole bubble** |
| (`--overhang-cap`) | 12 | flanking-region k-mers per side — *separate, not lifted by `--no-caps`* |

Source of truth in PanGenie: `external/pangenie-tools/src/stepwiseuniquekmercomputer.cpp`
`select_kmers()` (lines 46–93) — `max_kmers = is_biallelic ? 16 : 32`,
`max_alleles = max(nr_of_paths, 301)`. **These are bare constants with no
explanatory comment.** A new `--no-caps` flag (added 2026-05-31) lifts all three
in-bubble caps (overhang untouched).

## 2. Why does PanGenie cap? (not in the paper — inferred from the model)

The PanGenie paper (Ebler et al. 2022, *Nat Genet* — `papers/s41588-022-01043-w.pdf`)
describes how unique k-mers are *determined* (Methods, p8: "a set of k-mers for each
bubble that occurs at most once within a single allele sequence and are not found
anywhere outside the variant bubble") but **never states why the count is capped**.
The rationale is evident from PanGenie's model, not the text:

1. **Statistical (pseudo-replication).** PanGenie genotypes via an HMM whose
   per-bubble **emission probability multiplies a Poisson term per selected
   k-mer**, treating them as conditionally independent. But unique k-mers within
   one short allele are **highly correlated** — the same reads cover them. Adding
   many correlated k-mers artificially sharpens the emission likelihood
   (over-counts evidence → over-confident genotype calls). Capping per allele
   bounds this over-counting.
2. **Computational.** The cap bounds index size, memory, and HMM runtime. This
   matters most at **complex STR/VNTR bubbles** — the paper notes (p9) "most
   complex bubbles are located inside STR/VNTR regions", and those generate huge
   numbers of unique k-mers. `max(nr_paths, 301)` is the safety valve.

## 3. Why kMate (probably) does NOT need them

kMate is **not** PanGenie's HMM. It runs a per-founder **weighted Poisson EM** on
the k-mer-count simplex (`src/em_solver.py`), then projects through `var_pa`.
Crucially, production uses the **`inv_mb` weight** `ω_k = 1/m_b` (ALGORITHM.md
§4.2), where `m_b` = number of k-mers in k's bubble:

> "a bubble contributing `m_b` k-mers contributes total mass `m_b · (1/m_b) = 1`
> to the M-step — **one *locus* worth of evidence, independent of how many
> k-mers**." (composite likelihood, Lindsay 1988)

That is **exactly the statistical job PanGenie's caps were doing**, done more
precisely and continuously:

- **Under `inv_mb` (production):** a bubble's total contribution is normalized to
  1 *regardless of k-mer count*. So lifting the caps is **bubble-weight-neutral** —
  a bubble with 7,000 k-mers contributes the same mass as one with 16. The extra
  k-mers can only refine the *within-bubble average* (sampling-variance reduction);
  they **cannot** let a high-k-mer bubble dominate `h`. ⇒ caps are not needed for
  correctness; expected effect on accuracy is neutral-to-slightly-positive.
- **Under `uniform` (`ω_k ≡ 1`, the balanced-panel MLE mode, §10 M6):** there is
  **no** per-bubble normalization, so removing caps WOULD reintroduce
  pseudo-replication / bubble domination. **For uniform weighting the caps DO
  matter.** No-caps is only safe because production is `inv_mb`.

The remaining cost of no-caps is therefore **purely computational**, not statistical.

## 4. Index-level result (Chr1 canary)

`scripts/nocap_chr1_canary.sh` → `analysis/panel_qc/nocap_vs_capped_kmer_pa_chr1.txt`
(comparator `scripts/compare_kmer_pa_dirs.py`). No-cap index → K_pa, both `filt2inv`:

| | k-mers | nnz | per-founder ratio |
|---|---|---|---|
| capped (production) | 10,950,663 | 461,069,206 | 1.00× |
| **nocap** | **18,486,760** | 654,387,165 | **1.42×** |
| Δ | **+68.8%** | +41.9% | +42% k-mers/founder |

- **capped ⊆ nocap, exactly** (`capped-only = 0`; 100% subset). Caps are a pure
  truncation — they never pick a *different* k-mer, only fewer.
- **Carriers identical on shared k-mers:** 100.0000% per-cell, 100% identical
  columns. Removing caps does not perturb *which* founder carries *which* k-mer.
- Raw no-cap = 50.3M k-mers; `filt2inv` drops ac==1 (18.1M private singletons),
  ac==0 (13.4M), ac==F (0.33M), keeping 36.7%. Most of the explosion is private
  noise the filter removes — but **+7.5M genuinely informative (`2≤ac≤230`)
  k-mers** were being discarded by the caps.
- Build: ~5 h wall, 53 GB peak (savio4_htc condo).

## 5. Verdict & what would change it

**Removing caps is not necessary for production (inv_mb).** The caps' statistical
purpose is subsumed by `ω_k = 1/m_b`; what's left is a ~69%-larger index for, at
best, within-bubble variance reduction that the per-bubble weight largely caps
anyway. Keep the PanGenie-default caps as production default.

**The one thing that would flip this:** an empirical accuracy test. Run a benchmark
sample (SEEDMIX or the p80 control) through the EM with the no-cap Chr1 K_pa vs the
capped baseline (both `--kmer-weight inv_mb`) and compare AF RMSE / outliers. If
no-caps shows a real, replicated accuracy gain that justifies the index-size and
runtime cost, revisit; otherwise the caps stay. **Until that test is run, the
recommendation is: keep caps.**

Artifacts (Chr1 only; gitignored): index `panel/pangenie_index/pang_135_haploid_nocap/`,
K_pa `data/kmer_pa_231_arch3_nocap_filt2inv/` — archived 2026-07-07 to
`panel/pangenie_index/archive/kmer_pa_231_arch3_nocap_filt2inv/` (still regenerable at the
original `data/` path via `scripts/nocap_chr1_canary.sh` if the accuracy test above is ever run).

## 6. See also
- `ALGORITHM.md` §4.2 (the `inv_mb` de-replication weight) and §10 M6 (uniform vs inv_mb)
- `docs/PIPELINE_STATE.md` §0 (production index = in-house, capped; PG = comparator)
- `external/pangenie-tools/src/stepwiseuniquekmercomputer.cpp:46` (the caps)

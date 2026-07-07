# Founder-normalization fix — the per-founder M-step (`normalize="per_founder"`)

**Date:** 2026-07-06
**Scope:** Why kMate's EM M-step now normalizes each founder by its own
**full-panel** k-mer content (`normalize="per_founder"`, the new default) instead
of a single global scalar (`normalize="global"`, the legacy multinomial behaviour).
**Status:** Decided and wired. `per_founder` is the default everywhere; `global`
is kept for reproducing old runs but is deprecated.

Code-verified against `src/em_solver.py`, `src/block_em.py`,
`src/per_sample_per_chrom.py`, and the hapFIRE comparator
`external/HapFIRE/haplotype_generation.py`.

---

## TL;DR

- **Symptom.** In global mode the EM drove ~19–49 of the 231 equimolar seed-mix
  founders to numerically ~0 (worst $\sim 10^{-15}$ to $10^{-31}$), even though
  every founder is present at $1/231$ by design. This biased every
  $p_0$-anchored quantity — the founder-level GEA read spurious selection
  ($s$ = slope-from-$p_0$) on founders the estimator had absorbed. Per-record
  **allele-frequency accuracy was fine** ($r\approx 0.98$); only the founder
  **decomposition** $\mathbf{h}$ was wrong.
- **Root cause.** The M-step normalized each iteration by one *global* scalar,
  the same for every founder. At the fixed point this makes
  $\hat h_f \propto h^{\text{true}}_f\cdot K^w_f$, where $K^w_f$ is founder $f$'s
  own $\omega$-weighted **full-panel** k-mer content. Founders with more
  k-mers (the long-read/cactus founders) are over-credited; k-mer-poor founders
  (many short-read/PanGenie founders) are starved, and the multiplicative EM
  compounds the starvation to zero. This is a **completeness bias**, present even
  on noiseless synthetic counts.
- **Fix.** Divide each founder's M-step numerator by its own **full-panel** $K^w_f$
  before renormalizing to the simplex — the RNA-seq **effective-length**
  correction (RSEM/kallisto/salmon). This makes $\mathbf{h}^{\text{true}}$ an
  *exact* fixed point for any $K^w$ heterogeneity, and reduces byte-for-byte to
  the old `global` behaviour when $K^w$ is constant across founders. **$K^w_f$
  must be summed over the full estimation unit (ALL k-mers, including those with
  $c_k=0$ this run), not over the observed set** — see §3(b).
- **Do NOT add an artificial frequency floor.** hapFIRE avoids the bias
  structurally (many-block averaging), not via a floor; the principled fix is
  removing the bias (`per_founder`) plus averaging over units (per-chromosome /
  per-block).

---

## 1. Problem

kMate global-mode EM collapsed a large fraction of the 231 founders to
numerically zero on the equimolar seed-mix. All founders are present at exactly
$1/231$ by construction, so any founder driven to $\sim 10^{-15}$–$10^{-31}$ is a
pure estimation artifact. Because the seed-mix decomposition defines $p_0$ (the
founder frequencies at generation 0), and the founder-level GEA reads a selection
coefficient as the slope of evolved-vs-$p_0$ frequency, a collapsed founder
manifests as **spurious selection**. The whole founder-level GEA was affected.

Crucially the failure is confined to the **decomposition** $\mathbf{h}$: the
downstream per-record allele-frequency projection (§6 of `ALGORITHM.md`) was
still accurate ($r\approx 0.98$ against truth). It is the founder identities, not
the AF, that were wrong.

---

## 2. Root cause (code-verified)

The EM M-step (`src/em_solver.py`) normalized each iteration by a **global scalar**
`total_c` (the sum of all $\omega$-weighted counts) — identical for every founder.
Write the M-step numerator (before normalization) as

$$\mathrm{raw}_f \;=\; h_f \cdot \sum_{k} (K_{\mathrm{pa}})_{f,k}\,\frac{\omega_k\,c_k}{\mu_k(\mathbf{h})}.$$

Since $E[c_k/\mu_k]=1$ for every carried k-mer (whether or not it happens to draw
$c_k=0$), the expected numerator is $h_f$ times founder $f$'s own **full-panel
k-mer budget**,

$$K^w_f \;=\; \sum_{k} \omega_k\,(K_{\mathrm{pa}})_{f,k},$$

the $\omega$-weighted count of **all** k-mers that founder $f$ carries in the
panel/window. Dividing by one *global* constant therefore converges to

$$\hat h_f \;\propto\; h^{\text{true}}_f \cdot K^w_f.$$

Founders carrying more k-mers are systematically over-credited; k-mer-poor
founders are starved. On the 231 panel this maps directly onto the assembly
modality: the 78 long-read/cactus founders are k-mer-rich and get inflated, while
many of the 153 short-read/PanGenie founders are k-mer-poor and get starved. The
multiplicative EM update ($h^{(t+1)}\propto h^{(t)}\cdot\ldots$) compounds a
small per-iteration starvation into numerical zero.

This is a **completeness bias**, not a read-noise or coverage artifact: it
appears even on **noiseless synthetic counts**. The smoking gun — the per-founder
error correlates **+0.84 with $K^w_f$** under the old normalization.

---

## 3. The fix — `normalize="per_founder"` (now the default)

Divide each founder's M-step numerator by its own full-panel $K^w_f$, then
renormalize to the simplex:

$$B_f \;=\; \mathrm{total\_c}\cdot\frac{\mathrm{raw}_f / K^w_f}{\sum_{f'} \mathrm{raw}_{f'}/K^w_{f'}},$$

after which the existing Dirichlet / anchor pseudocounts are applied and the
vector is renormalized exactly as before. Because the bias was
$\hat h_f\propto h^{\text{true}}_f\cdot K^w_f$, dividing by $K^w_f$ makes
$\mathbf{h}^{\text{true}}$ an **exact fixed point for any $K^w$ heterogeneity**.

**This is the RNA-seq effective-length correction.** RSEM (Li & Dewey 2011,
verified against the paper) computes transcript abundance as
$\tau_i = (\theta_i/\ell_i)\big/\sum_j(\theta_j/\ell_j)$ — dividing each
transcript's estimated read share by its effective length. Under the mapping
{transcripts $\to$ founders, reads $\to$ k-mers, effective length
$\ell\to K^w$}, kMate's per-founder step is exactly the $/\ell$ term; kallisto and
salmon carry the same effective-length correction. The old kMate EM simply
**omitted** this step.

**Key implementation points.**

- **(a) Priors stay outside the division.** The Dirichlet $(\alpha-1)$ and the
  anchor `prior_weight · total_c · prior_h` pseudocounts are added *after* the
  $/K^w$ division, so the anchor pull is governed by `prior_weight` independent of
  $K^w$.
- **(b) $K^w$ is computed over the FULL estimation unit (all k-mers, incl.
  $c_k=0$), with a $10^{-12}$ floor.** This is the unbiased normalizer: because
  $E[c_k/\mu_k]=1$ for every carried k-mer regardless of whether it draws a zero
  count, $E[\mathrm{raw}_f]=h_f\,K^w_{f,\text{full}}$, so $\mathrm{raw}_f/K^w_{f,\text{full}}$
  is unbiased for $h_f$. Summing $K^w$ over the **observed** ($c_k>0$) columns
  instead is a **survivorship bias** — it shrinks specifically for k-mer-poor /
  discriminative founders that draw a bad-luck run of zero counts, reintroducing
  the very collapse the fix targets (on a noisy $\sim$0.3× seed-mix this alone
  absorbed $\sim$24 founders vs 0 with the full-panel $K^w$). Callers that
  pre-slice `kmer_pa` to the observed columns therefore MUST pass `kfw`, computed
  over the full unit: the global driver computes `kfw_full` before the nonzero
  slice; `block_em` passes a full-window `kfw`; the bootstrap
  (`h_uncertainty.bootstrap_cov_h`) computes `kfw_full` once outside the replicate
  loop. When the full matrix (with zero-count columns) is passed directly, the
  in-solver `kmer_pa @ w` fallback is already correct.
- **(c) Exact reduction to legacy.** When $K^w$ is constant across founders the
  per-founder division is a common factor that cancels in the simplex
  renormalization, so `per_founder` reduces **byte-for-byte** to the old `global`
  behaviour.

---

## 4. Evidence

All results are on the production `filt2inv` Chr1 panel unless noted.

### 4.1 Simulation (real p231 g0 sim reads)

- Per-founder error $\leftrightarrow K^w$ correlation: **+0.84 (global) → +0.07
  (per_founder)**.
- Founders absorbed to ~0: **32 → 0**.
- h-RMSE: **~50× lower** under `per_founder`.
- On a **planted-selection** (skewed) truth, the estimator slope (est ~ true) is
  **~0.6 under global** — it *compresses* real frequency differences — and
  **~1.0 under per_founder** (preserves them).
- **Regularization alternatives were rejected.** Dirichlet / anchor-to-$p_0$
  drive the slope to **~0.006**: a prior cannot distinguish a bug-collapsed
  founder from a truly-rare one, so it flattens genuine signal. The bias must be
  removed at the source, not papered over.

### 4.2 Per-record AF factorial ablation

16 configs = {4 filters} × {multinomial, per_founder} × {$\omega$: none, $1/m_b$},
on real g0 reads:

- **Every `per_founder` config beats every `multinomial` config** on AF-MAE.
- The **old production config** (filt2inv + multinomial + $\omega=1/m_b$) is the
  **worst row**: AF-MAE **0.0096**, 33 founders absorbed.
- The **best config** is filt2inv + per_founder + no-$\omega$: AF-MAE **0.0036**
  (2.7× better), **0 absorbed**, ~100× fewer AF outliers.
- Under `per_founder`, **dropping $\omega$ helps** (AF-MAE 0.0036 no-$\omega$ vs
  0.0045 with $\omega$; h-RMSE 0.0017 vs 0.0022).
- Dropping private-singleton k-mers (filt2inv) still helps as **denoising**
  (AF-MAE 0.0036 vs 0.0044 keeping privates; eff-founders 199 vs 164).
- filt2 ≈ filt2inv (the invariant ac=F drop is neutral).

### 4.3 Real seed-mix (8 replicates)

Genome-wide = mean over 5 chromosomes.

- Collapsed founders (aggregate $h<10^{-3}$) drop consistently, e.g. per sample
  **~21–27 → ~3–11**.
- The smallest genome-wide founder frequency under the fix is **~$2\times10^{-4}$**
  with **zero** founders below $10^{-6}$ — versus the old normalization driving
  4–5 founders to $10^{-20}\ldots10^{-30}$.
- Genome-wide, sample-averaged **kMate-vs-hapFIRE agreement improves from Pearson
  $r=0.38$ (old) to $r=0.55$ (fix)**; collapsed founders **21 → 3**.

### 4.4 Why hapFIRE does not show this bias

hapFIRE fits *haplotype* frequencies per LD block via least squares on a 0/1/2
design matrix — there is **no per-founder k-mer budget** — and it averages over
hundreds of independent blocks. Its ecotype solve (verified in
`external/HapFIRE/haplotype_generation.py`) is
$\min_h \lVert H^\top h - y\rVert$ s.t. $h\ge 0$, $\sum h = 1$, which **permits
exact zeros** and indeed outputs exact zeros for selected-out ecotypes. Its
robustness to spurious collapse is **structural (many-block averaging)**, not a
floor or prior.

**Takeaway:** do NOT add an artificial frequency floor to kMate. The principled
fix is `per_founder` (remove the bias) plus averaging over units (per-chromosome
/ per-block).

---

## 5. Decision — what changed in production

- **`normalize="per_founder"` is the default everywhere** (global and window
  modes), wired through `em_solver.solve_em`, `block_em.solve_em_per_block`, and
  the `per_sample_per_chrom.py` driver via `--normalize {per_founder,global}`.
  The legacy multinomial normalization is kept as `--normalize global` for
  reproducing old runs but is **deprecated**.

- **GLOBAL mode is the GrENE-Net production mode** and is what this fix updates.
  Both modes emit per-record SNP/SV AF (the driver always projects
  $\hat{\mathbf h}\to$ AF through `var_pa`); global mode uses **one founder mixture
  per chromosome**. Because the GrENE-Net cohort is **heavily selfing** there is
  little within-sample recombination mosaic to resolve, so global mode — not
  window mode — is the production estimator for *both* the seed-mix $p_0$ **and**
  the evolved samples. Global mode therefore produces **both** the per-sample
  founder $\mathbf h$ (the $p_0$-anchored, founder-level analysis) **and** the
  per-record SNP/SV `alt_freq` tables that ship downstream — in one pass.
  Production config: **per_founder + keep filt2inv + DROP $\omega$**
  (`--kmer-weight uniform`). This fixes the founder-$h$ collapse and updates the
  AF tables; it **supersedes** the old "$\omega=1/m_b$ production" recommendation.
  The full cohort (8 seed-mix + ~2,168 evolved) is being re-run under this config.

- **WINDOW mode** (per-window EM + global anchor + HMM smoothing) is for
  *recombinant* pools and is **not** the GrENE-Net production path (the selfing
  cohort uses global). For completeness: `per_founder` is the default there too
  and is the clear winner for **pure-local** per-window estimation (`local-only`
  mode: ~18–21% lower RMSE than multinomial). Under the full `star2` recipe the
  multinomial-vs-per_founder difference is marginal/mixed — the smoothing, a
  variance-reduction step, helps the noisier multinomial estimator slightly more
  (isolated and shown to be a real **bias-variance** effect, not a bug). Since it
  is marginal and `per_founder` wins the meaningful local-only use case,
  `per_founder` is the single default; the heavy anchor/smoothing `star2` tuning
  does not warrant a second normalization mode. Local-only is the meaningful
  window mode for detecting what is *actually locally present*.

---

## 6. Carry into the paper

- **The M-step must divide by each founder's effective k-mer content
  ($/K^w_f$).** State it as the effective-length correction shared with
  RSEM/kallisto/salmon; the omission of this term is what collapsed k-mer-poor
  founders to zero.
- **The collapse was a completeness bias, not read noise** — reproduced on
  noiseless counts, correlating +0.84 with per-founder k-mer budget under the old
  normalization; +0.07 after the fix.
- **AF accuracy vs founder decomposition are separable, but both improve.** Under
  the old normalization the founder decomposition $\mathbf{h}$ was catastrophically
  wrong (collapse) while per-record AF was only mildly off ($r\approx 0.98$) —
  because AF sums $\mathbf h$ over carriers and is robust to which of several
  near-identical founders absorbs the mass. The fix corrects **both**:
  $\mathbf{h}$ dramatically (and thus $p_0$-anchored selection), and per-record
  AF-MAE by ~2.7$\times$. Report which quantity a diagnostic actually tests.
- **A prior/floor is the wrong fix.** Regularizers collapse the est-vs-true slope
  to ~0.006 because they cannot distinguish a bug-collapsed founder from a truly
  rare one. Remove the bias, then average over units.
- **hapFIRE's zero-robustness is structural (many-block averaging), not a
  safeguard** — it has no floor and does emit exact zeros. Comparisons should not
  attribute kMate's fix to matching a hapFIRE floor.
- **Global mode is the GrENE-Net production estimator (selfing cohort)** and
  produces both the founder $\mathbf h$ and the per-record SNP/SV AF in one pass;
  window mode is not used in production here. The production weighting changed to
  per_founder + filt2inv + `--kmer-weight uniform` (drop $\omega=1/m_b$),
  superseding the earlier $\omega=1/m_b$-for-global note.

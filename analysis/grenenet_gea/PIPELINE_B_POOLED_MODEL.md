# Pipeline B v2 — pooled-trajectory selection coefficients + inverse-variance climate regression

**Status: FINAL modeling spec** (frozen 2026-06-05). Implementation to follow this exactly.

**Goal.** Find SVs whose allele frequency changes with climate, using the GrENE-net
replicate common-garden time series. Revised "pipeline B" so that: (i) replicate plots
are pooled *before* the nonlinear slope fit (noise cancels instead of producing spurious
slopes); (ii) each trajectory point carries an explicit uncertainty from **sampling
(flower counts + depth)** and **drift among plots**, so the slope and its SE come from one
weighted fit; (iii) each site is weighted in the climate regression by how reliable its
estimate is.

**Currency = allele frequency** (relative selection — the climate-adaptation question).
A frequency-based slope stays a valid relative selection coefficient even as a plot
shrinks; census size enters only as *noise*, which the sampling term down-weights.
(Absolute allele-copy growth, p × flowers, is a separate demography-aware analysis, parked.)

Running example: one SV, founding frequency **p₀ = 0.12**; one warm site (bio1 = 17.5 °C)
with **4 plots** reading 0.30 / 0.05 / 0.34 / 0.28 at gen 3.

---

## Notation (indices)

- **g** — a *site* (common garden), g = 1 … 20.
- **j** — a *plot* (independent replicate population) within a site, j = 1 … n_g.
- **t** — *generation*, t ∈ {0, 1, 2, 3}. Gen 0 = the shared founding seed mix.

Everything runs **separately for each SV**; equations below are for one SV.

---

## Step 1 — pool the plots into one site frequency

$$\bar p_{t,g} \;=\; \frac{\sum_j f_{jg}\, p_{tjg}}{\sum_j f_{jg}}$$

**where:**
- **p̄_{t,g}** — *output:* the SV's pooled allele frequency at site *g*, generation *t*, in [0,1]. *(e.g. 0.18)*
- **p_{tjg}** — allele frequency kMate estimated in plot *j*, site *g*, gen *t*. *(gen-3 plots = 0.30, 0.05, 0.34, 0.28)*
- **f_{jg}** — weight of plot *j* = flowers collected from it. *(e.g. 250)* [Default: flower-weighted; equal-weight is the alternative.]
- **j** runs over the **persistent plots** of site *g* — the plots present at **all of gen 1, 2, 3** (so the trajectory tracks a fixed set; this is the chosen handling of plot extinction, see Limitations). Gen 0 is fixed: p̄_{0,g} = p₀ = 0.12 for every site.

Result: one trajectory per site, e.g. `[0.12, 0.18, 0.22, 0.24]`.

*Estimand note:* because pooling precedes the logit, `s_g` is the log-odds slope of the
**expected site frequency**, not the mean of per-plot selection coefficients (logit is
nonlinear, so the two differ by Jensen). The expected-frequency slope is the intended target.

---

## Step 2 — uncertainty of each trajectory point (variance-components, on the log-odds scale)

Each timepoint of a site's trajectory has a variance with two sources: **binomial sampling**
(flowers + depth) and **drift among plots**. These are **not additive** — the observed
among-plot spread already *contains* the sampling noise — so we combine them as a
variance-components estimate (per-plot variance = the larger of the two), and we work
directly on the **log-odds scale** to avoid delta-method linearization of a large spread.

Per plot *j*, binomial sampling variance of `logit(p_{tjg})`:

$$v^{\text{samp}}_{tjg} \;=\; \frac{1}{N^{\text{eff}}_{tjg}\,\bar p_{t,g}(1-\bar p_{t,g})}, \qquad \frac{1}{N^{\text{eff}}_{tjg}} = \frac{1}{2\,\text{flowers}_{tjg}} + \frac{1}{\text{coverage}_{tjg}}$$

Among-plot (drift) variance, computed on the log-odds scale:

$$\sigma^2_{t,g} \;=\; \operatorname{Var}_j\!\big[\operatorname{logit}(p_{tjg})\big] \quad(\text{shrunk: } \tilde\sigma^2_{t,g}=\tfrac{(n_g-1)\sigma^2_{t,g}+d_0\,\sigma^2_{\text{pool}}}{(n_g-1)+d_0})$$

Variance of the **site-mean** log-odds (the value used in Step 3):

$$V_{t,g} \;=\; \frac{\max\!\big(\tilde\sigma^2_{t,g},\; \bar v^{\text{samp}}_{t,g}\big)}{n_g}$$

**where:**
- **N^eff_{tjg}** — effective chromosomes sampled in plot *j*: bounded by individuals (`2·flowers`) and reads (`coverage`). Few flowers (a crashing plot) → small N^eff → large sampling variance. *(flowers→genomes factor is tunable; flowers over-counts independent chromosomes for half-sibs; `coverage` is the per-plot genome-wide mean, not per-SV local depth — both documented approximations.)*
- **v^samp_{tjg}** — per-plot binomial sampling variance on the log-odds scale (binomial delta `1/[N p(1−p)]`). This is the *irreducible sampling floor*; it omits kMate's deconvolution uncertainty, so we **sanity-check** it against kMate's own per-record `se` on a sample of SVs.
- **σ²_{t,g}** — among-plot disagreement on the log-odds scale (drift + residual sampling). **v̄^samp_{t,g}** — mean of the plots' `v^samp`.
- **max(σ̃², v̄^samp)** — variance-components combination: uses the observed spread when drift is real, falls back to the binomial floor when plots happen to agree (no double-counting).
- **σ²_pool, d₀** — shrinkage target = per-SV pooled among-plot variance; prior df **d₀ = 4** (sites with n_g = 1 fully borrow the pool).
- **n_g** — number of (persistent) plots in site *g*.
- **V_{t,g}** — *output:* log-odds-scale variance of the site-mean point; weights the Step-3 fit.

**Gen 0:** the shared seed mix — no among-plot term. Its variance is the among-replicate
variance of `logit(p₀)` over the 8 SEEDMIX reps, divided by 8. (Caveat: these are
*technical* replicates, so this under-states each garden's true founding-draw uncertainty;
gen 0 therefore gets a high but not infinite weight, anchoring the slope near p₀.)

---

## Step 3 — within-site model: weighted slope of log-odds vs generation

For each site *g*, fit the line over t ∈ {0,1,2,3} by **weighted least squares**, weights ω_{t,g} = 1/V_{t,g}:

$$\operatorname{logit}(\bar p_{t,g}) \;=\; a_g \;+\; s_g\, t \;+\; \varepsilon_{t,g}, \qquad \operatorname{Var}(\varepsilon_{t,g}) = V_{t,g}$$

**where:**
- **logit(p̄_{t,g}) = ln[ p̄_{t,g} / (1 − p̄_{t,g}) ]** — *the outcome (y):* log-odds of the site frequency. *(logit 0.18 = −1.52)*
- **t** — *the predictor (x):* generation 0, 1, 2, 3.
- **a_g** — intercept = log-odds at gen 0 (nuisance). *(logit 0.12 = −1.99)*
- **s_g** — *the slope = the per-site selection coefficient (quantity of interest).* *(e.g. +0.30; positive = rising)*
- **ω_{t,g} = 1/V_{t,g}** — the weight of each generation in the fit: well-measured generations (small V) count more.

Both `s_g` and its standard error come from **this same weighted fit**:

$$s_g \;=\; \frac{\sum_t \omega_{t,g}\,(t-\bar t_\omega)\,\operatorname{logit}(\bar p_{t,g})}{\sum_t \omega_{t,g}\,(t-\bar t_\omega)^2}, \qquad
\operatorname{SE}(s_g)^2 \;=\; \frac{1}{\sum_t \omega_{t,g}\,(t-\bar t_\omega)^2}$$

- **t̄_ω = Σ_t ω_{t,g} t / Σ_t ω_{t,g}** — the weighted mean generation.
- **SE(s_g)** — *output:* how shaky `s_g` is, straight out of the weighted regression. *(e.g. 0.12)*

(If all V_{t,g} were equal, t̄_ω → 1.5 and these reduce to the simple `Σ(t−1.5)·logit / 5`
slope — so the unweighted version is just the equal-variance special case.)

*Caveat:* with 4 points, known variances, and 2 parameters there is no residual df, so
`SE(s_g)` reflects only the input point variances and **assumes the logit-linear-in-t model
is correct** (it does not capture trajectory curvature / non-constant s).

---

## Step 4 — between-site model: site `s` vs climate, weighted by reliability

One **random-effects meta-regression** across the **20 sites**:

$$s_g \;=\; \beta_0 \;+\; \beta_1\, \text{climate}_g \;+\; \eta_g, \qquad \text{weight } w_g = \frac{1}{\operatorname{SE}(s_g)^2 + \hat\tau^2}$$

**where:**
- **s_g** — *the outcome (y):* each site's slope from Step 3 (20 values). *(this site: +0.30)*
- **climate_g** — *the predictor (x):* site climate, bio1 (mean annual temperature), standardized across sites. *(17.5 °C)*
- **β₀** — intercept (nuisance).
- **β₁** — *the slope = the headline answer:* change in `s` per unit of climate. β₁ > 0 ⇒ allele rises faster in warmer sites = climate-dependent selection.
- **η_g** — residual (a site's departure from the climate line).
- **τ̂²** — **between-site heterogeneity** beyond climate (other environment, local selection), estimated by DerSimonian–Laird. Including it (random- not fixed-effect) stops a few high-precision sites from dominating β₁.
- **w_g = 1/(SE(s_g)² + τ̂²)** — *the weight:* site reliability, now folding in true between-site spread.
- **Robustness:** also report β₁ under `1/SE²`, `n_g`, and equal weights — if the estimate and the top blocks are stable across these, the weighting choice is not driving the result.

**Significance of β₁:** site-level permutation — reassign the `climate_g` labels across the 20
sites, refit, record β₁; repeat thousands of times → null distribution. Report a two-sided p
(any climate association) and a one-sided p (up-in-warm).

---

## Chain summary

| Step | what it produces | y | x | weights |
|---|---|---|---|---|
| 1 | site frequency `p̄_{t,g}` | — | — | flowers `f_{jg}` |
| 2 | point variance `V_{t,g}` (sampling + drift) | — | — | — |
| 3 | per-site slope `s_g` **and** `SE(s_g)` | log-odds of `p̄_{t,g}` | generation `t` | `1/V_{t,g}` |
| 4 | climate effect `β₁` + permutation p | `s_g` | site climate (bio1) | `1/SE(s_g)²` |

## Locked choices

1. **Persistent plots** (present at gen 1,2,3); **flower-weighted** pooling (Step 1).
2. **ε = 1e-3** logit clip.
3. **Frequency** currency; variance-components combination of binomial sampling + drift, on the log-odds scale (Step 2); absolute-copy-number / demography kept separate.
4. **Shrinkage** of drift variance toward the per-SV pool, prior df **d₀ = 4**.
5. **Weighted (GLS) within-site fit** (Step 3): slope and SE from one model, weights `1/V_{t,g}`.
6. **Random-effects** between-site weight `1/(SE²+τ̂²)`, τ̂² by DerSimonian–Laird (Step 4); β₁ also reported under 1/SE², n_g, equal weights.
7. **Site-level permutation** null for β₁; two-sided + one-sided up-in-warm.

## Diagnostics to run alongside

- **Sampling-floor check:** `√(v^samp)` vs kMate's per-record `se` on a sample of SVs (is the binomial floor realistic?).
- **Mean–variance check:** does `SE(s_g)` (or σ²) track `|s_g|` / frequency across sites? If high-effect sites are systematically noisier, inverse-variance weighting could fight the signal.
- **Weighting robustness:** β₁ and top blocks under the four weight schemes (Step 4).

## Deferred improvement — permutation-calibrated WZA (replaces the spline correction)

The downstream block step uses Booker's WZA with a **degree-2 polynomial** SNP-number
correction (`wza_script.py::adjust_WZA_with_spline`). That correction is fragile: in sparse
large-window tails the quadratic can predict a **negative SD → NaN block p** (confirmed by
the WZA author, T. Booker, by email; his fix = cap windows + degree-7 polynomial). We have
**not** hit it (0 NaN across our 7,249 blocks, max block 918 SVs), and we added a
**positive-SD floor** safety net in `adjust_WZA_with_spline` so it can never NaN.

**Better fix, deferred:** because pipeline B already has a **site-permutation null**, we can
calibrate the block WZA *by permutation* instead of the polynomial:
- compute the raw WZA score per block on the real data **and** on each climate-label
  permutation (re-aggregating the permuted per-SV statistics to blocks);
- block p = empirical tail of its own permutation distribution.

A large window has a large WZA under both real and permuted climate, so the permutation null
**self-calibrates for window size** (what the polynomial approximates) — no degree to choose,
no negative SD, never NaN — and it preserves within-block LD. It also unifies calibration
with the per-SV β₁ test (same permutation). NOT Fisher's method: Fisher assumes within-window
independence (false inside an LD block) and has no size correction, and is less powerful than
weighted-Z (Whitlock 2005; Booker's own testing). To implement: reuse the `beta_perm` array
from Step 4, map per-SV → block WZA per permutation, take the empirical tail.

## Known limitations (stated, not fixed here)

- **Plot extinction / survivorship (MNAR).** Restricting to persistent plots conditions on
  survival; if survival is allele- or climate-dependent, `s_g` is a survival-conditioned
  estimate and the bias can be climate-dependent. A dedicated missingness/survival analysis
  is deferred. (So "census enters only as noise" holds only absent survival-dependent loss.)
- **Drift autocorrelation.** The log-odds trajectory of a fixed plot set is a random walk, so
  its level-vs-time residuals are autocorrelated; the diagonal `V_{t,g}` ignores this. It is
  not modeled (4 points can't support a temporal covariance) and only affects the *weights*,
  not the permutation-based β₁ test.
- **N ≈ 20 sites.** β₁ rests on ~20 between-site contrasts; per-SV genome-wide power is
  inherently limited (handled downstream by block-level WZA aggregation).
- **Single climate axis (bio1)** for now; other bioclim variables and their multiple-testing
  burden are a separate sweep.

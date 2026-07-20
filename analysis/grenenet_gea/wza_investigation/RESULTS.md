# WZA investigation — results (gen1 SNP × bio1, Kendall-τ)

Inputs: `analysis/grenenet_gea/phase1_replication/results/kendall/kendall_snp_gen1_bio1.csv`
(1,989,384 SNP records, MAF≥0.05, 16,467 hapFIRE blocks). Code: `wza_core.py`,
`run_investigation.py`, `centromere_check` (inline). Outputs + figs:
`analysis/grenenet_gea/wza_investigation/results/`.

SNPs/block fed to WZA: median **14**, q75 53, q95 656, **max 9,158** (atomized kMate
records inflate vs the 3,028-SNP raw hapFIRE max).

---

## Is the heavy tail centromeric? — SHOWN, and the answer is "partly" (fig4)

Mapped every block to its genomic midpoint and distance to the nearest TAIR10
centromere (Chr1 15.0, Chr2 3.6, Chr3 13.8, Chr4 3.95, Chr5 11.75 Mb).

| block size | within 1 Mb | 2 Mb | 3 Mb | 5 Mb of a centromere |
|---|---|---|---|---|
| > 3000 SNPs (n=55) | 36% | 64% | **80%** | 91% |
| > 2000 SNPs (n=172)| 21% | 41% | 57% | 73% |
| > 1000 SNPs (n=537)| 9% | 23% | 35% | 53% |

**The extreme tail (>3000 SNPs) is mostly centromeric (80–91% within 3–5 Mb)** —
fig4 shows the vertical spike of huge blocks at each red centromere line. **But
the broader large-block set is not**: 43% of >2000-SNP blocks (74/172) sit >3 Mb
from any centromere — they are **arm blocks** (high-LD / low-recombination / sweep
regions). The biggest arm blocks: 2_729 (Chr2:6.7 Mb, 4715 SNPs), 5_724
(Chr5:8.5 Mb, 4490).

**CAM5's own large neighbor `2_1264` (2,866 SNPs) is at Chr2:11.4 Mb — 7.8 Mb from
the centromere. It is an ARM block, NOT pericentromeric** (the STATUS.md
"pericentromeric" label was wrong).

➜ **Consequence: a position-based centromere filter would NOT cleanly fix WZA.** It
misses ~half the problematic large blocks and risks dropping real arm signal. The
**SNP-count cap** targets the actual failure cause (count, not location) and is the
cleaner fix.

## Degree × cap sweep (sweep_summary.csv)

| variant | NaN p | neg-SD windows | Bonferroni hits |
|---|---|---|---|
| **deg2_nocap (canonical)** | **18** | **18** | 63 |
| deg7_nocap (phase-1) | 0 | 0 | 10 |
| **deg2_cap2000** | **0** | 0 | 25 |
| deg7_cap2000 | 0 | 0 | 6 |
| deg2_cap(q95=656) | 0 | 0 | 14 |
| deg2_cap(q75=53) | 0 | 0 | 3 |

- **Canonical deg-2 with no cap = the email bug**: 18 windows get NaN (negative
  predicted SD in the heavy tail). Capping at ANY level removes all NaN.
- deg-7 also removes NaN (the phase-1 workaround) but by over-fitting the tail.
- Bonferroni count is unstable across choices (3–63) — the correction choice
  materially changes the candidate list.

## Worry #3 — does WZA over-weight LOW-SNP windows (CAM5, 13 SNPs)? **NO — opposite**

Fraction of windows with corrected Z_pVal < 0.05, by window SNP count
(overweighting_by_snpbin.csv, fig3). Canonical **deg2_nocap**:

| SNPs/block | 2–5 | 6–10 | 11–20 | 21–50 | 51–100 | 101–200 | 201–500 | 501–1k | >1k |
|---|---|---|---|---|---|---|---|---|---|
| frac p<0.05 | **0.0005** | 0.009 | 0.044 | 0.139 | 0.230 | **0.281** | 0.236 | 0.072 | 0.056 |

Small windows are **massively UNDER-powered** (2–5 SNPs: 0.05% significant, 100×
below uniform), power **peaks at 100–200 SNPs**, then collapses in the tail (where
NaN lives). So canonical WZA does **not** favor small windows like CAM5 — it
penalizes them. CAM5 reaches significance because its 13 SNPs are genuinely,
coherently associated (raw Z = **17.75**, i.e. ~18σ pre-correction).

**The deg-7 correction is what boosts CAM5-sized windows.** deg7_nocap flattens the
curve and raises the 11–50-SNP bins relative to deg-2 (11–20 bin: 0.087 vs 0.044).
CAM5 (13 SNPs) lives exactly there.

### CAM5 (block 2_1265, 13 SNPs, raw Z=17.75) across variants (cam5_across_variants.csv)

| variant | predicted SD@13 | Z_pVal | rank /16179 |
|---|---|---|---|
| deg2_nocap | 3.93 | 1.5e-3 | 335 |
| **deg7_nocap (phase-1)** | 3.29 | **5.7e-5** | **28** |
| deg2_cap2000 | 3.67 | 5.1e-4 | 174 |
| deg2_cap(q95) | 3.31 | 7.2e-5 | 36 |

**Adjudication:** the *empirical* rolling SD of Z near 13 SNPs is **3.45**. deg-2
predicts 3.93 (slightly conservative), deg-7 predicts 3.29 (slightly
anti-conservative). So CAM5's honest corrected p is **~1e-4 to 1e-3** — a real
top-0.2–1% block, but **NOT the p≈3e-8 the phase-1 deg-7 run reported**; that
extra ~4 orders of magnitude is deg-7 tail over-fitting. CAM5 is a true signal,
just over-stated by the deg-7 correction.

> Caveat: WZA on raw Kendall does NOT correct population structure (LFMM does); part
> of CAM5's raw Z=17.75 could be climate-correlated structure. Separate issue.

## deg-2 vs deg-7, fit quality (fig6/fig7, RMSE to rolling support)
| | NO CAP RMSE(SD) | CAP2000 RMSE(SD) | min predicted SD over support |
|---|---|---|---|
| deg-2 | 2.06 | 1.60 | +2.58 (nocap) — but **<0 when EXTRAPOLATED past the ~4,700-SNP support** → NaN |
| deg-7 | **0.94** | **0.89** | follows the curve; **explodes upward past 4,500** (nocap) |

**deg-7 genuinely fits the empirical curve better** — the author was right a quadratic
is too rigid. The NaN is an **extrapolation** artifact (deg-2 turning negative past
the rolling support out to the 9,158-SNP blocks), NOT a within-range fit failure. The
**cap removes the extrapolation → no NaN for any degree** (fig7 right panel).

## Full cap-vs-no-cap, ALL blocks (compare_cap_nocap.py, fig8) — REVISES the "must cap" claim
| regime | NaN | GIF | Bonf | BH q<.05 |
|---|---|---|---|---|
| deg2 no-cap (+floor) | 0 | **1.166** | 81 | 295 |
| deg2 cap2000 | 0 | 1.061 | 25 | 176 |
| deg7 no-cap (phase-1) | 0 | **0.974** | 10 | 35 |
| deg7 cap2000 | 0 | 0.956 | 6 | 36 |

Agreement (Jaccard of BH-sig sets): deg7 cap↔nocap **0.78**; deg2 cap↔nocap 0.56;
deg2cap↔deg7cap **0.15**. Spearman(−log p) ≥0.93 everywhere (ranking preserved;
threshold-crossing is what moves).

**Concession: for the phase-1 pipeline (deg-7) the cap barely matters** — fig8 right
panel is a tight diagonal, GIF ~0.96 either way, CAM5 BH-sig in all regimes. Results
hold cap or no-cap (matches the project's own check). **The cap only earns its keep
under deg-2**, where fig8 left panel shows the 18 giant blocks (4,700–9,158 SNPs)
stranded at −log10 p ≈ 300: the **SD safety-floor hack crushes them to p≈1e-300**
(spurious genome-wide hits at centromeric/large-LD blocks) → GIF 1.17, 295 BH. Cap (or
dropping the floor) removes them. **The real lever is the DEGREE, not the cap**
(deg2 vs deg7 = 176 vs 36 BH, Jaccard 0.15); deg-7 is the better-calibrated, cap-robust
choice phase-1 actually used.

## Bottom line
1. WZA was a reasonable choice; the failure is the **window definition** (hapFIRE
   BigLD blocks: 1→9,158 SNPs) + **not using WZA's SNP cap**.
2. The heavy tail is **partly centromeric (extreme tail yes, broader large blocks
   no)** → centromere filtering is not the clean fix; **cap SNPs/window** is.
3. **The cap is the load-bearing fix; the polynomial degree is secondary.** With a
   cap, deg-7 (author's choice) is defensible (fits better) and deg-2 (canonical) is
   more stable — Bonferroni swings 25↔6 between them. Best adjudicated by the **local
   empirical null** (poly-free), which both only approximate.
4. WZA does **not** overweight low-SNP windows; small windows are *under*-powered.
   CAM5 polynomial-free truth p ≈ **1.6e-4** (real top-~0.1%): deg-2 too conservative
   (1.5e-3), deg-7 too liberal (5.7e-5) but closer; phase-1's p≈3e-8 is uncapped deg-7
   tail over-fit. Recommend **SNP cap (~2000 or q95) + validate degree vs empirical
   null**; drop the SD safety-floor hack.

## Pseudo-replication / panel-density (pseudoreplication_test.py, fig10) — answers "kMate≠phase1 even on Kendall"
WZA denom sqrt(Σpq²) assumes SNP INDEPENDENCE; correct correlated-Z denom = sqrt(pq'R pq).
Measured from the 326-pool AF matrix:
| block | N | mean\|r\| | M_eff | WZA-Z | corr-Z | inflation |
|---|---|---|---|---|---|---|
| 4_2519 (Chr4 CRK) | 1169 | 0.52 | **58** | 160.2 | 9.1 | **17.5×** |
| 2_1265 (CAM5) | 13 | 0.76 | 3 | 17.8 | 5.8 | 3.1× |

- Within-block raw-Z scales as **√N** (fig10 left tracks the perfect-LD line), so a denser
  panel inflates a high-LD block ~√(SNP-ratio): kMate 1169 vs phase-1 ~377 → √3.1≈1.76×;
  remaining gap to the observed 2.9× is kMate-vs-hapFIRE AF differences (~1.6×). BOTH, not just density.
- NOT "fabrication": the SNP-number correction standardizes Z vs the SNP-count-matched null
  (also √N-inflated), so it removes the AVERAGE pseudo-replication. It fails only for
  atypical-LD blocks (CRK mean|r|=0.52 >> genome avg) and across panels (correction trained
  on one panel's N-distribution doesn't transfer).
- Fix options (we HAVE the AF matrix): (1) correlation-aware weighted-Z sqrt(pq'Rpq) [principled,
  threshold-free, panel-robust]; (2) M_eff-recalibrated correction; (3) LD pruning [crude];
  (4) lead-SNP/block [conservative, = lib.collapse_to_blocks]. Even R-aware BOTH blocks stay sig
  (9.1, 5.8) — signals real but overstated; 4_2519 still > CAM5 → structure (LFMM) still open.

## BACK TO BASICS: is the raw AF→climate signal real? (raw_signal.py, fig11) — the deepest finding
The Kendall is SPATIAL: per-SNP τ between gen1 AF and pool bio1 across 326 pools.
"Increasing/decreasing in frequency" = with temperature across sites.

| block | N | pool τ (326 pools) | pool p | **site τ (31 sites)** | **site p** | coherence |
|---|---|---|---|---|---|---|
| 4_2519 (Chr4 CRK) | 1169 | -0.215 | 1.2e-8 | **-0.316** | **0.012** | 74% same sign |
| 2_1265 (CAM5) | 13 | +0.217 | 8.1e-9 | **+0.308** | **0.015** | 100% same sign |

**The signal IS real but MODEST.** CAM5 ALT freq rises with temperature (τ≈+0.3), Chr4
CRK falls (τ≈-0.3); both biologically coherent. **But the genome-wide-significant
p-values are pseudoreplication at TWO nested levels:**
1. **DOMINANT & UPSTREAM of WZA: pools-within-sites.** 326 pools but only **31
   independent climates**. Per-SNP Kendall across 326 pools inflates p ~6-7 orders
   (8e-9) vs the honest 31-site test (**p≈0.015**). WZA cannot fix this — it's upstream.
2. SNPs-within-LD-blocks (the WZA independence issue; 17.5×/3× from pseudoreplication_test).

**Resolves "was WZA the right choice":** the conceptual paradox is real — blocks are
DEFINED by LD (correlation) yet WZA's denominator assumes within-block independence.
CAM5 is 100% coherent (M_eff=3 of 13) = ONE haplotype signal, not 13 tests. WZA tackles
the SECONDARY (SNP) pseudoreplication while the PRIMARY (site) one dominates and is
untouched. For LD-defined blocks the honest unit is haplotype-level / lead-SNP / site-level
— which is exactly the project's two-stage site-permutation plan (memory: pseudoreplication
wall N≈20 sites not pools). NB phase-1's PUBLISHED headline used LFMM(K=16) which absorbs
site structure as latent factors, so CAM5's published p is on firmer ground than this raw
Kendall→WZA path (the most pseudoreplicated one). Honest effect size: a modest climate
cline, τ≈0.3.

## Genome-wide PC1 + site-permutation (pc1_genomewide.py, pc1_manhattan.ipynb, fig12) — THE HONEST RESULT
One test/block (PC1 = dominant haplotype axis, correlation-PCA) + site-permutation null
(climate shuffled across 31 sites). Removes BOTH pseudoreplications (SNP-count + pool-in-site).
gen1 × bio1, 16,467 blocks, 20k perms:
- **0 blocks pass BH q<0.05** (min bh_q ≈ 0.15; perm floor −log10(5e-5)≈4.3).
- CAM5 2_1265: PC1-VE 0.81, r=+0.40, perm_p 0.028, **bh_q 0.40**.
- Chr4 CRK 4_2519: PC1-VE 0.59 (not a clean single haplotype), r=−0.48, perm_p 0.0075, **bh_q 0.29**.
- The raw-WZA #1 peak (p≈1e-11) is NOT significant once corrected.
fig12: raw WZA deg-7 Manhattan (34 BH-sig peaks, Chr4 CRK at ~1e-11) vs PC1+site-perm
(flat, 0 BH-sig). NOT "no climate adaptation" — it's the weakest slice (1 axis, 1 gen, no
temporal Δp, no multi-axis/model aggregation, no LFMM). It shows the individual
genome-wide-significant p-values in the raw Kendall→WZA path were largely pseudoreplication;
honest per-block single-axis evidence is modest (nominal, not BH-sig). Real power needs the
multi-axis/temporal/site-permutation two-stage design (project plan).

## QQ / calibration (pc1_manhattan.ipynb, fig13) — λ alone is misleading; read the SHAPE
λ_GIF: WZA deg-2+floor **0.58**, WZA deg-7 **0.81**, PC1+site-perm **1.71**. The naive
"WZA inflated / PC1 calibrated" is WRONG by λ. Honest reading:
- WZA deg-2+floor: bulk DEFLATED (λ0.58) but a few giant floored blocks explode to
  −log10 p≈300 = catastrophic tail false-positives, not whole-genome inflation.
- WZA deg-7: bulk ~OK (λ0.81) but TAIL rides far above diagonal (−log10 p≈11) = the
  pseudo-replication peaks BH calls (34 hits).
- PC1+site-perm: NO runaway tail (capped at perm floor 4.3) — WZA peaks gone — but whole
  line sits mildly above diagonal (λ1.71) = a DIFFUSE excess. This is NOT pseudo-replication
  (site-perm removed it); it's residual genome-wide STRUCTURE + genuine POLYGENIC climate
  covariation (caveat #3), which site-perm can't separate — only LFMM can. Broad-and-weak,
  so BH still 0 hits. Diffuse λ>1 + 0 FDR = classic polygenic-or-structure signature.
Two failure modes shown cleanly: WZA = a tall artifactual TAIL; PC1+site-perm = no tail,
just an unresolved mild genome-wide shift needing LFMM.

## RECURRENCE (2026-07-03) — the same SD-floor bug resurfaced on clq0.9 blocks

The `phase1_replication/clq90` run initially assumed the finer clq0.9 blocks (max
5,168 SNPs vs the 9,158 here) were small enough that canonical deg-2 wouldn't need
this section's cap/floor fix at all. It still hit the exact bug documented above
(negative-SD extrapolation in the sparse large-block tail, silently floored to
fabricated significance) — just at a smaller absolute block size (~1,150+ SNPs snp,
~400+ nonsnp instead of ~4,700+). It drove most of that run's "cross-model
reproducible core" headline. Fixed there by applying this section's own
recommendation (drop the floor hack, cap within the polynomial's support) with
class-specific caps (1000/350) rather than reusing 2000. See
`../phase1_replication/STATUS_clq90.md` §0 BUG + FIX for the full writeup.
**Takeaway for future WZA runs on any new block definition: do not assume a smaller
max block size means the floor hack is safe — check the negative-SD count
directly** (`wza_core.apply_correction(..., sd_floor=False)`, count `neg_sd`) rather
than trusting "0 NaN" from the floored production script, which cannot distinguish
"clean fit" from "floor hack papering over a bad fit."

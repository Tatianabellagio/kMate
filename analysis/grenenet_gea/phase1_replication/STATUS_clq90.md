# Status — phase-1 GEA replication on **clq0.9 (r²≥0.9) LD blocks**

> ## ⚠ CORRECTION (2026-07-03) — the original deg-2 "0 NaN, workaround not needed"
> ## result below was WRONG. A real WZA fitting bug inflated ~10-25 large blocks
> ## per output to spurious p≈0. Fixed; corrected results replace the tables below.
> See **"§0 BUG + FIX"** right after this box for the full story; all numbers
> further down in this file are already updated to the corrected regime.

> ## ⚠ UPDATE (2026-07-21) — the deg-2 "fix" above was INCOMPLETE. Removing the
> ## SD-floor hack + capping stopped the *negative*-SD→p≈0 path, but deg-2's SD
> ## polynomial still **turns over and under-predicts** near the cap on these
> ## sparse-tail clq0.9/mcf90 blocks → it still fabricated `Z_pVal==0` for
> ## unremarkable (cap-size) blocks and NaN in the sparse SV class. Confirmed it
> ## fabricates **even under a signal-free permutation null** (climate shuffled
> ## across the 31 sites); the null SD-vs-SNP-count is genuinely
> ## monotone-up-then-plateau (isotonic-R²≈0.98). **PRODUCTION CORRECTION IS NOW
> ## `--sd-fit isotonic`** (monotone-non-decreasing SD + empirical-interp mean,
> ## `wza_script.adjust_WZA_with_spline`): 0 fabricated p==0 / 0 NaN across all
> ## 189 bio1+multiaxis outputs. deg-2/deg-7 kept as options only. Full validation:
> ## `../wza_investigation/wza_sd_fix_test.ipynb` (+ permutation null + independent
> ## audit). BH-sig counts drop sharply (deg-2 over-called the whole large-block
> ## tail): snp kendall 30→2, binomial 35→2, lfmm 77→10. Manhattan/compare default
> ## regime is now `isotonic`.

> ## ⚠ UPDATE (2026-07-27) — the "★ RESULTS" (CAM5 + cross-model recurrence) and §5
> ## (nonSNP-specific candidate genes) sections below were never re-verified against
> ## the isotonic fix above and are **STALE**. Checked directly against the current
> ## production outputs (`results/clq90/wza/wza_*_gen9_bio1_isotonic.csv`, 3-class):
> ## - **CAM5 no longer reaches BH-FDR in any model×class.** Best q is now **0.346**
> ##   (lfmm/snp) — not the q=0.066–0.083 "just below FDR" cited below (that was
> ##   deg-2, pre-isotonic). Recomputed directly from the isotonic files, not from
> ##   any cached compare table.
> ## - **The "reproducible core" (Chr1_2343/GAPC2, Chr4_6307/CRK13-16, Chr3_6144/
> ##   UMAMIT32) does not survive.** Zero blocks are BH-sig in all 3 models for ANY
> ##   class under isotonic (down from the claimed 63/66/35 significant snp blocks
> ##   per model under deg-2 to 1–21). Chr3_6144 is not significant anywhere now;
> ##   Chr1_2343/Chr4_6307 survive in only 2–3 of 9 model×class combos.
> ## - **§5's "437 blocks → 548 genes, heat-stress dominates" is a separate, also-
> ##   stale artifact**: `multiaxis/nonsnp_specific_genes.py` and `overlap_table.py`
> ##   hardcode `_deg2.csv` + the retired pooled 2-class (`snp`/`nonsnp`) split, never
> ##   updated when the pipeline moved to isotonic + 3-class (`snp`/`sv`/`smallindel`).
> ##   Recomputed on the current isotonic + 3-class outputs (`{cls}_specific_peaks_
> ##   {bonf,fdr}.csv`, already on disk from `multiaxis/_build_snp_vs_nonsnp_peaks_nb.py`):
> ##   union FDR = **327 blocks → 322 genes** (nonsnp-pooled alone: 217→212). Only
> ##   **133/547 (24%) of the old gene list survives**; the specific heat-stress genes
> ##   driving the old narrative (HSFA2, COR15A, VRN2) are **absent** from the
> ##   corrected list (HSBP1/AT4G15802 is the one survivor). No new functional-theme
> ##   read has been done on the corrected list yet — do not cite "heat-stress
> ##   dominates" until one is.
> ## Superseded scripts: `multiaxis/nonsnp_specific_genes.py`, `multiaxis/overlap_table.py`,
> ## `multiaxis/_build_nonsnp_genes_nb.py` (+ the notebook it built,
> ## `notebooks/nonsnp_specific_genes.ipynb`) — all deg-2/2-class. Current equivalents:
> ## `multiaxis/class_peaks_overlap_table.py` (overlap) and the `class_specific()`
> ## helper in `multiaxis/_build_snp_vs_nonsnp_peaks_nb.py` → `snp_vs_nonsnp_new_peaks.ipynb`
> ## (candidate genes), both isotonic + 3-class.

> # 🚨 CRITICAL (2026-07-27) — THE BLOCK PARTITION IN THIS RUN IS WRONG. EVERY
> # RESULT BELOW WAS COMPUTED ON A FRACTION OF THE DATA.
>
> Decision #1 of this run ("blocks don't tile the genome → drop gap variants") is a
> **methodological error, not a defensible trade-off.** Measured directly on the
> current per-record data, strict interval containment (`lib.assign_clq_blocks`,
> which returns `''` for anything between two BigLD islands) silently discards:
>
> | class | records | kept | **DROPPED** |
> |---|---|---|---|
> | snp | 2,016,071 | 1,214,880 | **39.7%** |
> | **sv** | 27,586 | 13,973 | **49.3%** |
> | smallindel | 676,203 | 403,719 | **40.3%** |
> | nonsnp | 703,789 | 417,692 | **40.7%** |
>
> **Half of all SVs never entered the climate GEA** — for the variant class this
> whole project exists to study.
>
> **Why it is wrong.** The per-variant models (kendall / lfmm / quasi-binomial) run
> on *every* record; WZA is only an aggregation step to control inflation. A variant
> not falling inside a BigLD island is not unusable — it just needs a window. And
> the dropped records are not in distant unlinked territory: the median dropped
> record sits **731 bp** from the nearest block edge, 57% are within 1 kb, and the
> median inter-block gap is **267 bp**. They were dropped only because they were not
> in the set used to *define* the blocks (BigLD is run on a MAF>0.05, call≥90%,
> one-variant-per-position subset — 291k of 700k records on Chr1).
>
> **HapFM — the software that builds these haploblocks — explicitly does not do
> this.** `HapFM/bin/utility_functions.py::convert_fine_genomewide_breakpoints`
> converts the common-set breakpoints back to genome-wide coordinates so the
> partition covers every variant (`left = common_index[prev_right] + 1`; first block
> starts at 0, last ends at `r-1`), then merges sub-2-variant blocks. Block
> *definition* (on good-quality common variants) and variant *assignment* (all of
> them) are deliberately separate steps in HapFM. **This pipeline did the first and
> skipped the second.**
>
> **Fix:** `analysis/grenenet_gea/blocks_tiling.py` implements HapFM's conversion in
> position space — identical r²-derived LD boundaries, **0.0% of records dropped**
> at both r²=0.9 and r²=0.5. All GEA results in this file predate that fix and must
> be regenerated on the tiling assignment before being used or cited.

> # ✅ SETTLED (2026-07-28) — WZA correction regime: **tiling · isotonic SD · `const` mean · NO CAP**
>
> ## The short version: this is upstream's setup with TWO forced changes
>
> | | upstream WZA | us | why |
> |---|---|---|---|
> | cap | **none** (`--sample_snps` default 0, and never passed in ANY invocation in their repo) | **none** | same as upstream |
> | SD fit | deg-2 | **isotonic** | deg-2's SD goes NEGATIVE on our block-size range |
> | mean fit | deg-2 | **const** | deg-2's mean extrapolates to −289 on our range → p==0 |
>
> Both deviations are forced by ONE fact: Booker validated on windows with tightly
> distributed SNP counts (10 kbp bins; 302 pine genes in the demo). **Our LD blocks span
> 2 → 20,365 variants.** The SNP-number correction is just an interpolation of a nuisance
> curve; across four orders of magnitude with a sparse tail, any polynomial extrapolates off
> a cliff. Nothing else about the method changes.
>
> Retracted along the way (recorded so it is not repeated): v1 of the decision notebook fit
> the diagnostics to raw p instead of the rank-transformed p the code corrects; `isotonic_auto`
> was briefly the pick for the mean and failed audit (below); and the cap question could have
> been closed in one grep of upstream's invocations instead of a parameter sweep.
>
> Evidence: `notebooks/cap_poly_decision.ipynb` (rebuilt by
> `multiaxis/_build_cap_poly_decision_nb.py`), 2 block defs × 4 classes, extended to
> 3 models × 2 defs × 2 classes for the mean fit.
>
> **⚠ First: v1 of that notebook fit every diagnostic to the WRONG STATISTIC.**
> `run_wza.py` passes the **raw** per-record p as `--summary_stat`; `wza_script.py:167`
> then **rank-transforms it genome-wide** (`csv["pVal"] = csv[stat].rank()/n`) before the
> z. v1's `wza_Z()` helper used raw p. That is not cosmetic — under the rank transform
> E[z]=0 by construction, so mean(Z)-vs-block-size is flat-and-noisy (≈0, dipping to −5 in
> the tail) instead of the smooth 34→66 rise seen on raw p. Two of v1's verdicts flip.
>
> **Checked against upstream** (`github.com/TBooker/WZA` @ `8191f58`, last commit 2023-12-13,
> cloned at `~/WZA`) — three things this settles:
> - **`--sample_snps` defaults to `0`, i.e. NO CAP** (`general_WZA_script.py:143`,
>   `max_SNP_count = int(1e6)`). So no-cap *is* the upstream default; capping is the opt-in.
>   (`-1` selects a 75th-percentile cap — on our blocks p75 is only 29 (clq0.9 snp) and would
>   resample 25% of all windows, so that option is unusable for LD-block windows anyway.)
> - **Upstream fits the MEAN with deg-2 too** (`general_WZA_script.py:86`), same as the SD.
>   The `np.interp` mean was a **kMate-local invention** introduced with the isotonic branch
>   on 2026-07-21 — it is not upstream behaviour. Smoothing the mean returns us to upstream's
>   class of fit, it does not depart from it.
> - The repo contains **no guidance on negative SDs or NaN p-values** anywhere, and no commit
>   after 2023-12. Booker's 2024-10-04 promise to "add something about this to the GitHub
>   repo" never materialised, so there is no upstream answer to defer to here.
>
> **1. No cap.** The cap has *two* effects and Booker's email only concerned the first:
> (a) it keeps the SD **polynomial** from extrapolating into the sparse tail and turning
> negative; (b) `wza_script.py:204-207` replaces any oversized block's Z with the **mean WZA
> over 100 random `cap`-sized subsets** and records `SNPs = cap`, shrinking |Z|, discarding
> the power of having many variants, and (since `.sample()` is unseeded) making capped blocks
> **irreproducible run-to-run**. Isotonic cannot produce a negative SD, so (a) is moot, and:
> - only **8–15 blocks per class** (≤0.29%) sit past rolling support — the last rolling
>   window is plotted at the *mean* size of its top 50 blocks, so support reaches to
>   3,171 (clq0.9 snp) / 4,044 (clq0.5 snp), not the (n−50)th order statistic;
> - where the tail is measurable, its empirical SD is **below** the last supported value
>   (ratio 0.80 / 0.87) and its mean is **more negative** (−5.1 vs −1.6) → isotonic's flat
>   clip is **conservative**, the safe direction;
> - cap vs no-cap moves almost nothing: **34 vs 34** Bonferroni blocks, 31 shared, 0
>   fabricated `p==0` and 0 NaN in *both* regimes.
>
> **2. Mean fit: `interp` → `isotonic_auto`.** The isotonic branch predicted the mean with
> `np.interp` — piecewise-linear through every rolling point, i.e. no smoothing at all.
> Out-of-sample it is the **worst of six candidates in 12/12** model × def × class cases,
> **37–41% worse** than any smoother, because under the rank transform the mean has almost
> no real N-dependence and `interp` tracks the noise. `isotonic_auto` / `deg2` / `const`
> tie within 1–5% *within support*; `const` is nominally best in 9/12. **`isotonic_auto`
> chosen** on two grounds:
> - the (small but real) trend's *sign is model-dependent* — Spearman ρ(mean, N) is
>   **negative** for kendall (−0.13…−0.23) and binomial (−0.01…−0.14) but **positive** for
>   lfmm (+0.07…+0.09) — so inferring direction from the data beats hardcoding it;
> - **only `isotonic_auto` is bounded past support** (§4b). Predicted mean(Z) at each class's
>   largest block vs what those blocks empirically average:
>
> | class (clq0.9 tile) | isotonic_auto | deg2 | const | empirical tail |
> |---|---|---|---|---|
> | snp | −5.5 | **−289.5** | −0.2 | −10.5 |
> | smallindel | −1.1 | **+16.7** | −0.1 | −5.8 |
> | nonsnp | −1.5 | **+10.0** | −0.1 | −5.8 |
>
> **This is a SECOND fabricated-significance path in the historical deg-2 regime, independent
> of the SD turning over** — a mean predicted ~280 too low makes (Z−mean)/sd enormous → p==0,
> on exactly the same large blocks where the SD was also failing. The two compound. For
> smallindel/nonsnp deg-2 gets the *sign* wrong (+17 vs −5.8), so those blocks could never
> reach significance — total power loss. `isotonic_auto` lands above the empirical tail mean
> in every case → p slightly too large → conservative, the safe direction. Upstream's deg-2
> mean is fine for gene/bp windows with tight SNP counts; it does not transfer to LD blocks
> spanning 2 → 20,365 variants.
>
> **3. deg-7 is disqualified, not merely second.** Booker suggested deg-7 *for his
> dataset*. On our rank-transformed statistic its out-of-sample SD RMSE blows up — worst
> case sv/clq0.5 **17.7 vs isotonic 0.13**; also 3.02 vs 0.41 (snp/clq0.9). deg-2 stays
> competitive only for sv (0.122 vs 0.125, i.e. a tie) and gives negative SDs in 7/8 cases.
> **Isotonic wins or ties everywhere and never goes negative.**
>
> **4. v1's cap rule was wrong** and is retracted: it set the cap at the (n−50)th order
> statistic, which gave **sv a cap of 29** against a max block of 318. That measured where
> rolling *support* thins, not where the *fit* fails — and isotonic's fit fails nowhere.
>
> **Production invocation** (`--sample-snps 0` = no cap; `--mean-fit` added 2026-07-28,
> defaults to `interp` only so old outputs stay reproducible):
> ```
> run_wza.py --sd-fit isotonic --mean-fit isotonic_auto --sample-snps 0 ...
> ```
> **Not yet re-run**: the 189 bio1+multiaxis outputs still carry the capped `interp` regime.
> Regenerating them on tiling + this regime is the outstanding task, and it also removes the
> cap-vs-assignment confound in the four block arms (a09s/a09t/a05s/a05t), where cap was
> varied together with the assignment rule and so could not be attributed.
>
> ### Permutation-null calibration (jobs 35975820 pool / 35976724 site, 20 perms × 4 classes)
> `wza_investigation/perm_null_regime.py` — bio1 site-permuted (signal-free), tiling blocks,
> four regimes. Confirms the regime choice **and** exposes a much larger problem it cannot fix.
>
> **Confirmed:**
> - **`bonf_past_support` = 0.00 in every class × every isotonic regime.** The flat SD
>   extrapolation past the support edge produced *zero* false positives in 20 reps × 4
>   classes. This closes the one risk the fit diagnostics could not.
> - **Cap contributes nothing:** proposed_nocap vs proposed_cap identical to 2 dp in all four
>   classes (snp 8.50/8.50, sv 0.75/0.75, smallindel 8.60/8.60, nonsnp 8.60/8.60).
> - **deg-2/upstream is 2–3× worse** (snp bonf 28.2 vs 8.5, FDR 167.8 vs 58.6) plus ~17 NaN
>   blocks and 0.15–0.55 fabricated p==0 per replicate. Proposed beats old production on FDR
>   in all 4 classes.
>
> **🚨 THE DOMINANT PROBLEM IS CALIBRATION, NOT THE FIT.** Nominal Bonferroni expectation is
> **0.05 blocks**; observed is **8.5** (~170× anti-conservative), and **~58 blocks pass BH
> q<0.05 under a null containing nothing.** Near-identical across all four regimes, so no
> mean/SD fit choice touches it.
>
> **Pseudoreplication is NOT the cause** — the obvious suspect, ruled out. Collapsing to
> flower-weighted site means (n=31 instead of n=352 pools, `--level site`) barely moves it.
>
> **➡ DECISION (T. Bellagio, 2026-07-28): production stays POOL-LEVEL — use all 352 pools.
> Do not collapse to site means.** `--level site` in `perm_null_regime.py` is a *diagnostic
> only*; it was never a proposed production change. The measurement below is what justifies
> the decision: site collapse discards ~90% of the observations (352 → 31) and buys
> essentially no calibration (snp 8.50 → 8.15 null Bonferroni blocks), so it is all cost and
> no benefit. Do not re-open this.
>
> | class | blocks | median blk | bonf pool→site | FDR pool→site |
> |---|---|---|---|---|
> | **sv** | 4,520 | 4 | **0.75 → 0.50** | **1.65 → 0.85** |
> | smallindel | 41,731 | 7 | 8.60 → 8.45 | 52.7 → 47.6 |
> | nonsnp | 41,908 | 7 | 8.60 → 8.50 | 54.3 → 47.8 |
> | snp | 56,846 | 12 | 8.50 → 8.15 | 58.6 → 54.0 |
>
> So the residual is the **WZA aggregation itself**: within-block LD makes the null block-Z
> heavier-tailed than the Normal that `1 - norm.cdf(Z, mean_pred, sd_pred)` assumes. The
> rolling-window correction only standardizes the first two moments; it cannot fix tail shape,
> and Bonferroni at ~57k blocks probes p≈9e-7 (~4.8 SD) — exactly where a fat tail bites.
>
> **Consequences:**
> 1. **The non-SNP candidate-gene lists are heavily contaminated.** At ~8.6 false Bonferroni
>    blocks per single (model × axis) scan, the "428 genes" / "322 genes" unions over 60 scans
>    have an artifact count of the same order as the list. **Do not cite those counts.**
> 2. **SV is the one well-calibrated class** (0.50–0.75 vs 8.5), ~11× better, because SV blocks
>    are tiny (median 4, max 318) so they accumulate little LD-correlated signal. Convenient,
>    given SVs are the point of the project — the SV-specific results are the most trustworthy.
> 3. Miscalibration tracks **block size**, which is also why clq0.5 (coarser) is riskier
>    than clq0.9 — relevant to the still-open blockdef choice.
>
> **Proposed fix (not implemented):** replace the Normal reference with an **empirical null**
> for the block statistic — calibrate observed standardized Z against its site-permutation
> distribution within block-size strata, rather than against N(0,1). The permutation machinery
> already exists. Resolving p to ~9e-7 needs either many more permutations or a tail model
> (e.g. generalized Pareto) fitted to the permutation tail.

Re-run of the phase-1 GrENE-Net GEA (Kendall-τ + LFMM K=16 + binomial → WZA) on
kMate AF, but with **finer LD blocks** (clq0.9 BigLD islands, 58,376 blocks) in
place of the coarse phase-1 hapFIRE blocks (16,674). Motivation: the coarse blocks
were too big (up to 9,158 SNPs), which broke the WZA polynomial correction and
diluted single-haplotype signals (see `../wza_investigation/`). Owner: T. Bellagio.
Started 2026-07-03.

Self-contained subtree: **all code in this folder** (`run_binomial.py` +
`reblock.py` + `compare_clq90.py` + `plot_clq90_manhattan.py` + the two sbatch
drivers), **all outputs under `analysis/grenenet_gea/phase1_replication/results/clq90/`**.

> **Consolidated 2026-07-06:** the session log `SESSION_20260703.md` is merged into this file —
> its Kendall-inflation decision, the 20-axis extension, and the nonSNP-specific candidate-gene
> deliverable are appended below (§3–§5). Original in git (`5cefcfa`).

---

## What changed vs the phase-1-block run (3 asks)

| dimension | phase-1-block run | **clq0.9 run** |
|---|---|---|
| **LD blocks** | hapFIRE, 16,674 (nearest-SNP, tiles genome) | **clq0.9 BigLD, 58,376** (`lib.assign_clq_blocks`, interval, LD islands) |
| **variant classes** | snp · smallindel · sv (separate) | **snp** vs **nonsnp = sv ∪ small-indel together** |
| **generation** | gen9 (last-gen) | gen9 (last-gen) — unchanged |
| **models** | kendall · lfmm(K16) · binomial | same three |
| **climate** | bio1 | bio1 |

## Key methodological decisions (this run)

1. ~~**Blocks don't tile the genome → drop gap variants.**~~ **❌ WRONG — RETRACTED
   2026-07-27, see the CRITICAL banner at the top of this file.** The original text
   read: *"clq0.9 blocks are LD islands with inter-block gaps; ~40 % of records fall
   in a gap (SNP 39 %, nonsnp 40 %) and are dropped — WZA is only defined on a
   window. In-block: SNP 1,205,512 records over 57,587 blocks; nonsnp 413,044 over
   42,677. `lib.assign_clq_blocks` returns `''` for gap variants; `reblock.py` drops
   them."*
   The stated rationale ("WZA is only defined on a window") does not justify
   discarding the record — it justifies *giving it a window*, which is exactly what
   HapFM's `convert_fine_genomewide_breakpoints` does and what this run omitted.
   Cost: 39.7% of SNPs, **49.3% of SVs**, 40.7% of non-SNP records never tested.
   Use `blocks_tiling.py` (0.0% dropped, same LD boundaries) instead.

2. **WZA polynomial: canonical deg-2 (PRIMARY), WITH a per-class SNP cap.**
   clq0.9 blocks are much smaller than the old hapFIRE blocks (SNP max 5,168 vs
   9,158; nonsnp max 963), so the original plan was to drop the deg-7/cap
   workaround entirely. **That was wrong — see §0 BUG + FIX below.** Even at this
   smaller scale, the sparse large-block tail (>~1,150 SNPs snp / >~400 nonsnp)
   still extrapolates the deg-2 SD prediction negative; a pre-existing "safety
   floor" hack silently converted that into fabricated significance instead of the
   NaN it should have produced. Fixed by removing the floor hack and re-adding a
   **class-specific cap (1000 snp / 350 nonsnp)** — smaller than clq0.9's max block
   size but large enough to keep essentially all real large-block signal (only the
   representative-1000/350-SNP-subset Z is used for those blocks, not a truncation
   of the test itself). We keep **deg7cap2000 as a sensitivity** only, unaffected
   by this bug (it already capped at 2000, above where clq0.9's tail turns bad, so
   it never used the floor hack in practice — confirmed by re-checking its NaN
   count is unchanged pre/post fix).

3. **Binomial over-precision → per-variant quasi-binomial (K=16 LF).** The pool
   binomial uses N = flowers×2 per pool, treating each genome as independent; with
   355 pools but only **31 independent climates**, per-record Wald p massively
   underflow, which WZA can only crush to a hard z-cap. **PRODUCTION MODEL
   (2026-07-03, §0B below): per-variant quasi-binomial with K=16 LFMM latent
   factors** — fit `[alt,ref] ~ const + z(bio1) + LF1..LF16` at full N, estimate the
   variant's Pearson dispersion φ = PearsonChi2/df_resid (clip φ≥1), scale the
   climate Wald z → z/√φ, recompute p from N(0,1). This is now
   `phase1_replication/run_quasibinom.py` (`pval_quasi`) — ported verbatim from
   the retired `gea_newpanel/run_quasibinom_latent.py` — reblocked onto our
   clq0.9 partition (blocks_mcf90). The consolidated 3-class driver is
   `run_3class_quasibinom.sbatch`.
   Raw per-record GIF **8.6→2.6–2.8**, no floored/underflowed p (min ~1e-13 vs the
   old ~1e-300 pile-up), and it does NOT flatten (35 snp / 27 nonsnp BH-sig blocks).
   - **Superseded: effective-N = site** (the previous production model, kept in
     `clq90/archive_binomial_effNsite_pre_quasibinom_20260703/`). It set GLM
     `var_weights = 1/(pools per site)` to scale information to ~31 climates. It
     worked but was an ad-hoc deflation whose *raw* per-record GIF was still ~8.6
     (the deflation mostly rescales the tail, not the bulk calibration).
   - **Reconciles the old "quasi-binomial backfires (φ≈0.2, under-dispersion)"
     note:** that earlier test measured φ on the *already-deflated* effective-N
     counts and without the K=16 factors, so little residual dispersion was left to
     scale. On full N with K=16 the residual is genuinely over-dispersed (φ median
     ≈4.4) and the correction deflates as intended. Not a contradiction — different
     count scale.
   - Binomial is the **only** model that re-runs (this is a model change, not a
     block change); kendall/lfmm per-record stats are block-independent and reused.
   - See **§0B** for the full 3-way comparison (quasi vs beta-binomial vs ACER
     effective-N) that led to this choice, and why ACER was rejected.

## Pipeline (how it runs)

```
run_clq90_binomial.sbatch   # effective-N=site binomial, snp + nonsnp, gen9 → clq90/binomial/ (RAW, full)
run_clq90_wza.sbatch        # (dependency=afterok:binomial)
   STEP 1  reblock.py       # reassign clq0.9 block, drop gaps → clq90/wza_in/{model}_{cls}_...csv
                            #   kendall snp reused; kendall nonsnp = concat(sv,smallindel)
                            #   lfmm snp/nonsnp reused (lfmm_nonsnp already existed)
                            #   binomial snp/nonsnp from the effective-N run
   STEP 2  run_wza.py       # deg2 (primary) + deg7cap2000 (sensitivity) → clq90/wza/
   STEP 3  NaN/neg-SD audit per regime
compare_clq90.py            # BH/Bonferroni, cross-class & cross-model recurrence, CAM5 rank → clq90/compare_clq90.csv
plot_clq90_manhattan.py     # 3×2 (model×class) Manhattan (basic env) → clq90/manhattan_clq90_deg2.png
```

Reblocked WZA inputs live in a **dedicated `clq90/wza_in/`** dir (not the raw
source dirs) so the raw binomial is never clobbered and other r² thresholds can be
re-blocked from the raw per-record files.

## §0 BUG + FIX (2026-07-03) — the deg-2 "0 NaN" claim was masking a real fitting bug

**Symptom:** T. Bellagio flagged the deg-2 Manhattan plots as "weird" — a handful of
points per panel pinned exactly at −log10 p = 300, disconnected from everything
else. That pile was clipped-for-plotting float64 underflow (`Z_pVal.clip(1e-300)`
in `plot_clq90_manhattan.py`), which was already noted as "cosmetic" in the first
pass of this doc — **that framing was wrong.**

**Root cause (same failure mode as `../wza_investigation/RESULTS.md`, which this
run's original write-up should have cross-checked and didn't):** `wza_script.py`'s
deg-2 SNP-number correction fits `Z` vs block-size (`SNPs`) with a rolling-window
mean/SD + degree-2 polyfit, then `Z_pVal = 1 - norm.cdf(Z, mean_pred, sd_pred)`. For
the **sparse tail of large blocks** — clq0.9 SNP-class blocks with >~1,150 SNPs
(22/57,587), nonsnp with >~400 SNPs (7-11/42,677) — the quadratic **extrapolates
past its rolling-window support and predicts a negative SD**. A pre-existing
"safety floor" hack (`wza_script.py`, meant only to avoid a NaN crash) clipped that
negative SD to **the smallest positive SD predicted anywhere on the curve** — a
value from a totally different block-size regime, not a real estimate for that
block. Result: fabricated, often-literal-zero p-values **independent of the true
Z** — e.g. block Chr5_4452 (SNP class, N=1726, raw **Z=0.23**, i.e. an
unremarkable/null WZA score) was assigned **Z_pVal=0.0** (nominally the single most
significant block in the genome) purely from this extrapolation artifact.
Confirmed present in all 6 model×class deg-2 outputs (7–24 affected blocks each).

**This was not a cosmetic tail — it drove the headline "reproducible core."**
Cross-checking the bad-block list against the old §"Cross-model / cross-class
recurrence" table below: **7 of the original 10** cross-model-significant SNP
blocks (Chr1_7103, Chr2_1283, Chr2_1623, Chr2_2819, Chr3_6064, Chr3_6144, Chr4_1721)
were exactly the negative-SD-floor blocks — the "most reproducible" claim was
mostly the bug, not signal.

**Fix (matches `wza_investigation`'s own bottom line: "cap SNPs/window... drop the
SD safety-floor hack"), applied here:**
1. **Removed the SD safety-floor hack from `wza_script.py`.** Negative predicted SD
   now correctly yields `Z_pVal = NaN` (excluded from FDR) instead of a fabricated
   value — matches `wza_core.py`'s (the investigation's) no-floor default.
2. **Re-added a per-class SNP cap** (`--sample_snps`) so no block is ever evaluated
   outside the polynomial's well-supported domain: **1,000 for snp, 350 for
   nonsnp**. Chosen empirically (`wza_core.py` cap×degree sweep, cross-checked
   against all 3 models): the largest cap that gives **0 negative-SD blocks in
   kendall, lfmm, AND binomial simultaneously** (snp breaks again at cap≥1,100;
   nonsnp at cap≥400 — note production's unseeded `--sample_snps` resampling landed
   slightly worse than the seeded sweep at exactly 400, so nonsnp settled at 350 for
   real margin). This is a **class-specific, data-driven** cap, not a reused
   `deg7cap2000`-style round number — clq0.9 blocks are much smaller than the old
   hapFIRE blocks (max 5,168 vs 9,158) so the old cap=2000 would NOT have been low
   enough to avoid the tail here.
3. Re-ran **STEP 2 (WZA) only** for all 6 model×class combos (STEP 1 reblocking and
   the binomial effective-N per-record stats are untouched — the bug was purely in
   the correction stage, not upstream). Verified **0 negative-SD / 0 NaN in all 6
   corrected outputs** by independently recomputing the rolling-fit sign, not just
   trusting the script's own report.
4. Regenerated `compare_clq90.csv`, `manhattan_clq90_deg2.png`, and
   `significant_genes_clq90_gen9_bio1_deg2.csv` from the corrected WZA outputs.
   Old (buggy) artifacts kept for the record in `clq90/wza/archive_deg2_prefloor_bug/`
   and `clq90/archive_prefloor_bug/`.

**Sanity check that the fix didn't disturb real small-block signal:** CAM5
(Chr2_4332, ~13–83 SNPs, far below either cap) is essentially unchanged — kendall
snp p=8.06e-5/q=0.066 now vs 8.5e-5/q=0.059 before (§ below) — as expected, since
capping/floor-removal only touches the large-block tail.

## ★ RESULTS (gen9, bio1, CORRECTED) ★  (`clq90/compare_clq90.csv`, `manhattan_clq90_deg2.png`)

> **[SUPERSEDED 2026-07-27 — see the banner at the top of this file.]** Everything in
> this section is deg-2 (2-class), the regime later shown to fabricate significance on
> sparse large blocks. Under the current isotonic + 3-class production output
> (`clq90/wza/wza_*_gen9_bio1_isotonic.csv`), **CAM5 does not reach BH-FDR in any
> model×class (best q=0.346)** and **the cross-model "reproducible core" below does
> not survive** (0 blocks BH-sig in all 3 models, for any class). Kept for the record;
> do not cite the numbers below as current.

### deg-2 (PRIMARY, capped 1000/350, no floor hack) — BH q<0.05 block counts
| model | snp (57,587 blk) | nonsnp (42,677 blk) |
|---|---|---|
| kendall  | 63 | 15 |
| lfmm     | 66 | 24 |
| binomial | 50 | 33 |

(For reference, the pre-fix — WRONG — counts were snp 74/79/58, nonsnp 18/27/38;
the corrected counts are lower everywhere, consistent with removing fabricated
large-block hits. deg7cap2000 sensitivity numbers are unaffected by this bug — it
already capped at 2000 and never hit the floor — and are unchanged: snp 17/16/1,
nonsnp 5/8/2.)

### CAM5 (AT2G27030) — gene span splits across 6 clq0.9 blocks (Chr2_4327…4332)
The phase-1 3′-end signal localises to **Chr2_4332** in every model (finer blocks
correctly isolate it). Best per (model,class), deg-2 (corrected):

| model | class | CAM5 block | rank | p | BH-q |
|---|---|---|---|---|---|
| kendall  | snp | Chr2_4332 | 70/57,587  | 8.1e-5 | **0.066** |
| lfmm     | snp | Chr2_4332 | 104/57,587 | 1.5e-4 | **0.083** |
| binomial | snp | Chr2_4332 | 190/57,587 | 6.6e-4 | 0.198 |
| (nonsnp best: lfmm rank 122, q=0.148 — CAM5 is a SNP signal, weaker in nonsnp) |

**Reading (unchanged from before the fix — CAM5 lives in a small block, was never
affected by the bug):** under the finer blocks + honest canonical deg-2, CAM5 is a
strong top-0.15–0.4 % block but sits **just below BH-FDR** (q≈0.07–0.08 in
kendall/lfmm snp). Cleaner than the coarse-block run: 57k tests (vs 16k) is a
heavier multiple-testing burden, and deg-2 doesn't over-fit the tail the way the
phase-1 deg-7 did (which had inflated CAM5 to q=0.017 / p≈3e-8 — shown to be tail
over-fit in `../wza_investigation/`, honest poly-free p≈1.6e-4, matching here).

### Cross-model / cross-class recurrence (deg-2, CORRECTED — the real robust core)
- **SNP blocks BH-sig in ALL 3 models (4):** Chr1_2343, Chr1_6539, Chr3_6144,
  Chr4_6307. (Down from the pre-fix 10 — Chr1_7103, Chr2_1283, Chr2_1623, Chr2_2819,
  Chr3_6064, Chr4_1721 dropped out; those were the negative-SD-floor artifacts.)
- **nonsnp blocks BH-sig in ALL 3 models (2):** Chr1_2343, Chr4_6307. (Down from 4 —
  Chr2_2819 and Chr3_6144 dropped from the ALL-3-models nonsnp set, though Chr3_6144
  remains BH-sig in binomial+kendall+lfmm at the SNP level and 5/6 combos overall.)
- **Chr1_2343 (GAPC2) and Chr4_6307 (CRK13-16 cluster)** are now the only blocks
  significant in every model AND every class — the actually-reproducible core.
  Chr3_6144 (AT3G30320/UMAMIT32) is close behind (5/6 combos). See
  `significant_genes_clq90_gen9_bio1_deg2.csv` for the full gene-annotated table
  (176 BH-sig blocks / 199 genes, union across all 6 combos).

> Manhattan note: the −log10 p = 300 pile-up from the bug is gone (verified: 0
> negative-SD blocks in all 6 corrected outputs). Any block still near the 1e-300
> plotting floor now reflects a genuinely well-fit extreme tail probability, not an
> artifact — e.g. capped-at-1000 blocks with raw Z in the teens-to-30s can
> legitimately underflow float64 there.

## Jobs (done)
- 35500487 `clq90_binom` — raw effective-N binomial (snp+nonsnp). COMPLETED.
- 35508732 `clq90_wza` — reblock + WZA deg2/deg7cap2000 + NaN audit. COMPLETED (26 min).
  (35500488 was an earlier attempt that hit a relative-path bug in run_wza.py, since
   fixed by abspath-ing the output; 35500485 was the pre-cleanup attempt, cancelled.)
  **This job's deg-2 output was later found buggy — see §0 above.**
- (interactive, 2026-07-03) SD-floor bug fix: cap×degree sweep (`wza_core.py`,
  reused from `../wza_investigation/`) + corrected WZA re-run (STEP 2 only, all 6
  model×class combos, cap 1000/350) + `compare_clq90.py` + `plot_clq90_manhattan.py`
  + `build_significant_genes_clq90.py` re-run. Ran interactively on a compute node
  (not sbatch — each combo ~1-2 min with the resampling cap).

## NEXT / open
- ~~Annotate the recurrent blocks → genes~~ **DONE** —
  `build_significant_genes_clq90.py` → `significant_genes_clq90_gen9_bio1_deg2.csv`
  (176 blocks / 199 genes on the corrected regime). Top hits: Chr1_2343 (GAPC2),
  Chr4_6307 (CRK13/14/15/16), Chr3_6144 (AT3G30320, UMAMIT32/AT3G30340).
- Fold the same floor-removal + per-class cap fix into `wza_script.py`'s use
  elsewhere in this repo (e.g. the phase1-block `deg7cap2000` sensitivity path was
  unaffected here only because it already capped at 2000 — any other caller of
  `wza_script.py`'s deg-2 path without a cap should be audited).
- Optional: other r² thresholds (clq0.5/0.7 tsvs exist in `blocks_mcf90/`).
- Optional: bio2–19 (only bio1 here, as in the phase-1-block run). **DONE — see §4 below (20 axes).**

---

# §3–§5 (merged from SESSION_20260703.md) — Kendall decision, multi-axis, candidate genes

## §3. Kendall-τ inflation → left AS-IS (no honest correction preserves signal)
Pool-level Kendall is the worst-inflated model: GIF ~9.2, **52% of the genome at p<0.05** (1e-16
floor). Corrections tested (`clq90/kendall_fix_test/`, `notebooks/kendall_fix_compare.ipynb`):

| method | GIF | max −log10p | verdict |
|---|---|---|---|
| pool-level (current) | 9.2 | 16 (floor) | inflated |
| site-collapse (31 sites) parametric | 2.26 | 2.9 | residual structure, coarse |
| MSR structure-preserving null | 1.3 | 4.0 (perm floor) | flat-ceiling |
| pool + genomic control | 1.0 | 2.4 | flattens |
| site + genomic control | 1.0 | 2.9 | flattens |

**Every honest correction flattens Kendall** — the tall peaks *were* the pseudoreplication (355
fake-independent pools; only ~31 independent climate sites). Decision: keep Kendall as the
raw/structure-uncorrected phase-1 reference; **LFMM is the model to trust for calibrated climate
signal** (structure-corrected via latent factors, keeps resolution). No production change to Kendall.

## §4. Multi-axis extension — all 20 axes (bio1–19 + PC1)
Ran the full clq0.9 pipeline (kendall + lfmm + quasi-binomial → reblock → WZA deg-2, caps 1000/350)
for **20 climate axes × 3 models × 2 classes = 120 WZA outputs**.
- PC1 = PC1 of 19 standardized bioclim (47% var; r=+0.83 temp, −0.81 precip; warm-dry↔cool-wet),
  added as `pc1` col to `class_matrices/gen9.pools.csv`.
- Code `multiaxis/` (axes.sh, ma_{kendall,lfmm,quasibinom,wza}.sbatch, reblock_multiaxis.py,
  run_lfmm_nogif.R); outputs `clq90/multiaxis/`. Jobs COMPLETED (40/40 each): kendall 35517083,
  quasibinom 35517085, lfmm-no-gif 35517329, wza 35517330. Summary `multiaxis/multiaxis_summary.csv`.
- **LFMM `calibrate="gif"` DROPPED (per user):** gif divides by λ=median(z²)/0.456, only valid when
  λ>1; where K=16 over-corrects (λ<1, e.g. bio19 λ=0.931 snp) it would *inflate*. `run_lfmm_nogif.R`
  writes RAW p + logs λ per axis (range 0.931–2.536) → raw LFMM p is inflated on high-λ axes; λ
  recorded if a guarded recalibration is wanted later.

## §5. Deliverables — cross-axis overlap + nonSNP-specific candidate genes

> **[SUPERSEDED 2026-07-27 — see the top-of-file banner.]** Both (a) and (b) below were
> built on deg-2 WZA + the retired pooled 2-class (snp/nonsnp) split. Corrected
> (isotonic + 3-class) numbers: union FDR = **327 blocks → 322 genes** (nonsnp-pooled
> alone: 217→212, vs the 437/548 below); only 133/547 (24%) of the gene list survives;
> HSFA2/COR15A/VRN2 (the genes driving the "heat-stress dominates" read) are **absent**
> from the corrected list, HSBP1 is the one survivor. No corrected functional-theme
> read has been done yet — **do not cite "heat-stress dominates" as current.**
> Current scripts: `class_peaks_overlap_table.py` (a), `class_specific()` in
> `_build_snp_vs_nonsnp_peaks_nb.py` (b). Text below kept for the record only.

**(a) Cross-axis SNP↔nonSNP hit-block overlap** (`multiaxis/overlap_table.py` →
`overlap_snp_nonsnp_by_axis.csv`): consistent across axes — **SNP hits ≈ 2× nonSNP, Jaccard
~0.22–0.32, 23–59 nonSNP-only blocks/axis** (most nonSNP signal not shared with SNP).

**(b) nonSNP-specific hit blocks → genes** (`multiaxis/nonsnp_specific_genes.py`, Ensembl Plants
REST → `nonsnp_specific_genes{,_table}.csv`, notebook `notebooks/nonsnp_specific_genes.ipynb`):
**437 blocks** BH-sig in nonSNP not SNP (305 never a SNP hit anywhere), **548 genes**.
- **Functional read — HEAT-STRESS dominates; flowering/circadian essentially ABSENT.**
  Heat (11 genes, in CALIBRATED models not just Kendall): **HSBP** (10 axes, all 3 models),
  **HSFA2** (master thermotolerance TF), **HSP20-like small-HSP cluster** Chr1_9404 (5 axes, all 3),
  HEAT-repeat Chr5_2542, HSP17.4, DNAJ. Drought/ABA (5): **RAS1** (5 axes), ERD, senescence.
  Cold (1): COR15A. Flowering/circadian: only VRN2 (1 axis) — no FT/FLC/CO/GI/CCA1/TOC1/PRR/ELF/PIF4/PHYB.
- Takeaway: the SV/indel-borne climate signal points at heat-shock/stress-tolerance machinery, NOT
  the flowering-time pathway that dominates classic Arabidopsis climate GWAS, across temp + precip axes.

Open (from the session log): nothing git-committed at the time (now snapshotted); LFMM raw p inflated
on high-λ axes (λ logged); optional cross-axis overlap heatmap + calibrated-only nonSNP gene table.

## §6. Corrected functional-theme read (2026-07-27) — "heat-stress dominates" DOES NOT hold up

Re-ran the symbol/description/functional-category annotation (`multiaxis/candidate_genes_corrected.py`,
identical Ensembl Plants REST lookup + identical `CATS`/`CURATED` keyword scheme as the retired
`nonsnp_specific_genes.py`/`_build_nonsnp_genes_nb.py`, so directly comparable) on the current
isotonic + 3-class union of class-specific blocks (`{nonsnp,sv,smallindel}_specific_peaks_{bonf,fdr}.csv`):
**345 blocks → 248 gene-bearing blocks → 336 unique genes** (union over classes/tiers; close to,
not identical to, the 327/322 FDR-only count in the banner above because this pass also folds in
the Bonferroni tier and uses per-gene block dedup).

| category | genes | vs. old (deg-2/2-class) claim |
|---|---|---|
| heat | **5** | 11 |
| drought_ABA | **2** | 5 |
| cold | **1** | 1 (different gene) |
| flower_circ | **0** | 1 (VRN2) |
| other (no keyword match) | **328 (97.6%)** | — |

Only **HSBP/AT4G15802** survives from the old heat list (HSFA2, the HSP20 cluster, RAS1, COR15A,
VRN2 are all gone from the corrected candidate set); the one cold hit is a different gene
(COR413-PM2, not COR15A). **8/336 genes (2.4%) landing in any stress category reads as background
rate, not a theme** — consistent with the GWAS-thread GO-enrichment finding (`VAREXP_SELECTION_HANDOFF.md`)
that named stress genes turn up at background frequency in these non-SNP candidate lists generally.
**Conclusion: do not cite a functional theme (heat-stress or otherwise) for the multiaxis climate-GEA
non-SNP-specific gene list** — the corrected list has no coherent signal beyond individual anecdotes.
Full annotated table: `multiaxis/candidate_genes_corrected.csv`.

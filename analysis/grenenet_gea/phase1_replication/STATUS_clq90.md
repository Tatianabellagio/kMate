# Status — phase-1 GEA replication on **clq0.9 (r²≥0.9) LD blocks**

> ## ⚠ CORRECTION (2026-07-03) — the original deg-2 "0 NaN, workaround not needed"
> ## result below was WRONG. A real WZA fitting bug inflated ~10-25 large blocks
> ## per output to spurious p≈0. Fixed; corrected results replace the tables below.
> See **"§0 BUG + FIX"** right after this box for the full story; all numbers
> further down in this file are already updated to the corrected regime.

Re-run of the phase-1 GrENE-Net GEA (Kendall-τ + LFMM K=16 + binomial → WZA) on
kMate AF, but with **finer LD blocks** (clq0.9 BigLD islands, 58,376 blocks) in
place of the coarse phase-1 hapFIRE blocks (16,674). Motivation: the coarse blocks
were too big (up to 9,158 SNPs), which broke the WZA polynomial correction and
diluted single-haplotype signals (see `../wza_investigation/`). Owner: T. Bellagio.
Started 2026-07-03.

Self-contained subtree: **all code in this folder** (`run_binomial.py` +
`reblock.py` + `compare_clq90.py` + `plot_clq90_manhattan.py` + the two sbatch
drivers), **all outputs under `results/grenenet_gea/phase1_replication/clq90/`**.

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

1. **Blocks don't tile the genome → drop gap variants.** clq0.9 blocks are LD
   islands with inter-block gaps; **~40 % of records fall in a gap** (SNP 39 %,
   nonsnp 40 %) and are dropped — WZA is only defined on a window. In-block: SNP
   1,205,512 records over 57,587 blocks; nonsnp 413,044 over 42,677.
   `lib.assign_clq_blocks` returns `''` for gap variants; `reblock.py` drops them.

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
   climate Wald z → z/√φ, recompute p from N(0,1). This is `gea_newpanel/
   run_quasibinom_latent.py` (`pval_quasi`), reblocked onto our clq0.9 partition.
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

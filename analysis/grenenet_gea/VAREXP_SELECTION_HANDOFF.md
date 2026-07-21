# SNP vs non-SNP variance partition of ecotype selection — result + handoff

> ## ✅ NUMBERS ARE CURRENT — post-Kf_w / `--unit chrom`, regenerated & verified 2026-07-21
> Every quantity in this doc has been updated to the **post-fix founder h** (full-panel Kf_w
> correction + `--unit chrom`). The full pipeline was re-run end-to-end on 2026-07-21
> (`build_sample_h_cache` → `ecotype_fitness` → `build_selection_trait` → `build_class_grms` →
> `varexp_selection` → `class_split_gwas`) and diffed against the previously-committed outputs:
> `selection_s_matrix.npz`, `class_grms.npz`, `varexp.csv`, and `class_gwas_summary.json` all
> reproduced **bit-identically** (they were already regenerated post-fix on 2026-07-07, committed
> `3c34ec8`), confirming the numbers below are stable. The founder-collapse is fixed: **all 231 of
> 231 founders are retained** (the old `p0>1e-3` drop-floor that dropped ~19 twin-absorbed founders
> is now a floored no-op; 9977 gets a defined near-zero h instead of being absorbed).
>
> The one genuinely-stale artifact was the *Tier-2* ecotype-fitness GWAS `ecotype_fitness/gwas/
> gwas_z.npz` (dated 2026-07-02, built on the pre-fix phenotype and never rebuilt when the
> phenotype was corrected on 07-07). It was regenerated 2026-07-21: the generalist axis reproduces
> (corr 0.92–1.0), but the climate-zone axes (cold/mid/hot) drifted substantially (corr 0.27–0.57)
> — the Kf_w fix finally propagating into Tier-2. That is a *different* analysis from this doc
> (see `ecotype_fitness.py`), noted here only so the staleness is on record.

Session date: 2026-07-03 (original); numbers current as of **2026-07-21**. Branch `add-kmate`.
Env: `kmate` (`source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate`). Always
`hostname` first — compute only on `n*.savio*`, never `ln00X` (hook-enforced). Heavy steps
(`build_class_grms`, `class_split_gwas`) run via **sbatch** (`run_class_grms.sbatch`,
`run_class_split_gwas.sbatch`), not interactively — the 16 GB interactive-session cgroup OOM-kills
`class_split_gwas` (needs ~13 GB+).

> **Consolidated 2026-07-06:** the same-day session log `SESSION_2026-07-03_CLASS_SPLIT_GWAS.md`
> (class-split GWAS peaks, overlap tables, Manhattan PNGs, non-SNP-only candidate genes) is
> merged into this file — see "Downstream" section below. Original preserved in git (commit
> `5cefcfa`).

**Bottom line (both threads agree):** at every resolution — genome-wide kinship, SNP-untagged
kinship, per-marker GWAS peaks, per-site and multi-site scans, all bioclim axes — SNPs and
non-SNP markers (indels+SVs) tell essentially the same story. **The non-SNP layer adds almost
nothing to the polygenic selection signal**; what little it flags uniquely is scattered,
marginal, and mostly sub-Bonferroni. (Distinct from the per-variant *temporal* SV insertion
signal in `SV_TEMPORAL_PURGING_SUMMARY.md` — different question and unit.)

## The question (reframed this session)

Prior GWAS work asked "which individual SVs/haploblocks are under selection?" → **null**
(SVs = passengers; block ×2.07 was a unit artifact; haplotype enrichment collapsed under
rotation nulls). User reframed to the phase-1-vs-kMate story:

> **How much of the ecotype selection-response variance do SNPs vs non-SNP variation
> (indels+SVs) explain, and how much do we GAIN by adding the non-SNP layer?**

Design decisions (all locked with the user):
- **Estimand**: report BOTH marginal (each class alone) AND conditional gain (non-SNP beyond SNP).
- **Class split**: 2-way SNP vs non-SNP primary; 3-way SNP/indel/SV secondary.
- **Trait flavour**: relative (compositional selection) primary; census secondary.

## Trait decision (important — do NOT use raw h / Δh)

Checked distributions empirically. Raw founder freq `h` / `Δh` are unusable for GWAS:
skew **+12**, excess kurtosis **+237**, and ~10 founders carry **52%** of the across-founder
variance (a "231-founder" GWAS is really ~10 clades). The fix = **selection coefficient in
log-odds space**: `s = logit-slope of h over gens 0→3` → skew +1, variance spread across founders.

**Locked trait** = per-founder, per-site **logit-slope selection coefficient** `s[site×founder]`:
- built by `build_selection_trait.py` → `analysis/grenenet_gea/varexp/selection_s_matrix.npz`
  (`S[30×231]`, `sites`, `founders`, `bio1`, `p0`, `analyzable`(bool), `n_plots`, `freq_last`).
- GLOBAL-mode chrom-averaged `h` (reused from `ecotype_fitness/sample_global_h.npz`, GLOBAL_MODE_DECISION).
- founding reference `p0` = mean over 8 SEEDMIX reps (ESTIMATED, not forced uniform 1/231 — the
  twin-absorption identifiability bias cancels in the slope; verified seedmix averaging is correct,
  mean is exactly 1/231 by closure, spread is real & reproducible across all 8 reps).
- **presence filter**: founders are kept via a small additive **p0 floor** (`P0_FLOOR`), NOT a
  drop threshold, so **all 231 of 231 founders are analyzable** (post-fix, current). The old
  `p0>1e-3` drop-floor kept only 212 (dropping ~19 twin-absorbed founders like 9977 at h≈1e-15);
  after the full-panel Kf_w fix the collapse is gone (9977 now h≈1.5e-4) and the floor is a no-op —
  every founder gets a defined near-zero slope instead of being excluded.
- **NO reliability/cross-chrom-SD weighting** (user decision: chrom-averaging already regularizes).
- raw `s` primary; **RINT** (rank-inverse-normal per site) as sensitivity.

## What was built

1. `build_selection_trait.py` → `selection_s_matrix.npz`. **DONE.** (trait skew −0.88, excess kurtosis +1.4, 30 sites, 231 founders — well-behaved, variance spread across founders.)
2. `build_class_grms.py` + `run_class_grms.sbatch` → `analysis/grenenet_gea/varexp/class_grms.npz`.
   **DONE** (job 35513509). Per-class founder GRMs (231×231), standardized `Z=(g−p)/√(p(1−p))`
   — **0/1 founder presence = haploid coding** (inbred lines; NOT the diploid 2p form). Markers:
   MAC≥5, call≥0.9. `K_snp` (1.75M), `K_indel` (512k), `K_sv` (12.8k), `K_nonsnp` (525k),
   `K_all`, `K_snp_matched` (SNPs subsampled to non-SNP count & MAC spectrum, 525k). Aligns to
   trait founders by ID in the modeling step (var_pa founder ORDER differs from seedmix order).
3. `varexp_selection.py` → `varexp.csv`, `varexp_persite.npz`, `varexp_meta.json`. **DONE.**
   Estimators: (1) marginal single-GRM AI-REML h²+SE+LRT, (2) joint 2-GRM REML partition + LRT for
   non-SNP-beyond-SNP, (3) GBLUP repeated k-fold CV predictive R² & gain (NFOLD=6 NREP=10, SEED=0
   → deterministic, reproduces bit-identically). Traits = 4 aggregate axes (w_global + bio1-tercile
   zones) × {raw,RINT} + 30 per-site (raw). Light (all modeling on 231×231 from cache, ~15 min).

## KEY FINDING (genome-wide GRM) — non-SNP is redundant with SNP

On `w_global` (raw): marginal h² = **0.672 (snp) vs 0.673 (nonsnp) vs 0.668 (matched) vs 0.646
(sv) vs 0.673 (all)** — indistinguishable across classes; joint 2-GRM **LRT p_add_nonsnp = 0.5**
(n.s.); CV **gain_nonsnp = −0.003** (R² snp 0.368, nonsnp 0.362, both 0.366). RINT is the same
story (h² snp 0.644 / nonsnp 0.643, gain −0.003). Driver: **corr(K_snp, K_nonsnp) = 0.998.** The
selection trait is strongly kinship-predictable (winners are whole clades), but SNP and non-SNP GRMs
are the *same genealogy*, so non-SNP adds nothing at the polygenic level. Consistent with
SVs=passengers + SNPs tagging ~99.7% of non-SNP markers at r²≥0.2.
(Numbers dropped from the pre-fix run's h²≈0.938 to ≈0.67 with the corrected `--unit chrom` /
Kf_w founder h and all 231 founders retained — the class-*equality* conclusion is unchanged.)

**Caveat = the reframe:** a genome-wide GRM answers "is non-SNP *relatedness* different from SNP
relatedness?" (no, by construction) — NOT "what does the non-SNP layer add that SNPs can't see?"
Aggregate kinship washes out the local, SNP-untagged signal that is the actual kMate-gain question.

## NEXT STEP — DONE (2026-07-03): SNP-untagged non-SNP GRM — still no gain

Built `K` from ONLY indels/SVs that SNPs cannot tag, and partitioned selection variance against
`K_snp`. Answer: **still no gain — the null survives isolating the most SNP-blind possible
markers.** This preempts the obvious rebuttal to the earlier genome-wide result ("of course
K_snp and K_nonsnp look alike, they're both just genome-wide relatedness") by showing the null
even when non-SNP markers are restricted to ones with essentially zero local LD to any SNP.

**What was built:**
1. `build_nonsnp_tagging.py` (new — extends `build_sv_snp_ld.py`'s SV-only r² scan to the exact
   marker set that feeds `K_nonsnp`, i.e. indel+SV, MAC≥5, call≥90%, same filter as
   `build_class_grms.py`, so results are an exact subset by column index, not a position-matched
   approximation). For each of the 525,043 non-SNP markers, max founder-genotype r² against any
   panel SNP within ±50kb (SNP side unfiltered — max chance to tag). Run per-chrom directly on
   the compute node (~35–120s/chrom, no sbatch needed). Output:
   `analysis/grenenet_gea/varexp/nonsnp_tagging_chr{1..5}.npz` (col_idx, pos, cls, best_r2,
   best_snp_pos, n_snp_window).
   **Finding en route:** tagging is far tighter than expected — even at a lenient r²≥0.2 bar,
   **99.7% of non-SNP markers genome-wide are tagged.** Only 1,383/525,043 (0.26%) have NO
   panel SNP within ±50kb reaching r²≥0.2. At r²≥0.5: 91.8% tagged (42,140 untagged, 8%). At
   r²≥0.8 (matches the old SV-only headline): ~73% tagged.
2. `build_untagged_grm.py` — builds two untagged GRMs at both thresholds (reuses
   `build_class_grms._ZZt`): `K_untagged_r02` (1,383 markers, 1,360 indel + 23 SV,
   corr with K_snp = **0.364** — genuinely decorrelated) and `K_untagged_r05` (42,140 markers,
   1,068 SV, corr with K_snp = 0.950 — still fairly tied to it). ->
   `analysis/grenenet_gea/varexp/untagged_grms.npz`. (These GRMs are panel-genotype-derived, so
   independent of the Kf_w founder-h fix; the h²/gain numbers below are on the post-fix trait.)
3. `varexp_untagged.py` — reruns marginal h², joint 2-GRM REML+LRT, GBLUP CV gain (imports the
   estimators straight from `varexp_selection.py`) for `K_snp` + `K_untagged_{r02,r05}`, across
   the 4 aggregate axes (w_global + 3 bio1-tercile zones) × {raw, RINT}. ->
   `analysis/grenenet_gea/varexp/varexp_untagged.csv`. (Minor bug, harmless: the summary printer
   crashes on `df.thresh` — pandas resolves `.thresh` to a builtin, not the column; fixed to use
   `df["thresh"]`, but the CSV itself is already complete before the crash — no rerun needed.)

**Result (w_global, raw, r²<0.2 — the strict/primary threshold):** marginal h²(untagged) =
**0.708** (SE 0.115, still large) — because even locally SNP-invisible markers still reconstruct
the SAME founder clade structure (kinship is a genome-wide, not local, property; unlinked
markers anywhere still tag the genealogy). But conditional on `K_snp`: joint LRT **p_add = 0.10**
(n.s.), and CV **gain = +0.010** (r² snp 0.368 vs both 0.378) — a trivial ~1%-of-variance bump,
not significant. **Robust to the null across all 4 axes × 2 transforms × 2 thresholds (16 cells):**
CV gain is uniformly small-positive (+0.003 to +0.010, max +0.010 at w_hot/raw/r02) and every
joint-LRT p_add is **≥ 0.056** (smallest at w_global/RINT/r02, still n.s.; r05 cells run 0.28–0.5).
No cell reaches significance — the "no local, SNP-independent non-SNP signal" conclusion holds,
though the gains sit consistently just-positive rather than sign-flipping as in the pre-fix run.

**Interpretation:** the "no kMate gain" conclusion is now a genuinely strong, non-tautological
result. It's not that K_nonsnp merely LOOKS like K_snp (the genome-wide confound) — even markers
SNPs cannot locally tag still only recover the same founder genealogy, and genealogy is already
~67–71% predictive (marginal h²) of the selection trait via SNPs alone. There is no local, SNP-independent
non-SNP signal riding on top of clade structure. This is the strongest form of "SVs/indels are
passengers, not an independent selection currency" reached in this whole investigation thread —
consistent with [[project-gea-sv-selection-currency]]'s hap-cluster-level null and the parallel
SuSiE fine-mapping session's independent passenger conclusion.

## Peaks + per-marker GWAS comparison (2026-07-03, same session) — SNP-only vs non-SNP-only scan

Distinct from the kinship/variance-partition question above: does an actual **single-marker**
GWAS (not a GRM) run separately on SNP-only vs non-SNP-only markers find the same **peaks**, at
both per-site and multi-site (JOINT) resolution? New script `class_split_gwas.py`:

- Same locked trait (`s`, 231 analyzable founders, 30 sites), rank-normalized per site.
- **ONE shared LOCO kinship correction per chromosome** (built from ALL classes pooled,
  MAC≥12/call≥90%, total-minus-chromosome trick — no repeated whole-genome loads) — deliberately
  shared across both scans so any difference in results is attributable only to the TEST-marker
  class, not to also varying the correction (justified: K_snp/K_nonsnp are near-identical per
  the variance-partition result above).
- Test markers MAC≥5/call≥90% (matches `class_grms.py`), split SNP (1,752,846 markers) vs
  non-SNP (525,043, indel+SV pooled) vs strict SV-only (12,789). Per (chrom, class): one
  eigendecomposition + one U.T-rotation of the marker block, reused across all 30 sites (only
  REML-delta + GLS re-runs per site) — full genome × 30 sites in ~5 min, run via
  `run_class_split_gwas.sbatch` (the interactive-session cgroup OOM-kills it at ~13 GB).
- Per-site z genomic-control calibrated, stacked into per-class Z[M,30], run through the same
  Bolormaa JOINT(30df)/GLOBAL(1df)/CLIMATE(1df) meta as `founder_gwas_multisite.py`.
- **Numerical gotcha hit and fixed:** ~7–10 markers per chromosome (out of ~450k) are
  near-FIXED (nearly all 231 founders carry the allele, not rare — MAC treats both tails
  symmetrically so this looks like a rare-variant floor case but isn't). A near-fixed genotype
  vector is ~collinear with the model intercept after LOCO rotation, so `Saa*Sbb - Sab^2`
  computes as pure float64 rounding noise (~1e-16 relative, i.e. exactly machine epsilon — not
  a small real value) instead of the true 0. Dividing by that near-zero noise produced |z|~1e19
  for those ~7 markers, which alone corrupted `corr(Z)` (and hence Bolormaa for EVERY marker —
  lambda_JOINT collapsed to 0.000, "Bonferroni-significant" for 38% of markers). Fixed by
  masking markers with relative det ≤1e-8 (clean 10-order-of-magnitude gap between the ~7
  degenerate markers at 1e-16 and the smallest legitimate marker at ~1.4e-6) as untestable
  before the cross-site meta — standard GWAS practice (monomorphic/near-collinear-after-
  conditioning exclusion), negligible marker loss (<0.002%).

**Result — variance explained is nearly identical two independent ways:**
mean per-site genomic-control λ = **1.015 (SNP) vs 1.008 (non-SNP) vs 1.004 (SV)** — essentially
the same mild inflation, from real single-marker testing (not REML), corroborating the earlier GRM
h² result (median per-site h²: 0.862 snp vs 0.855 nonsnp; median CV R²: 0.459 vs 0.454, from
`varexp_persite.npz`). The cross-site-meta λ_JOINT is 1.249 (snp) / 1.202 (nonsnp) / 1.065 (sv).

**Result — peaks substantially concordant, not identical:** the single strongest JOINT hit
genome-wide sits in the **exact same 20kb window (chr2 ~2.34Mb, window idx 117) for both classes
independently**. Genome-wide, per-20kb-window max(-log10 p_joint) correlates Spearman ρ=0.549
(p≈0, n=5,483 shared windows) between SNP and non-SNP. Top-hit window overlap (Jaccard): 0.25 at
the top 0.1%, 0.32 at top 0.5%, 0.20 at top 1% (strongest peaks coincide; the tail diverges).
SV-vs-SNP is far weaker (ρ=0.126, top-0.1% Jaccard 0.14) — SV's sparse 12.8k-marker scan doesn't
converge on the SNP peaks the way non-SNP does. Both classes: **0 significant CLIMATE hits**
(consistent with the established climate-null), but real JOINT hits — 751 FDR<0.05 / 49 Bonferroni
(SNP) vs 160 / 24 (non-SNP) vs 10 / 2 (SV), roughly tracking the marker-count ratios. Outputs:
`analysis/grenenet_gea/varexp/class_gwas_{snp,nonsnp,sv}.npz` + `class_gwas_summary.json`.

**Net:** sharpens the kinship-level conclusion — it's not only that SNP/non-SNP panels explain
the same *total* variance, an actual per-marker scan on each panel independently converges on
the same top locus and same broad regions (with real but imperfect concordance further down the
hit list). Running this GWAS with SNPs alone tells essentially the same story as SVs/indels alone.

## Design A (2026-07-06): Zhou-2022 (Nature) analog — SV h²/prediction of founder ORIGIN CLIMATE

Motivation: Zhou et al. 2022 (`papers/Zhou et al. 2022 - Nature.pdf`, repo `YaoZhou89/TGG`) show a
graph pangenome (SNP+indel+SV) captures more trait heritability (0.41 vs 0.33 linear, +24%) and
higher genomic-prediction accuracy than a single linear reference, driven by SVs resolving
incomplete LD / allelic heterogeneity. User asked whether kMate's short-read→long-read (78 cactus
long-read + 153 PanGenie) panel can show the same.

**Key design decision (why this is NOT the null above):** in a POOL experiment there is no
independently-measured per-ecotype field fitness — the only per-ecotype "fitness" is `s` (the pooled
logit-slope), which is exactly the trait the variance-partition above already found null. To make a
genuinely Zhou-faithful individual-level test, the phenotype must come from OUTSIDE the pool.
Chosen phenotype = per-founder **origin climate bio1–19** (worldclim), keyed to founders by
`ecotypeid`. Genome predicts provenance via local adaptation; independent of the pool.

- Phenotype source: `gea_grene-net/key_files/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv`
  (1871 ecotypes × bio1-19; **all 231 GRM founders matched, complete bioclim**).
- Genotype "two ways" = short-read→long-read contrast INSIDE the arch3 panel, reusing the prebuilt
  GRMs: SNP arm `K_snp` vs pangenome arm `K_all`/`+K_sv` (the indel+SV layers need the long-read
  graph). Marginals K_indel/K_sv/K_nonsnp, control K_snp_matched, mechanism `K_untagged_r02`.
- Methods: reuse `varexp_selection.py` estimators verbatim (AI-REML marginal h²+SE+LRT, joint 2-GRM
  partition + LRT for the SV/non-SNP layer beyond SNP, repeated 6-fold×10 GBLUP CV predictive R²).
  Headline gain = CV R²(SNP+SV) − R²(SNP), across bio1-19 × {raw, RINT}.
- Script `varexp_bioclim.py` → `analysis/grenenet_gea/varexp/bioclim_varexp.csv` +
  `bioclim_varexp_summary.json`. Env kmate. GRM corrs at the matched set: **corr(K_snp,K_sv)=0.955,
  corr(K_snp,K_nonsnp)=0.998, corr(K_snp,K_untagged_r02)=0.364.**
- **Prior expectation (honest):** aggregate GRM gain likely small — origin climate is heavily
  predicted by genome-wide relatedness (structure = geography), and K_snp≈K_all, same as the `s`
  result. The place a positive could still appear is single-SV LASSO (allelic heterogeneity at
  big-effect climate loci) — deferred to a follow-up (needs raw marker genotypes vs bioclim, not
  just GRMs). RESULT: see below once `logs/varexp_bioclim.out` completes.

**Step 1 (GRM/GBLUP) status:** `varexp_bioclim.py` running (bg). Reports, per bio1-19 × {raw,RINT}:
marginal h² (snp/nonsnp/indel/sv/all/snp_matched/untagged_r02), joint 2-GRM LRT for SV- and
non-SNP-beyond-SNP AND untagged-beyond-SNP, and GBLUP CV gains **`gain_sv` = R²(snp+sv)−R²(snp)**
and **`gain_nonsnp` = R²(snp+nonsnp)−R²(snp)** (non-SNP layer as a whole, per user request — indels
are the bulk of it). Early (bio1): R²(snp)=+0.56, gain_sv=−0.005, gain_nonsnp=−0.003 → tracking the
predicted null (origin climate ~56% genome-predictable from SNPs, SV/non-SNP layer adds ~0).

**Step 2 (multilocus LASSO) — BUILT & QUEUED (hold until step 1 confirms null):**
`lasso_bioclim.py` + `lasso_bioclim.sbatch` (savio4_htc/co_moilab, 32c/96G/12h). Python/sklearn
ElasticNetCV (l1_ratio=1 = LASSO; glmnet not in r_env). Per bio-trait, arms = LD-pruned SNP
backbone (r²<0.9, ±50kb, cap 150k) vs +SV vs +all-nonSNP (non-SNP kept UNPRUNED so an untagged SV
can't be pruned away); SNP backbone identical across arms so any gain = the non-SNP layer, not
marker count. Reports CV R² gain + every non-SNP marker with nonzero coef, its SNP-tag r² (from
`nonsnp_tagging_chr*.npz`) and gene — an untagged (r²<0.2) selected marker = the incomplete-LD
smoking gun (Zhou Fig.3 mechanism). Reuses `build_class_grms._load_chrom/_classes` for the exact
marker filter. LD-prune + tagging loader unit-tested OK. → `lasso_bioclim{.csv,_selected.csv,_summary.json}`.

## Downstream — overlap tables, Manhattans, candidate genes (merged from 2026-07-03 session log)

**Overlap tables** (`notebooks/class_gwas_multitrait.ipynb`, `class_gwas_persite.ipynb`, both
rebuilt & executed 2026-07-21) mirror the phase-1 kendall/lfmm/binomial `overlap_df` at the
**clq0.9 block level** (one min-p lead per LD block), with both Bonferroni and FDR counts. Current
JOINT block counts: **Bonferroni 14 (snp) / 10 (nonsnp) / 2 (sv)**, of which 5 non-SNP blocks are
not SNP-significant; **BH-FDR 8,269 (snp) / 1,000 (nonsnp) / 17 (sv)**, of which 258 non-SNP blocks
are not SNP-significant. GLOBAL + CLIMATE_bio1: ~0. Temperature bioclim vars show some FDR-only SNP
blocks that mostly vanish under Bonferroni, don't replicate in non-SNP, and aren't corrected across
the ~20 collinear climate axes → **FDR-tail noise, not real climate adaptation** (the established
climate-null). Bug fixed en route (still in place): `lib.collapse_to_blocks` picks the block lead by
MAX |stat|; passing `-p` made `abs(-p)` pick the *least* significant marker — switched to `-log10(p)`.

**Plots:** `plot_class_gwas_pngs.py` → PNGs in `analysis/grenenet_gea/varexp/gwas_plots/`
(multitrait = contrasts × 3 classes × {manhattan,qq}; per-site = gardens × 3 classes ×
{manhattan,qq}). Threshold lines drawn at the **block level** (fix: originally marker-level, so the
FDR line silently vanished when marker-level FDR was empty); legend shows significant-block count.
These PNGs are built from `class_gwas_{snp,nonsnp,sv}.npz`, which reproduced **bit-identically** on
the 07-21 rerun, so the existing plots remain valid; **not re-rendered this pass** (re-run the
script if the 30-vs-31-garden enumeration in filenames matters).

**Non-SNP-only candidate genes** (`nonsnp_only_genes.py` → `nonsnp_only_{blocks,genes}.csv`;
`nonsnp_only_genes_describe.py` → `..._described.csv` via Ensembl Plants + UniProt REST) —
**regenerated 2026-07-21 on the current `class_gwas_*.npz`.** Across all 52 contrasts (multitrait
JOINT/GLOBAL + CLIMATE×[bio1-19,PC1] + 30 per-site gardens), clq0.9 blocks non-SNP-FDR-significant
but NOT SNP-FDR-significant → **431 blocks (25 at Bonferroni) → 969 genes** (465 in-block, 504
flank-only, ±2 kb promoter). Top surfaced symbols on the current data include **FTSH1, MYB112,
MYB89, NAC093, SUMO5, and PAP-family** genes.
> ⚠️ **The old themed narrative (ADS2 cold-desaturase, GIGANTEA, HSP70, ERA1) is SUPERSEDED** — it
> came from the pre-fix (2026-07-03) run and does **not** reappear as the headline set on current
> data. The regenerated gene list above is comprehensive but un-interpreted; the *biological theme*
> assignment (which hits are cold/heat/ABA-relevant, which are sturdy vs fragile-tail) is for the
> PI to re-review before anything is cited.

**Caveat carried throughout:** this is the *fragile tail* by construction (blocks one class calls and
the other doesn't) — mostly FDR-only, per-site, and/or flank hits, inheriting the bioclim-null caveat.
A hypothesis-generating list, not confirmed loci.

## Deliverables still TODO
- **PI re-review of the regenerated candidate-gene themes** (`nonsnp_only_genes*` rerun 07-21 on
  current data → 431 blocks / 969 genes; the pre-fix ADS2/GI/HSP70/ERA1 list is superseded; current
  top symbols FTSH1/MYB112/MYB89/NAC093/SUMO5/PAP-family) — assign biological themes / sturdy-vs-tail.
- Cross-axis multiple-testing correction on the bioclim CLIMATE hits (20 independent BH runs on
  collinear temperature vars) to confirm they are noise.
- Pull the actual indel/SV variants inside whatever candidate blocks survive the regeneration above
  (size, freq, position vs gene); annotated locus plot for the sturdiest in-block JOINT hit.
- Annotate/inspect the shared chr2~2.34Mb top locus (what gene/region, near a known candidate?).
- Figure notebook (`basic` env — matplotlib HANGS in `plotting`; compute→npz, render elsewhere):
  per-class VE across the 30 sites + SNP→SNP+nonSNP gain panel, plus the untagged-GRM gain panel,
  plus a SNP-vs-nonSNP Manhattan overlay.
- class_split_gwas.py now emits the **3-way SNP/indel-nonSNP/SV split** (sv arm added); a folded
  3-way (SNP/indel/SV) breakdown *table* in the notebook is still TODO.
- Census-flavour trait as robustness (currently relative only).
- Optional robustness: per-site (not just aggregate-axis) untagged-GRM gain scan, mirroring
  `varexp_selection.py`'s per-site loop — not yet done for the untagged GRMs.
- Annotate/inspect the shared chr2~2.34Mb top locus (what gene/region, is it near a known
  candidate) — not yet done.

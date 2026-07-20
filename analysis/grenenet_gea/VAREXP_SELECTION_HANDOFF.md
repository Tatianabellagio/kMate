# SNP vs non-SNP variance partition of ecotype selection — result + handoff

> ## ⚠️ NUMBERS PREDATE THE Kf_w / `--unit chrom` FIX (2026-07-08) — need regeneration
> Every quantity in this doc (h²=0.938, all joint LRT p-values, corr(K_snp,K_nonsnp)=0.998,
> per-site λ=1.040, the chr2 ~2.34Mb top peak, all GWAS hit counts) was computed on **PRE-FIX
> founder h** — collapsed-founder global-EM, before the full-panel Kf_w correction and the move to
> `--unit chrom`. They **must be regenerated** and should not be cited as current. The cited npz
> inputs (`ecotype_fitness/sample_global_h.npz`, `varexp/selection_s_matrix.npz`) have since been
> regenerated under the fix. In particular the founder-collapse (below) is fixed: founder 9977 is now
> rescued to h≈1.5e-4 (was ~1e-15, "absorbed") and ~230 of 231 founders are retained.

Session date: 2026-07-03. Branch `add-kmate`. Env: `kmate`
(`source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate`). Always `hostname`
first — compute only on `n*.savio*`, never `ln00X` (hook-enforced).

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
  (`S[31×231]`, `sites`, `founders`, `bio1`, `p0`, `analyzable`(bool), `n_plots`, `freq_last`).
- GLOBAL-mode chrom-averaged `h` (reused from `ecotype_fitness/sample_global_h.npz`, GLOBAL_MODE_DECISION).
- founding reference `p0` = mean over 8 SEEDMIX reps (ESTIMATED, not forced uniform 1/231 — the
  twin-absorption identifiability bias cancels in the slope; verified seedmix averaging is correct,
  mean is exactly 1/231 by closure, spread is real & reproducible across all 8 reps).
- **presence filter**: analyzable = `p0 > 1e-3`. **[PRE-FIX: kept 212 of 231, dropping ~19
  twin-absorbed founders like 9977 at h≈1e-15.]** After the full-panel Kf_w fix the founder collapse
  is gone: **9977 is rescued to h≈1.5e-4 and ~230 of 231 founders are retained**, so the `1e-3` floor
  is now a near-no-op rather than a mass exclusion. This filter — and every downstream count — needs
  to be recomputed on the post-fix h.
- **NO reliability/cross-chrom-SD weighting** (user decision: chrom-averaging already regularizes).
- raw `s` primary; **RINT** (rank-inverse-normal per site) as sensitivity.

## What was built

1. `build_selection_trait.py` → `selection_s_matrix.npz`. **DONE.** (skew −0.98, kurt 5.7, 31 sites, 212 founders.)
2. `build_class_grms.py` + `run_class_grms.sbatch` → `analysis/grenenet_gea/varexp/class_grms.npz`.
   **DONE** (job 35513509). Per-class founder GRMs (231×231), standardized `Z=(g−p)/√(p(1−p))`
   — **0/1 founder presence = haploid coding** (inbred lines; NOT the diploid 2p form). Markers:
   MAC≥5, call≥0.9. `K_snp` (1.75M), `K_indel` (512k), `K_sv` (12.8k), `K_nonsnp` (525k),
   `K_all`, `K_snp_matched` (SNPs subsampled to non-SNP count & MAC spectrum, 525k). Aligns to
   trait founders by ID in the modeling step (var_pa founder ORDER differs from seedmix order).
3. `varexp_selection.py` → `varexp.csv`, `varexp_persite.npz`, `varexp_meta.json`.
   Estimators: (1) marginal single-GRM AI-REML h²+SE+LRT, (2) joint 2-GRM REML partition + LRT for
   non-SNP-beyond-SNP, (3) GBLUP repeated k-fold CV predictive R² & gain. Traits = 4 aggregate axes
   (w_global + bio1-tercile zones) × {raw,RINT} + 31 per-site (raw). **RUNNING in background**
   (pid 2192846, `logs/varexp_selection.out`, NFOLD=6 NREP=10; ~15-20 min). Check it finished:
   `ls analysis/grenenet_gea/varexp/varexp.csv`.

## KEY FINDING (genome-wide GRM) — non-SNP is redundant with SNP

On `w_global`: marginal h² = **0.938 (snp) vs 0.937 (nonsnp) vs 0.938 (matched)** — identical;
joint 2-GRM **LRT p_add_nonsnp = 0.32** (n.s.); CV **gain = −0.002** (R² snp 0.491, nonsnp 0.486,
both 0.489). Driver: **corr(K_snp, K_nonsnp) = 0.998–0.999.** The selection trait is ~94%
kinship-predictable (winners are whole clades), but SNP and non-SNP GRMs are the *same genealogy*,
so non-SNP adds nothing at the polygenic level. Consistent with SVs=passengers + SNPs-tag-86%-of-SVs.

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
   corr with K_snp = **0.396** — genuinely decorrelated) and `K_untagged_r05` (42,140 markers,
   corr with K_snp = 0.953 — still fairly tied to it). ->
   `analysis/grenenet_gea/varexp/untagged_grms.npz`.
3. `varexp_untagged.py` — reruns marginal h², joint 2-GRM REML+LRT, GBLUP CV gain (imports the
   estimators straight from `varexp_selection.py`) for `K_snp` + `K_untagged_{r02,r05}`, across
   the 4 aggregate axes (w_global + 3 bio1-tercile zones) × {raw, RINT}. ->
   `analysis/grenenet_gea/varexp/varexp_untagged.csv`. (Minor bug, harmless: the summary printer
   crashes on `df.thresh` — pandas resolves `.thresh` to a builtin, not the column; fixed to use
   `df["thresh"]`, but the CSV itself is already complete before the crash — no rerun needed.)

**Result (w_global, raw, r²<0.2 — the strict/primary threshold):** marginal h²(untagged) =
**0.889** (SE 0.112, still huge!) — because even locally SNP-invisible markers still reconstruct
the SAME founder clade structure (kinship is a genome-wide, not local, property; unlinked
markers anywhere still tag the genealogy). But conditional on `K_snp`: joint LRT **p_add = 0.35**
(n.s.), and CV **gain = +0.0004** (r² snp 0.484 vs both 0.485) — trivial and not a real effect
size. **Fully robust**: across all 4 axes × 2 transforms × 2 thresholds (16 cells), CV gain never
exceeds |0.004| and flips sign axis-to-axis (noise, not signal); every joint-LRT p_add is ≥0.24
except one w_cold/RINT/r02 cell at p=0.051 (not significant after 16 looks, and its raw-scale
twin is p=0.24 — not corroborated).

**Interpretation:** the "no kMate gain" conclusion is now a genuinely strong, non-tautological
result. It's not that K_nonsnp merely LOOKS like K_snp (the genome-wide confound) — even markers
SNPs cannot locally tag still only recover the same founder genealogy, and genealogy is already
~89–94% predictive of the selection trait via SNPs alone. There is no local, SNP-independent
non-SNP signal riding on top of clade structure. This is the strongest form of "SVs/indels are
passengers, not an independent selection currency" reached in this whole investigation thread —
consistent with [[project-gea-sv-selection-currency]]'s hap-cluster-level null and the parallel
SuSiE fine-mapping session's independent passenger conclusion.

## Peaks + per-marker GWAS comparison (2026-07-03, same session) — SNP-only vs non-SNP-only scan

Distinct from the kinship/variance-partition question above: does an actual **single-marker**
GWAS (not a GRM) run separately on SNP-only vs non-SNP-only markers find the same **peaks**, at
both per-site and multi-site (JOINT) resolution? New script `class_split_gwas.py`:

- Same locked trait (`s`, 212 analyzable founders), rank-normalized per site.
- **ONE shared LOCO kinship correction per chromosome** (built from ALL classes pooled,
  MAC≥12/call≥90%, total-minus-chromosome trick — no repeated whole-genome loads) — deliberately
  shared across both scans so any difference in results is attributable only to the TEST-marker
  class, not to also varying the correction (justified: K_snp/K_nonsnp are near-identical per
  the variance-partition result above).
- Test markers MAC≥5/call≥90% (matches `class_grms.py`), split SNP (1,750,565 markers) vs
  non-SNP (536,432, indel+SV pooled). Per (chrom, class): one eigendecomposition + one
  U.T-rotation of the marker block, reused across all 31 sites (only REML-delta + GLS re-runs
  per site) — full genome x 31 sites in ~5.5 min interactively, no sbatch needed.
- Per-site z genomic-control calibrated, stacked into per-class Z[M,31], run through the same
  Bolormaa JOINT(31df)/GLOBAL(1df)/CLIMATE(1df) meta as `founder_gwas_multisite.py`.
- **Numerical gotcha hit and fixed:** ~7–10 markers per chromosome (out of ~450k) are
  near-FIXED (202–207/212 founders carry the allele, not rare — MAC treats both tails
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

**Result — variance explained is identical two independent ways:**
mean per-site genomic-control λ = **1.040 (SNP) vs 1.040 (non-SNP)** — literally the same
inflation, from real single-marker testing (not REML), corroborating the earlier GRM h² result
(median per-site h²: 0.951 snp vs 0.952 nonsnp; median CV R²: 0.549 vs 0.545, from
`varexp_persite.npz`).

**Result — peaks substantially concordant, not identical:** the single strongest JOINT hit
genome-wide sits in the **exact same 20kb window (chr2 ~2.34Mb) for both classes independently**.
Genome-wide, per-20kb-window max(-log10 p_joint) correlates Spearman ρ=0.545 (p≈0, n=5,477
shared windows) between classes. Top-hit window overlap: 67% Jaccard at the top 0.1%, decaying
to ~29% at top 0.5–1% (strongest peaks coincide; the longer tail diverges more). Both classes:
0 significant CLIMATE hits (consistent with the established climate-null), but real JOINT
hits — 497 FDR<0.05 (SNP) vs 149 (non-SNP), roughly tracking the ~3.3x marker-count ratio, not
obviously different per-marker power. Outputs:
`analysis/grenenet_gea/varexp/class_gwas_{snp,nonsnp}.npz` + `class_gwas_summary.json`.

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

**Overlap tables** (`notebooks/class_gwas_multitrait.ipynb`, `class_gwas_persite.ipynb`) mirror the
phase-1 kendall/lfmm/binomial `overlap_df` at the **clq0.9 block level**, with both Bonferroni and
FDR counts. JOINT: 21 Bonf / 1364 FDR (snp) vs 9 / 104 (nonsnp), 82 shared, 79% of non-SNP hits are
also SNP hits. GLOBAL + CLIMATE_bio1: 0/0. Several temperature bioclim vars show FDR-only SNP hits
(bio10=140, bio11=97, bio13=72, bio16=64) but they mostly vanish under Bonferroni, don't replicate
in non-SNP, and are not corrected across the 20 collinear climate axes → **FDR-tail noise, not real
climate adaptation.** Bug fixed en route: `lib.collapse_to_blocks` picks the block lead by MAX |stat|;
passing `-p` made `abs(-p)` pick the *least* significant marker — switched to `-log10(p)`.

**Plots:** `plot_class_gwas_pngs.py` → **318 PNGs** in `analysis/grenenet_gea/varexp/gwas_plots/`
(132 multitrait = 22 contrasts × 3 classes × {manhattan,qq}; 186 per-site = 31 gardens × 3 classes ×
{manhattan,qq}). Threshold lines drawn at the **block level** (fix: originally marker-level, so the
FDR line silently vanished when marker-level FDR was empty); legend shows significant-block count.

**Non-SNP-only candidate genes** (`nonsnp_only_genes.py` → `nonsnp_only_{blocks,genes}.csv`;
`nonsnp_only_genes_describe.py` → `..._described.csv` via Ensembl Plants + UniProt REST): for each of
53 contrasts, clq0.9 blocks non-SNP-FDR-significant but NOT SNP-FDR-significant → **53 blocks (8 at
Bonferroni) → 126 genes** (71 in-block + 55 flank-only, ±2 kb promoter). Themed hits (all
hypothesis-level):
- **ADS2 (AT2G31360)** — Δ9 fatty-acid desaturase (membrane cold-acclimation). **Sturdiest: in-block,
  multitrait JOINT** (the one well-calibrated contrast).
- **GIGANTEA (GI, AT1G22770)** — clock/photoperiod/freezing regulator. Eye-catching but statistically
  weak: garden-4 only, ±2 kb *flank* of block Chr1_4247, one marker just over FDR, well under Bonferroni.
- Cold: SEX1/GWD, PI-4KBETA2. Heat: HSP70 (AT4G16660), HIP1, a Clp-N chaperone.
- ABA/drought: **ERA1 (AT5G40280)** farnesyltransferase β (only Bonferroni-level themed hit), XERICO, NAC032.

**Caveat carried throughout:** this is the *fragile tail* by construction (blocks one class calls and
the other doesn't) — mostly FDR-only, per-site, and/or flank hits, inheriting the bioclim-null caveat.
A hypothesis-generating list, not confirmed loci.

## Deliverables still TODO
- Cross-axis multiple-testing correction on the bioclim CLIMATE hits (20 independent BH runs on
  collinear temperature vars) to confirm they are noise.
- Pull the actual indel/SV variants inside the ADS2 and GI blocks (size, freq, position vs gene);
  annotated locus plot for ADS2 (as done for GI in `gwas_plots/GI_locus_site4_nonsnp.png`).
- Annotate/inspect the shared chr2~2.34Mb top locus (what gene/region, near a known candidate?).
- Figure notebook (`plotting` env — matplotlib HANGS in `plotting`; compute→npz, render elsewhere,
  or use `basic`): per-class VE across 31 sites + SNP→SNP+nonSNP gain panel, plus the new
  untagged-GRM gain panel, plus a SNP-vs-nonSNP Manhattan overlay.
- 3-way (SNP/indel/SV) breakdown table (K_indel, K_sv already in class_grms.npz; class_split_gwas.py
  currently only does the 2-way SNP/nonSNP split).
- Census-flavour trait as robustness (currently relative only).
- Optional robustness: per-site (not just aggregate-axis) untagged-GRM gain scan, mirroring
  `varexp_selection.py`'s per-site loop — not yet done for the untagged GRMs.
- Annotate/inspect the shared chr2~2.34Mb top locus (what gene/region, is it near a known
  candidate) — not yet done.

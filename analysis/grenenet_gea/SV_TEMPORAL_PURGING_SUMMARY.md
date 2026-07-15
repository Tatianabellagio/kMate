# SV temporal purging — thread summary (2026-07-02)

> ## ⚠️ CORRECTED (2026-07-15, second pass) — split verdict, not a uniform null
> **[2026-07-08 pass, still correct]** The de-trended metric (`s_distribution_by_site`, the per-p0-bin
> ALL-class-median subtraction) shows the **OVERALL / whole-distribution median shift is essentially
> ZERO: median +0.0011, negative at only 15/31 sites, n.s.** No genome-wide, frequency-independent SV
> purifying excess on this specific (central-tendency) statistic. The old headline number ("~22/31
> sites, three methods agree") for THIS statistic no longer holds.
>
> **[2026-07-15 audit, corrects the 2026-07-08 banner]** This banner previously ALSO dismissed
> `s_histogram`, `s_vs_climate`, `ecdf`, `ecdf-difference`, and `shiftfunction` as sharing the same
> confound — **that was checked directly and is wrong.** All of these (plus climate-slope β) describe
> the same thing in their own takeaways: not a whole-distribution shift, but a **heavier purged TAIL,
> concentrated at hot sites**. That is a different statistic than the de-trended median-shift metric,
> and it was directly tested: a 10th-percentile SV-vs-matched-SNP tail gap, de-trended the same way as
> `s_distribution`, barely moves (hot-site median −0.071→−0.055; correlation with bio1 strengthens,
> ρ=−0.47→−0.51, p=0.0075→0.0032). **The tail/hot-site signal is real and is not the artifact that
> explains the median-shift null.** Do not archive or dismiss `s_histogram`/`s_vs_climate`/`ecdf`/
> `ecdf-difference`/`shiftfunction` as redundant with `s_distribution_by_site` — they test a different,
> still-standing question. See the "Audited 2026-07-15" sections below for the climate-slope β detail
> and the direct tail de-trending test.

**Question:** In GrENE-net, comparing per-variant allele-frequency change across generations, do
non-SNPs (indels + SVs) show more selection than SNPs? (Pure temporal, per-variant — no haploblocks.)

## Estimator

Per variant per site, **s = plot-replicate selection coefficient** = mean over the site's ~10–12
replicate plots of the OLS slope of `logit(p)` on generation (gen 0 = shared founding p0).
`s < 0` = declining (purged); `s > 0` = rising (favoured). Using plots as replicates separates
selection from drift (among-plot SE). Code: `_temporal_s_plots_snp_vs_nonsnp.py::plot_replicate_sz`.

**Audited** (`_audit_s_classes.py`): class masks correct & disjoint (SNP `dlen=0`, indel 1–50 bp,
SV >50 bp); saved per-class arrays are genuinely different data; an independent from-scratch recompute
of `s` reproduced the function **exactly (Δ=0)**. The estimator is correct.

All class comparisons use a **frequency-matched SNP null**: SNP/indel resampled to the SV `p0`
distribution, so "SVs are rarer" is not a confound.

## Result

- **indels ≡ SNPs** at every frequency and in every view → the non-SNP *category* (98% indels) is also
  ≈ null. (`_nonsnp_temporal_category.py`, `_temporal_selection_snp_vs_nonsnp.py`: fold ≈ 1.0.)
- **[SUPERSEDED]** The raw/single-tail views below suggested "SVs carry a small excess of purifying
  (negative-s) selection in the purged tail that grows with climate harshness." **This does not
  survive the Kf_w / `--unit chrom` regeneration.** On the frequency-de-trended
  `s_distribution_by_site` metric the SV excess-vs-baseline is ≈0 (median **+0.0011**, negative at
  only **15/31 sites**, sign test n.s.). The apparent tail effect was a frequency /
  founder-projection confound, not a genome-wide purifying excess.

### [AUDITED 2026-07-15 — dismissal was WRONG, based on a mislabeled statistic] Cleanest test — per-variant climate slope β = d s / d climate  (`_compute_s_climate_slope.py`)
> Was marked "superseded," citing "per-axis correlations weaker/mixed than stated below: bio18 corr
> = +0.54, bio1 corr = −0.16." **That citation is a mislabeled statistic, not the SV-specific
> finding.** `+0.54/−0.16` is `s_purging_intensity.png`'s `corr(intensity, ...)` — the GENERAL
> matched-SNP purging-rate-vs-climate correlation (a confound check, computed on `msnp_neg` alone) —
> not the SV-specific **sign-excess** (`sv_neg − msnp_neg`) vs climate correlation the historical
> table below actually reports. Recomputed the real sign-excess directly from the current (post-Kf_w)
> `s_climate_slope_sign_by_site.csv`: **bio1 ρ=+0.51 (p=0.0035), bio18 ρ=−0.55 (p=0.0013)** — nearly
> unchanged from the historical bio1 ρ=+0.55 / bio18 ρ=−0.65, both still highly significant, same
> sign. **This arm's SV-specific climate-graded signal is essentially intact post-Kf_w and was never
> actually shown to be superseded.**
>
> **Direct de-trending test run 2026-07-15** (now Section 4 of `temporal_s_consolidated.ipynb`,
> `_build_temporal_s_consolidated_nb.py` — the standalone script this was first run in has been
> folded into the notebook and removed): applied the
> exact same per-p0-bin ALL-class-median subtraction that killed `s_distribution_by_site`'s median
> shift, then redid the sign-excess-vs-climate correlation on the residualized `s`. **The signal
> survives — it does not weaken:** bio1 ρ +0.38→+0.42 (p=0.019→0.019), bio18 ρ −0.54→−0.59
> (p=0.0015→0.0004); median excess flips from −0.0037 to **+0.0016** with more sites positive
> (10/31→19/31). This is the opposite of what you'd see if the p0/logit-boundary artifact explained
> the signal. **Conclusion: this arm's climate-graded SV signal is not the same artifact that killed
> the median-shift metric, and is not currently explained by any confound found in this thread.** The
> standing hitchhiking / global-mode caveat (see below) still applies — this doesn't prove direct SV
> selection vs. linkage to a selected haplotype — but the "superseded" label is retired.

A variant on a *constantly* purged haplotype (hitchhiking) has s<0 everywhere → β≈0; the p0/logit
artifact is a per-variant constant → cancels in the slope. So β isolates climate-**differential**
selection with the two big confounds removed.

| climate axis | SV vs matched-SNP | ins vs del | sign-excess vs climate |
|---|---|---|---|
| **bio1** (temperature) | KS p = 5×10⁻¹⁶ | **insertion-only** (ins β=−0.010, p=2×10⁻³⁹; del null) | ρ = +0.55 |
| **bio18** (dry-summer precip = aridity) | KS p = 4×10⁻¹² | both (ins p=7×10⁻¹⁶, del p=6×10⁻⁴) | **ρ = −0.65** |

- **Mechanism = insertion-vs-reference polarity, not SV length.** For temperature, only insertions are
  purged; deletions null. Small indels at SV-like frequency also show it → the real axis is
  *insertions of any size relative to Col-0/TAIR10 reference*, climate-graded.
- **Climate = hot + arid (Mediterranean/harsh-summer);** aridity (bio18) is the stronger correlate.

### Confound ruled out — NOT "harsh sites purge everything more"
1. The excess is a **within-site** SV − matched-SNP difference → site-wide intensity is differenced out.
2. Overall purging intensity (fraction of all matched SNPs purged, 0.55→0.80) is **not** higher at
   harsh sites: it correlates **+0.60 with bio18 (higher at WET summers)**, flat with temp
   (bio1 ρ=−0.24, ns). The SV excess trends the **opposite** way (arid) and **anti-correlates** with
   overall intensity (ρ=−0.40) → appears where general purging is weakest.
   *(Correction: an earlier claim that overall purging rises with aridity was backwards.)*
3. General purging cannot produce the insertion/deletion asymmetry.

## Method lessons (reusable)
- **logit-slope Jensen/boundary artifact:** under symmetric drift, rare variants get a negative median
  logit-slope, common a positive one (zero-selection Wright-Fisher sim reproduces the exact NEG→POS
  vs-p0 trend; median *linear* Δp stays ~0). It **cancels** in the frequency-matched class comparison —
  the p0→s trend is not "selection stronger on common variants."
- **Don't over-filter MAF:** the SV signal lives in rare variants; MAF≥0.10 erased it. Frequency-match
  instead.
- **Overlaid histograms can't resolve a small shift; finer bins make it worse.** Use ECDF-difference
  (x = s, y = ΔCDF) or selection-difference-vs-independent-axis (climate), not selection on both axes.

## Replicate-based arm (session 2, 2026-07-03) — PARALLELISM + PicMin
Uses the ~10–12 replicate **plots** within each site as parallel populations (each an independent
pool-seq → independent allele-frequency-change). Consistency across plots = drift control.
- **Parallelism** ρ = mean²/mean(slope²) across plots (AF-vapeR rank-1 eigenvalue analog; 0=drift,
  1=fully parallel). `_compute_parallelism.py`. **Responder** = ρ in top decile within its p0-bin
  (class-agnostic). **Repeatability** = # sites where responder.
- **SV enrichment among parallel responders** (freq-matched, bootstrap p=0.005): responder rate
  **1.16×**, repeatable ≥⅓ sites **1.72×**, strongly-repeatable ≥½ sites **3.65×** (escalates with
  stringency). **Insertion-driven** (insertions 0.098 vs deletions 0.038). **Holds at ALL MAF** (rare
  1.76×, mid 1.28×, common 1.29× — so no MAF filter needed; unlike the climate-purging median it
  survives common MAF).
- **Real PicMin** (`_picmin.py`: empirical per-site p vs SNP in p0-bin → Beta order statistics over
  the 31 site-lineages → min-over-orders → uniform-null calibration → BH-FDR): SV **1.23×** (q<0.1) /
  **1.33×** (q<0.05) enriched among repeated-adaptation loci; **insertion-driven** (ins 22–25% vs del
  ~13% ≈ indel ≈ SNP). Absolute significant frac is high (13–20%) because this founder-projection
  system has pervasive parallel sorting — the SV-vs-SNP **relative** contrast is the signal.
- **Per-site climate grid** (`parallelism_by_site.ipynb`, analog of s_ecdf_difference_by_site): SV
  parallelism excess is **climate-graded** — bio1 ρ=+0.43 (p=0.016), **bio18 ρ=−0.70 (p=0.000)**;
  the matched-SNP baseline parallelism is ~flat (+0.10 / −0.19) → the SV-specific gap opens at hot/
  arid sites (aridity-dominant), same as the climate-slope β arm.
- **AF-vapeR NOT run** (window eigen-method; SVs ~0–1 per window → can't give a per-SV parallelism;
  would reduce to SV-window colocalization = spatial-null problem). PicMin is the right replicate tool.
- **I/O GOTCHA (important for any per-plot analysis):** memmap random-row access of `pool_gen*_af.npy`
  is ~1 MB/s (page-fault thrash) on this filesystem; use direct **seek + np.fromfile** (`read_rows`
  in `_compute_parallelism.py`) = ~150 MB/s. Sequential `dd` is 560 MB/s.

**[AUDITED 2026-07-15 — mixed result, not a uniform "superseded"]** This section previously claimed
the "three methods agree" headline does not survive the Kf_w / `--unit chrom` rerun because all three
arms "share the frequency / founder-projection confound" that de-trended `s_distribution_by_site`
removes. That was asserted by analogy and never tested directly on these statistics. Checked against
the actual pre- vs post-Kf_w rerun numbers (git history, commit `3c34ec8`):
- **Parallelism (bulk responder rate)**: shrank ~40% (insertions 0.098→0.059) — consistent with, but
  not proof of, the same confound. Not re-verified null.
- **PicMin**: SV-vs-matched-SNP fold shrank from 1.27× to 1.07× — closer to parity, but not re-verified
  null with a dedicated significance test on the post-Kf_w numbers.
- **Climate-slope β (sign-excess vs climate, the actual SV-specific number)**: essentially
  **unchanged** (bio1 ρ +0.55→+0.51 p=0.0035, bio18 ρ −0.65→−0.55 p=0.0013 — recomputed directly
  from current data). The doc previously cited "+0.54/−0.16, weaker/mixed" here — that was a
  **mislabeled, unrelated statistic** (the general purging-intensity-vs-climate confound check, not
  the SV-specific sign-excess) and should not have been used as evidence of anything.
- **Parallelism-by-site climate-gradient excess**: numbers **essentially unchanged** (bio1
  +0.43→+0.44, bio18 −0.70→−0.68).

The last two being untouched by the exact fix that killed `s_distribution`'s signal means their
dismissal is **not supported by the evidence actually available** — and for climate-slope β the
original dismissal was actively wrong (mislabeled statistic), not just unverified. **Update
2026-07-15: ran the direct de-trending test on climate-slope β's sign-excess (see the climate-slope
β section above) — the signal SURVIVES de-trending (if anything slightly strengthens).** So:
`s_distribution_by_site`'s own de-trended metric is directly verified null; parallelism/PicMin
shrank but aren't confirmed null either way; **climate-slope β is now a verified, not-yet-explained
real signal**; the parallelism-by-site climate-gradient remains open (same de-trending test not yet
run on it). Do not cite this section as "all three methods died." The historical multi-method text is
retained above for the record.

## Caveats (unresolved)
Global-mode kMate AF is a linear projection of per-sample founder h → cannot separate SV-specific
selection from **hitchhiking** on climate-purged haplotypes; and insertion-polarity could be a founder
insertion-**calling / reference-mapping artifact**. Both need an **independent (local-mode / vg) SV
allele frequency** + ancestral polarization to resolve.

**Note on local-mode kMate as "the fix" (2026-07-03):** `GLOBAL_MODE_DECISION.md` (written the prior
session, same investigation) argues against using local/window-mode `h` for anything selection-adjacent
in this system — panel-incompleteness artifacts are spatially-coherent, deterministic, and
**replicate-reproducible**, which would masquerade as exactly this SV-specific/climate-graded signal,
and the natural "does it reproduce across replicates" filter runs backwards (real drift-driven signal
gets filtered out, artifacts get kept in). So local-mode is not a clean fix for the hitchhiking caveat
without first working through that confound (e.g. gen0 empirical null, leave-one-founder-out).

**Calling-quality proxy for the insertion-polarity artifact (2026-07-03, partial answer).** No outgroup
genome/alignment exists anywhere in this repo (checked: no *A. lyrata*/*arenosa*/*rubella* files, no
AA/ancestral annotation, no non-founder whole-genome alignment) — true ancestral polarization would be a
multi-step pipeline build from scratch, not available this session. As a cheaper substitute, tested
whether the purging signal is explained by **within-panel call confidence** instead: the founder VCF
(`panel/arch3/chr*/merged_231_chr*_final.vcf.gz`, Minigraph-Cactus) carries per-record `F_MISSING`
(fraction of 231 founders missing a genotype) and `MA` (alleles missing in panel haplotypes); `CONFLICT`
turned out to be unpopulated genome-wide (0/2.25M records) and was dropped. Two checks
(`_extract_vcf_callqual.sh` + `_sv_callqual_artifact.py`):
- **(A) Genome-wide, common SVs (MAC≥12, n=7,940): insertions are NOT worse-called than deletions** —
  if anything slightly *better* (F_MISSING median 0.0498 ins vs 0.0584 del, MWU p=3e-19; MA identical,
  p=0.75). No support for "insertions are just harder to call."
- **(B) The climate-slope β purging signal does not track call quality.** Splitting the 2,488 scored
  insertions into well- vs poorly-supported (F_MISSING≤median & MA=0, ~50/50 split), the bio1 and bio18
  KS-vs-matched-SNP purging signal is essentially the **same effect size in both tiers**
  (bio1 median β: −0.01034 well vs −0.01015 poor, both p≪1e-15; bio18 similarly both tiers p≪1e-12).
  Continuous check confirms it: Spearman(F_MISSING, β) = −0.012 (p=0.56), Spearman(MA, β) = −0.010
  (p=0.64) among insertions — a calling artifact predicts a real negative correlation (worse call →
  more spuriously negative β) and there is none.
- **Net: the insertion-polarity signal is NOT a within-panel calling-confidence artifact.** This rules
  out the cheapest version of the artifact hypothesis but is **not** true ancestral polarization — it
  can't distinguish real derived-insertion selection from a *systematic* reference/mapping bias that
  the Minigraph-Cactus graph itself wouldn't flag as low-confidence (e.g. a consistent bias in how
  insertions vs deletions are represented/genotyped against TAIR10, present even in well-supported
  calls). That systematic-bias version of the caveat still needs a true outgroup for a full resolution.
  Artifacts: `_extract_vcf_callqual.sh` (bcftools extraction), `_sv_callqual_artifact.py` (join + test),
  `analysis/grenenet_gea/archive/window_hapfreq_retired/sv_adaptive/results/{vcf_callqual_chr*.tsv,sv_callqual_artifact.npz}`.

## Artifacts
- **Estimator / compute:** `_temporal_s_plots_snp_vs_nonsnp.py`, `_compute_s_climate_slope.py`,
  `_compute_s_dist_by_stratum.py`, `_audit_s_classes.py`.
- **Other tests:** `_nonsnp_temporal_category.py`, `_temporal_selection_snp_vs_nonsnp.py`,
  `_temporal_sel_drift_maf.py`, `_temporal_s_enrich_initqty.py`.
- **Replicate arm (untouched this session, separate open thread):** `_compute_parallelism.py`
  (parallelism ρ + per-site z/ρ; has the fast `read_rows`), `_picmin.py` (real PicMin),
  `_build_parallelism_nb.py`, `_build_picmin_nb.py`, `_build_parallelism_by_site_nb.py`.
- **Calling-quality artifact check:** `_extract_vcf_callqual.sh`, `_sv_callqual_artifact.py` (see
  Caveats section above).
- **Notebook (2026-07-15, consolidated — the 7 former separate `s_*` notebooks and their `_build_*_nb.py`
  scripts, plus the standalone de-trending audit script, were merged and removed after the audit found
  none of them redundant with each other):** `analysis/grenenet_gea/notebooks/temporal_s_consolidated.ipynb`,
  built by `_build_temporal_s_consolidated_nb.py` (reads `s_dist_by_stratum.npz/.csv` +
  `s_climate_slope.npz`/`_sign_by_site.csv`; no new compute). Sections: (1) whole-distribution
  median shift, de-trended — null; (2) tail-specific views (histogram/ECDF/ECDF-difference/
  shift-function/vs-climate) — real; (3) climate-slope β — real; (4) direct de-trending audit of
  (2)-(3) — survives. `basic` env.
  **Replicate arm (separate, untouched):** `parallelism.ipynb`, `picmin.ipynb`,
  `parallelism_by_site.ipynb`.
- **Figures / data** (`analysis/grenenet_gea/archive/window_hapfreq_retired/sv_adaptive/results/`): `s_climate_slope_{bio1,bio18}.png`,
  `s_purging_intensity.png`, `s_ecdf_difference_by_site.png`, `s_vs_climate.png`, `s_climate_slope.npz`,
  `s_climate_slope_sign_by_site.csv`, `temporal_*.csv`; **replicate arm:** `parallelism.npz`,
  `picmin.npz`, `parallelism.png`, `picmin.png`, `parallelism_by_site{,_summary,_excess}.png`.

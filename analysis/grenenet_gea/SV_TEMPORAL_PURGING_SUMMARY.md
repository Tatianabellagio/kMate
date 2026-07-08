# SV temporal purging — thread summary (2026-07-02)

> ## ⚠️ SUPERSEDED (2026-07-08) — the signal did NOT survive
> The estimator was regenerated under **`--unit chrom` + full-panel Kf_w**
> (`--normalize per_founder`). On the rigorous frequency-de-trended metric
> (`s_distribution_by_site`) the **SV excess-vs-baseline is essentially ZERO:
> median +0.0011, negative at only 15/31 sites (a minority; sign test ≈ p=1.0,
> n.s.).** There is **no genome-wide SV purifying excess**; any residual is
> confined to a few of the hottest gardens. The old conclusion below — "~22/31
> sites negative, a modest excess of purifying selection, climate-graded, three
> methods agree" — **NO LONGER HOLDS.** The raw / single-tail methods
> (s_histogram, s_vs_climate, ecdf, shiftfunction) still show apparent purging
> only because they keep the frequency / founder-projection confound; the
> **de-trended `s_distribution_by_site` is the authoritative view (signal null).**
> The historical text is retained below for the record, marked superseded.

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

### [SUPERSEDED] Cleanest test — per-variant climate slope β = d s / d climate  (`_compute_s_climate_slope.py`)
> **This arm is frequency-confounded and no longer supports a climate-graded purging claim.** On the
> de-trended `s_distribution_by_site` the signal is null (see banner). The per-axis correlations are
> also weaker/mixed than stated below: aridity **bio18 corr = +0.54**, temperature **bio1 corr = −0.16**
> (so "rises toward harsh sites" is false on the bio1 axis). Table retained for the record.

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

**[SUPERSEDED] The "three methods agree" claim does not survive the Kf_w / `--unit chrom` rerun.**
The three arms (climate-slope β, ρ-parallelism, PicMin) all share the frequency / founder-projection
confound that the de-trended `s_distribution_by_site` removes — and on that authoritative metric the
SV excess-vs-baseline is ≈0 (median +0.0011, 15/31 sites, sign test n.s.). There is **no genome-wide
SV purifying/climate-graded excess**; any residual is confined to a few of the hottest gardens. The
historical multi-method text is retained above for the record only.

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
  `results/grenenet_gea/sv_adaptive/{vcf_callqual_chr*.tsv,sv_callqual_artifact.npz}`.

## Artifacts
- **Estimator / compute:** `_temporal_s_plots_snp_vs_nonsnp.py`, `_compute_s_climate_slope.py`,
  `_compute_s_dist_by_stratum.py`, `_audit_s_classes.py`.
- **Other tests:** `_nonsnp_temporal_category.py`, `_temporal_selection_snp_vs_nonsnp.py`,
  `_temporal_sel_drift_maf.py`, `_temporal_s_enrich_initqty.py`.
- **Replicate arm:** `_compute_parallelism.py` (parallelism ρ + per-site z/ρ; has the fast `read_rows`),
  `_picmin.py` (real PicMin), `_build_parallelism_nb.py`, `_build_picmin_nb.py`,
  `_build_parallelism_by_site_nb.py`.
- **Calling-quality artifact check:** `_extract_vcf_callqual.sh`, `_sv_callqual_artifact.py` (see
  Caveats section above).
- **Notebooks** (`analysis/grenenet_gea/notebooks/`): `s_climate_slope.ipynb` (β for bio1 & bio18 +
  purging-intensity), `s_ecdf_difference_by_site.ipynb`, `s_ecdf_by_site.ipynb`, `s_vs_climate.ipynb`,
  `s_histogram_by_site.ipynb`, `s_shiftfunction_by_site.ipynb`, `s_distribution_by_site.ipynb`;
  **replicate arm:** `parallelism.ipynb`, `picmin.ipynb`, `parallelism_by_site.ipynb`. `basic` env.
- **Figures / data** (`results/grenenet_gea/sv_adaptive/`): `s_climate_slope_{bio1,bio18}.png`,
  `s_purging_intensity.png`, `s_ecdf_difference_by_site.png`, `s_vs_climate.png`, `s_climate_slope.npz`,
  `s_climate_slope_sign_by_site.csv`, `temporal_*.csv`; **replicate arm:** `parallelism.npz`,
  `picmin.npz`, `parallelism.png`, `picmin.png`, `parallelism_by_site{,_summary,_excess}.png`.

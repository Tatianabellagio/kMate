# SV climate-GEA — literature review & recommended analysis plan

*GrENE-net structural-variant climate adaptation. Synthesis of four research sweeps,
2026-06-04. All claims traced to primary sources (refs at end). Skeptical by design:
flags where evidence is weak or where our design departs from published benchmarks.*

---

## 0. Bottom line up front

1. **SVs are well-established climate-adaptation targets** — strongest for **inversions and
   TE insertions**; in *Arabidopsis* the clean cases cluster on flowering-time / vernalization
   (FRI/FLC) and stress loci [Kang 2023; Fransz 2016; Baduel 2021; Todesco 2020].
2. **SVs add independent signal SNPs miss — but only for some classes.** Biallelic deletions are
   largely SNP-tagged; **inversions / multi-allelic / CNV SVs carry ~30–50% independent signal**
   (only 54% of adaptive human SVs in strong SNP LD [Yan 2021]; only ~17% of *Arabidopsis* SVs
   tagged at LD>0.6 [Kang 2023]). Frame SVs as **complementary**, not a replacement.
3. **Our specific study is genuinely novel.** No published **pool-seq, multi-SV-class,
   common-garden TEMPORAL climate-GEA** exists. SV-GEA to date is static individual-genotype
   pangenomes [Kang 2023; Murray 2025] or inversion-only pool-seq clines [Kapun 2016]; GrENE-net's
   own pipeline is **SNP-only** [Exposito-Alonso 2026]. We occupy an unfilled combination.
4. **Our statistical wall is real and correctly diagnosed.** GIF ≈ √(n_pools/n_sites) is the
   textbook **pseudoreplication** signature; the climate predictor is a **Level-2 (site-constant)**
   variable, so the effective N is **~20–30 sites, not pools**. LFMM's GIF calibration *revealed*
   this power limit; it didn't create it. The `(1|site)` degeneracy (84% p=1) is the standard
   multilevel result (a group random intercept competes with a group-level fixed effect)
   [Snijders & Bosker 2012; Gelman & Hill 2007].
5. **The way forward is units + the temporal axis, not a fancier GEA.** Test climate at the
   **site level** (means-as-outcomes / two-stage); lean on **Δp = within-garden change**, which
   differences out the shared seed-mix starting structure; use a **replicate-aware,
   drift-corrected temporal test (ACER)** with a **site-level permutation null**; **aggregate**
   for polygenic / biological-set signal; and **validate by cross-garden parallelism**. Accept
   N≈20 and don't over-correct (structure is weak; collinearity VIF≈1.19 is not the binding
   constraint).

---

## 1. SVs and climate / local adaptation

**Strongest evidence = inversions & TEs.** *Arabidopsis*: the Chr4 1.17-Mb *Vandal*-derived
inversion associates with drought fecundity **but sits in LD with Col-FRIGIDA** — the canonical
warning that an SV's climate signal can be a flowering-time hitchhiker [Fransz 2016]; TE insertions
into/near *FLC* and stress genes are repeatedly adaptive [Baduel 2021; Quadrana 2016; Nat Plants
2024]; the 32-genome pan-genome links SV-bearing genes to temperature variables (BIO2/BIO7) and
names candidate adaptive SVs (HPCA1/KNAT3 TE insertions) [Kang 2023] — though those are
candidate-/enrichment-level, not frequency-resolved. *FRI* deletions underlie the classic
flowering-time latitudinal cline [Johanson 2000]. Other systems: sunflower **37 inversion
haploblocks × soil/climate** [Todesco 2020]; *Drosophila* In(3R)Payne latitudinal cline [Kapun
2016]; *Mimulus* life-history inversion [Lowry & Willis 2010]; stickleback *Eda* + inversions.

**SV vs SNP (do SVs add signal?).** Class-dependent. Biallelic deletions ≈ SNP-redundant;
duplications/inversions/multi-allelic SVs poorly tagged → independent targets [Saitou 2021; Yan
2021; Mérot 2020]. Contested: some domestication/pangenome studies find adaptation
"predominantly SNP-driven" with SVs purged as deleterious. **Defensible claim:** SVs are
complementary, independent at a non-trivial minority (~30–50%) of loci.

**Interpretive trap for us:** flowering-time confounding (FRI/FLC, Chr4 inversion). Any
temperature-associated SV near vernalization/flowering loci needs an explicit "causal vs linked"
caveat.

---

## 2. How novel is a pool-seq SV climate-GEA? (SV-specific methodology)

**SV-GEA is an emerging frontier, not an established subfield.** The GEA/selection apparatus was
built and validated on **SNPs**; when SVs are run through it, researchers **borrow SNP frameworks
unchanged** (RDA, LFMM, BayPass, Fst). There is **no SV-tailored GEA or temporal-selection
method**. Direct SV-GEA cases are recent and sparse: Eucalyptus RDA → CHILL1 [Murray 2025];
*Arabidopsis* SV-GWAS vs 21 climate vars [Kang 2023]; TE-copy-number vs WorldClim [Baduel 2021];
all **individual-genotype, static** panels. Pool-seq SV work is **inversion-only clines** [Kapun
2014/2016].

**The unoccupied combination = ours:** pool-seq SV allele frequencies (incl. CNV/TE/PAV/indel, not
just inversions) → **climate GEA AND temporal selection** in a **replicated common-garden time
series**. GrENE-net itself is SNP-only [Czech 2022; Exposito-Alonso 2026]. → "to our knowledge,
the first" is defensible (caveat: bioRxiv keyword search was limited, so phrase as "we are not
aware of," not "first ever").

**SV-specific pitfalls the literature insists on (and our QC must address):**
- **Defining an SV "allele frequency" in pools** is a known hard problem; the equal-allele-sampling
  assumption of standard pool-seq AF estimators is **false for SVs** → Kessner-Turner-Novembre's
  **EM estimator** [Kessner 2013] is the key precedent; our k-mer/kMate approach must own this.
- **Short-read SV calls: up to 91% false positive** [David 2024]; the persuasive orthogonal QC is
  **showing SV frequencies recapture known population structure** (and ideally SNP-validated
  GrENE-net signals).
- **Repeat/centromere artifacts & LD-clustering of co-located SVs** producing **identical
  frequency trajectories** — *exactly our Chr3 identical-β clusters.* Standard fix: mask
  repeats/centromeres; **collapse co-located SVs into haploblocks** before counting hits.
- **Low-frequency & multi-allelic SVs are least reliable** → MAF / min-evidence filters; treat
  rare-SV hits skeptically (our Chr3:6.36 Mb, p₀≈0.07 cluster).
- **Genotyping-uncertainty propagation into GEA is an unmet gap** — a potential methodological
  differentiator if we do it.

---

## 3. GEA models & what to control for

**Method menu & how each corrects** (full table in refs): **LFMM2** (K latent factors from the
genotype matrix, GIF recalibration) [Caye 2019]; **RDA / partial RDA** (multivariate; pRDA
*conditions out* supplied structure/geography and **partitions variance** into pure-climate /
pure-structure / confounded) [Forester 2018; Capblancq 2021]; **BayPass / BayEnv2** (population
covariance Ω, env tested conditional on Ω) [Gautier 2015; Coop 2010]; **BayeScEnv** (Fst-env);
**GEMMA** (kinship random effect) [Zhou & Stephens 2012]; **Gradient Forest** (nonparametric
turnover, no significance test); **WZA** (window weighted-Z) [Booker 2024]. Benchmark ranking:
LFMM ≈ BayPass > BayeScEnv ≈ BayEnv2 [Gautier 2015]; RDA best multivariate TP/FP [Forester 2018].

**What MUST be controlled** [Rellstab 2015; Hoban 2016; Lotterhos & Whitlock 2014/2015]: population
structure, demography/allele-surfing, isolation-by-distance/spatial autocorrelation, relatedness,
and the **climate–structure collinearity** (over-correction erases true signal). **For GrENE-net
the classical confounders are largely moot** (one seed mix → little deep structure; no range
expansion in 3 generations). Our real threats, in order: **(1) pseudoreplication / Level-2
identification, (2) drift among plots mimicking signal, (3) climate–structure collinearity** —
and most off-the-shelf GEA tools silently mishandle (1).

**Choosing K:** sNMF cross-entropy [Frichot 2014], scree/Tracy-Widom, and the practical decider
**GIF/p-value-histogram calibration** [Caye 2019; François 2016]. **When GIF won't reach 1** (us):
that means inflation is **not low-rank structure** — adding K won't fix a pseudoreplication
(variance-of-the-mean) GIF and **will erase real signal**. Fix the *unit* first.

**Pseudoreplication & the Level-2 predictor** [Hurlbert 1984; Millar & Anderson 2004; Snijders &
Bosker 2012]: effective N = n_sites. `log(p_t/p₀) ~ climate + (1|site)` is **unidentified** (random
intercept = site mean = the thing climate explains) → p piles at 1 (our 84%). **Remedies:**
**means-as-outcomes** (aggregate to per-site value, regress ~20–30 sites on climate); **two-stage**
(per-site effect with `(1|plot)`, then between-site meta-regression); cluster-robust SEs only with
CR2/wild-bootstrap given <40 clusters.

**WZA** boosts power "when a small number of demes are sampled" [Booker 2024] — relevant — but is
**awkward for sparse SVs** (few SVs per window); aggregate over **biological sets** (SVs in
flowering/stress genes) instead.

---

## 4. Selection-scan statistics — temporal is our lever

**Classical families** (Fst-outlier [BayeScan, OutFLANK, pcadapt], haplotype/EHH, SFS/CLR) are
**mostly inapplicable**: pool-seq has no phase, SVs are sparse, 3 generations is too short for
sweeps, and SFS/CLR are demography-fragile. Only the **environmental-association** family fits the
gradient design — and we've already used its best-calibrated member (LFMM2).

**Temporal / E&R — the relevant family:**
- **CMH across replicates** is the E&R workhorse but its binomial null **ignores drift + pool-seq
  variance** → wildly anticonservative [Spitzer 2020]. *This is our Kendall/GIF inflation
  re-expressed.*
- **ACER** [Spitzer 2020] = adapted χ²/CMH that **injects drift (Ne, generations) + pool-seq
  variance** into the null. **The single best off-the-shelf, replicate-aware, pool-seq-native
  temporal test for us.** (github.com/MartaPelizzola/ACER)
- **Buffalo & Coop temporal covariance** [Buffalo & Coop 2020] — under drift, covariance of Δp
  between intervals = 0; positive covariance = selection. A **genome-wide diagnostic** (needs ≥3
  timepoints, which we have), complements ACER. (cvtkpy)
- **WFABC / ApproxWF** — per-locus s from a trajectory; stretched by a 3-generation series.

**Climate-*dependent* directional selection** ("rises-in-warm / falls-in-cold") = the antagonistic
signature [Exposito-Alonso 2019; Machado 2021]. **Proper test:** per-SV relationship of Δp
direction/magnitude vs site climate, with the **null built by permuting climate labels across the
~20–30 gardens** (site-level permutation [Machado 2021]) — the only generally valid null when demes
share an environment. **Condition on founding frequency + climate-of-origin** or you detect sorting
of pre-existing structure, not in-situ selection.

**Calibration:** GIF is "an honesty tool, not a power tool." Use **FDR** (not Bonferroni) — but only
on a **valid** null (ACER drift-aware or site-permutation), never on the inflated Kendall p's.

**Power reality:** standard GEA designs use **~40 demes** and still find modest per-locus power
[Lotterhos & Whitlock 2015; Booker 2024]; below ~20 independent environments, per-SV genome-wide
climate hits are **not realistically expected** — consistent with our calibrated-LFMM null.
**Replication buys temporal, not environmental, df** — so a replicate-aware temporal test (ACER)
beats per-SV GEA; the climate-correlated ceiling is still ~20 sites.

---

## 5. Recommended analysis plan for GrENE-net SV climate-GEA

Ordered, and consistent across all four sweeps:

1. **QC first (pre-empts the #1 reviewer attack).** Show SV frequencies **recapture known
   accession/garden structure** [David 2024]; mask repeats/centromeres; **collapse co-located SVs
   into LD haploblocks** (kills the Chr3-type identical-trajectory artifacts); apply MAF/min-count
   filters; flag multi-allelic/copy-number SVs as low-confidence.
2. **Make the unit the site.** Replace per-pool tests with **site-level means-as-outcomes** or a
   **two-stage** model (`(1|plot)` within site → between-site regression on climate). Accept
   **N≈20–30**.
3. **Primary per-SV scan = ACER** (drift- + pool-aware replicate temporal test) on Δp across
   replicate gardens — the correctly-calibrated replacement for the naive Kendall.
4. **Climate-correlation test = Δp-vs-climate with a SITE-LEVEL PERMUTATION null** (Machado
   template); condition on founding freq + climate-of-origin.
5. **Genome-wide diagnostic = Buffalo–Coop temporal covariance** (selection vs drift, ≥3
   timepoints).
6. **Aggregate for polygenic / set-based signal** — partial-RDA variance partition (climate vs
   structure vs confounded) + biological-set weighted-Z (flowering/stress genes), not fixed-window
   WZA.
7. **Don't over-correct structure** (weak; VIF≈1.19). Report the confounded fraction as undecidable.
8. **Validate by cross-garden parallelism** (selection synchronized across climate-matched gardens)
   — a stronger argument than any single scan p-value, and GrENE-net's intended design strength.
9. **Expect & quantify SNP-redundancy** — per SV hit, report LD/r² with the best nearby SNP and the
   SNP-GEA signal there; frame SVs as complementary.

**Honest headline for the poster/paper:** *"Naive per-SV SV climate-GEA looks striking but is
pseudoreplication; with ~20 independent climates the calibrated genome-wide signal is weak. The
power lies in the temporal axis — a replicate-aware, drift-corrected scan of allele-frequency
change validated by cross-garden parallelism — not in chasing individually-significant SVs."*

---

## Key references

**SVs & adaptation:** Mérot 2020 TREE 35:561; Wellenreuther/Catanach 2019 Mol Ecol 28:1331; Fransz
2016 Plant J 88:159; Kang 2023 Nat Commun 14:6259; Baduel 2021 Genome Biol 22:138; Quadrana 2016
eLife 5:e15716; Johanson 2000 Science 290:344; Todesco 2020 Nature 584:602; Kapun 2016 MBE 33:1317;
Lowry & Willis 2010 PLoS Biol 8:e1000500; Saitou 2021 MBE msab313; Yan 2021 eLife 10:e67615; David
2024 GBE 16:evae049; Pokrovac & Pezer 2022 Front Genet 13:1060898; Murray 2025 bioRxiv
2025.12.18.695066; Kessner 2013 MBE 31:1138.

**GEA methods & confounders:** Frichot 2013 MBE 30:1687; Caye 2019 MBE 36:852; Forester 2018 Mol
Ecol 27:2215; Capblancq & Forester 2021 MEE 12:2298; Gautier 2015 Genetics 201:1555; Coop 2010
Genetics 185:1411; Zhou & Stephens 2012 Nat Genet 44:821; Booker 2024 Mol Ecol Resour 24:e13768;
Rellstab 2015 Mol Ecol 24:4348; Hoban 2016 Am Nat 188:379; Lotterhos & Whitlock 2014 Mol Ecol
23:2178 & 2015 24:1031; Frichot 2014 Genetics 196:973; François 2016 Mol Ecol 25:454.

**Selection scans / temporal / pseudoreplication:** Spitzer 2020 Ann Appl Stat 14:202 (ACER);
Buffalo & Coop 2020 PNAS 117:20672; Foll 2015 Mol Ecol Resour 15:87 (WFABC); Ferrer-Admetlla 2016
Genetics 203:831 (ApproxWF); Machado 2021 eLife 10:e67577; Whitlock & Lotterhos 2015 Am Nat 186:S24
(OutFLANK); Luu 2017 Mol Ecol Resour 17:67 (pcadapt); Hurlbert 1984 Ecol Monogr 54:187; Millar &
Anderson 2004 Fish Res 70:397; Snijders & Bosker 2012 *Multilevel Analysis* 2e; Gelman & Hill 2007;
Devlin & Roeder 1999 Biometrics 55:997 (GIF).

**Study context:** Exposito-Alonso 2019 Nature 573:126; Czech 2022 (GrENE-net pilot); Exposito-Alonso
& GrENE-net consortium 2026 Science (adz0777) / bioRxiv 2025.05.28.654549.

*Sourcing caveat: a handful of DOIs/author lists were confirmed by title/venue via search snippets,
not read end-to-end — verify before citing verbatim in print.*

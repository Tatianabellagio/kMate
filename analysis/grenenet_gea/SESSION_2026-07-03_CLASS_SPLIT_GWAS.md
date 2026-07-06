# Session log — 2026-07-03 · SNP vs non-SNP: variance, peaks, and candidate genes

Branch `add-kmate`. Env: `kmate` (compute) / `basic` (notebooks — matplotlib hangs in `plotting`).
All heavy work on a Savio compute node (`n0037.savio4`), never the login node.
Sits inside the "does kMate's non-SNP (indel+SV) layer add adaptive signal beyond SNPs" thread —
see `VAREXP_SELECTION_HANDOFF.md` and memory `project-gea-varexp-snp-nonsnp`.

The through-line of the whole day: **at every resolution we looked — genome-wide kinship,
SNP-untagged kinship, per-marker peaks, per-site and multi-site GWAS, across all bioclim axes —
SNPs and non-SNP markers tell essentially the same story. The non-SNP layer adds almost nothing;
what little it flags uniquely is scattered, marginal, and mostly sub-Bonferroni.**

---

## 1. SNP-untagged non-SNP GRM — the "kMate gain" null, made non-tautological

Prior result: genome-wide K_snp vs K_nonsnp are ~redundant (corr 0.998), non-SNP adds no
selection-trait variance. Obvious objection: "of course — both are just genome-wide relatedness."
This session killed that objection.

- `build_nonsnp_tagging.py` → `results/grenenet_gea/varexp/nonsnp_tagging_chr{1..5}.npz`. For each
  of the 525k non-SNP markers feeding K_nonsnp (MAC≥5, call≥90%), max founder-LD r² vs any panel
  SNP within ±50 kb. **99.7% of non-SNP markers are SNP-tagged at r²≥0.2; only 1,383 (0.26%) are
  SNP-invisible.** (At r²≥0.5: 8% untagged.)
- `build_untagged_grm.py` → `untagged_grms.npz`. GRMs from only the SNP-untagged markers:
  `K_untagged_r02` (1,383 markers, corr with K_snp = **0.396** — genuinely decorrelated) and
  `K_untagged_r05` (42,140 markers).
- `varexp_untagged.py` → `varexp_untagged.csv`. Marginal h² / joint 2-GRM REML+LRT / GBLUP CV gain
  vs K_snp, over 4 axes (global + 3 bio1-tercile zones) × {raw, RINT} × {r02, r05}.

**Result:** even the strict SNP-invisible GRM (r²<0.2) has marginal h²=0.889 (unlinked markers still
recover the same founder clades) but **conditional on K_snp: joint LRT p=0.35, CV gain +0.0004
(trivial).** Robust across all 16 cells (gain never |>0.004|, flips sign). → The "no gain"
conclusion is now strong, not a construction artifact: even markers SNPs *cannot* tag only
re-derive the same genealogy, which SNPs already predict ~89–94% of.

---

## 2. Class-split single-marker GWAS — do SNP / non-SNP / SV scans find the same peaks?

`class_split_gwas.py` → `class_gwas_{snp,nonsnp,sv}.npz` + `class_gwas_summary.json`.
Single-marker LOCO-EMMAX on the locked selection trait `s` (212 founders), **one shared LOCO
kinship correction per chromosome** (so only the tested marker class differs), per-site z stacked
into a Bolormaa multi-trait meta → **JOINT** (31 df, any-site), **GLOBAL** (1 df, generalist),
**CLIMATE** (1 df, bio1-differential). Three test sets: SNP (1.75M), non-SNP/indel+SV (536k),
strict SV-only (13.8k). Efficient: one eigendecomposition + rotation per (chrom,class) reused
across all 31 sites; full genome × 31 sites in ~5.5 min.

**Numerical bug found + fixed (important):** ~7–10 markers/chrom are near-*fixed* (202–207/212
founders carry the allele — MAC treats both tails symmetrically so it looks like a rare case but
isn't). After LOCO rotation their genotype vector is ~collinear with the intercept, so the GLS
determinant `Saa·Sbb − Sab²` computes as pure float64 rounding noise (~1e-16 relative — machine
epsilon, vs >1e-6 for every real marker). Uncaught, these blow up to |z|~1e19 and wreck the
cross-site corr(Z) for *every* marker (λ_JOINT collapsed to 0, "38% Bonferroni-significant").
Fixed by masking markers with relative det ≤1e-8 as untestable (clean 10-order gap; <0.002% loss).
Standard near-monomorphic-after-conditioning GWAS hygiene.

**Findings:**
- **Variance explained is identical across classes** (independent confirmation of §1's kinship
  result via a completely different, per-marker method): mean per-site genomic-control **λ = 1.040
  for both SNP and non-SNP**; median per-site h² 0.951 vs 0.952; median CV R² 0.549 vs 0.545.
- **Peaks substantially concordant:** the single strongest JOINT hit genome-wide is the **same
  20 kb window (Chr2 ~2.34 Mb) for all three classes independently.** SNP-vs-nonSNP per-window
  −log10p correlates Spearman ρ=0.545 (p≈0); top-hit window Jaccard 0.67 at top-0.1% → 0.29 at
  top-1%. SNP-vs-SV weaker (ρ=0.14) — expected, SV set is ~135× smaller.
- **CLIMATE is null** for bio1 (and GLOBAL is null) in every class — consistent with the
  established climate-gradient null.

---

## 3. Notebooks — overlap tables (no Manhattans embedded)

Built the peak-overlap table the user asked for (mirroring the phase-1 kendall/lfmm/binomial
`overlap_df`), at the **clq0.9 block level**, with **both Bonferroni and FDR** counts:

- `notebooks/class_gwas_multitrait.ipynb` — rows = JOINT, GLOBAL, and **CLIMATE for all 19 bioclim
  vars + PC1 of all bioclim** (re-derived cheaply from the saved per-site Z; no GWAS rerun).
  Columns: `n_snp_bonf / n_snp_fdr / n_nonsnp_bonf / n_nonsnp_fdr / n_shared_fdr / jaccard_fdr /
  pct_*`. JOINT: 21 Bonf / 1364 FDR (snp) vs 9 / 104 (nonsnp), 82 shared, 79% of non-SNP hits are
  also SNP hits. GLOBAL + CLIMATE_bio1: 0/0. **Several temperature bioclim vars show FDR-only SNP
  hits (bio10=140, bio11=97, bio13=72, bio16=64) — but these mostly vanish under Bonferroni, don't
  replicate in non-SNP, and are NOT corrected across the 20 collinear climate axes → treated as
  FDR-tail noise, not real climate adaptation.**
- `notebooks/class_gwas_persite.ipynb` — same table but one row per garden (per-site scans have no
  cross-site CLIMATE contrast). Real hits at a few gardens (site 4: 5 Bonf / 409 FDR snp, 24 shared
  with non-SNP), many gardens null. Plus per-site λ-by-class plot and per-site top-hit table.

Bug fixed en route: `lib.collapse_to_blocks` selects the block lead by MAX |stat|; passing `-p`
made `abs(-p)=p` pick the *least* significant marker. Switched to `-log10(p)`. (Verified counts
then matched the manual BH.)

---

## 4. All Manhattan + QQ plots → PNG folder (not in notebooks)

`plot_class_gwas_pngs.py` → **318 PNGs** in `results/grenenet_gea/varexp/gwas_plots/`:
- multitrait (132): JOINT/GLOBAL/CLIMATE×(bio1–19+PC1) = 22 contrasts × 3 classes × {manhattan, qq}
- per-site (186): 31 gardens × 3 classes × {manhattan, qq}

Every Manhattan draws **both threshold lines at the BLOCK level** (matching the overlap tables —
this was a fix: originally marker-level, so the FDR line silently vanished whenever marker-level
FDR found nothing even though block-level FDR had hits). Legend shows the significant-block count
on each line, with an explicit "(0 blocks)" entry when FDR is empty.

---

## 5. Non-SNP-only candidate genes + API descriptions

The scientifically pointed question: **what does the non-SNP scan flag that SNPs miss, and are any
of those genes interesting?**

- `nonsnp_only_genes.py` → `nonsnp_only_{blocks,genes}.csv`. For every one of the 53 contrasts,
  the clq0.9 blocks non-SNP-FDR-significant but NOT SNP-FDR-significant (and a Bonferroni-strict
  subset). Blocks → TAIR10 genes within span **plus ±2 kb promoter flank**. Result: **53 blocks
  (8 at Bonferroni) → 126 genes (71 in-block + 55 flank-only).**
- `nonsnp_only_genes_describe.py` → `nonsnp_only_genes_described.csv`. Descriptions via **Ensembl
  Plants REST** (batch POST /lookup/id: symbol + description + biotype) and **UniProt REST**
  (per-gene: protein name + curated `cc_function`; 36/126 have one). Compute node reaches both.
- `_build_nonsnp_only_genes_nb.py` → `notebooks/nonsnp_only_genes.ipynb`: Bonferroni subset,
  keyword theme screen (flowering/cold/heat/circadian), full annotated table.

**Interesting genes by theme (all hypothesis-level — see caveats):**
- **Circadian + flowering:** **GIGANTEA (GI, AT1G22770)** — clock/photoperiod/freezing regulator.
  BUT: garden-4 only, ±2 kb *flank* of block Chr1_4247 (~1.5 kb upstream; block's real in-block
  genes are AT1G22780/90/800), one marker at −log10p≈4.82 — just over FDR (4.52), well under
  Bonferroni (5.95). See `gwas_plots/GI_locus_site4_nonsnp.png` (annotated genome-wide + Chr1 zoom).
- **Cold:** **ADS2 (AT2G31360)** — Δ9 fatty-acid desaturase (membrane cold-acclimation). **Best-
  motivated: in-block, multitrait JOINT** (the one well-calibrated contrast). Also SEX1/GWD
  (starch mobilization / freezing), PI-4KBETA2.
- **Heat:** HSP70 (AT4G16660), HIP1 (HSP70-interacting), a Clp-N chaperone.
- **ABA/drought (not requested, but the only Bonferroni-level themed hit):** **ERA1 (AT5G40280)**,
  farnesyltransferase β; plus XERICO, NAC032 (stress/senescence TF, JOINT).

**Caveat carried throughout:** this is the *fragile tail* by construction — blocks one class calls
and the other doesn't. Most are FDR-only, per-site, and/or flank hits. The CLIMATE-contrast genes
inherit the bioclim-null caveat (§3). ADS2 (in-block, JOINT) is the sturdiest; GI is the most
eye-catching but statistically weak. A hypothesis-generating list, not confirmed loci.

---

## File index (all created/modified today)

Scripts (`analysis/grenenet_gea/`): `build_nonsnp_tagging.py`, `build_untagged_grm.py`,
`varexp_untagged.py`, `class_split_gwas.py`, `plot_class_gwas_pngs.py`, `plot_gi_locus.py`,
`nonsnp_only_genes.py`, `nonsnp_only_genes_describe.py`, and notebook builders
`_build_class_gwas_multitrait_nb.py`, `_build_class_gwas_persite_nb.py`,
`_build_nonsnp_only_genes_nb.py`.

Notebooks (`analysis/grenenet_gea/notebooks/`): `class_gwas_multitrait.ipynb`,
`class_gwas_persite.ipynb`, `nonsnp_only_genes.ipynb`.

Outputs (`results/grenenet_gea/varexp/`): `nonsnp_tagging_chr{1..5}.npz`, `untagged_grms.npz`,
`varexp_untagged.csv`, `class_gwas_{snp,nonsnp,sv}.npz`, `class_gwas_summary.json`,
`nonsnp_only_{blocks,genes,genes_described}.csv`, `gwas_plots/` (318 PNGs + `GI_locus_site4_nonsnp.png`).

## Open / optional next steps
- Cross-axis multiple-testing correction on the bioclim CLIMATE hits (currently 20 independent BH
  runs on collinear temperature vars) to confirm they're noise.
- Pull the actual indel/SV variants inside the ADS2 and GI blocks (size, freq, position vs gene).
- Same annotated locus plot for ADS2 (the sturdiest candidate) as was done for GI.
- Census-flavour trait + 3-way SNP/indel/SV variance breakdown (still open from the earlier handoff).

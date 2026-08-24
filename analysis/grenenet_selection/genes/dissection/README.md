# genes_expl — high-resolution locus dissection of GEA new-peak genes

Per-gene "better-resolution" follow-up of hits flagged in the multiaxis LFMM-WZA
(`phase1_replication/multiaxis/newpeak_dotgrid_lfmm_sv.ipynb`). Each locus gets a
cam5-style combined figure — GEA Manhattan + founder-haplotype panel + founder-LD
triangle, sharing a genomic x-axis — so we can see *where the signal really lives*
relative to the genes the block was attributed to.

Pattern mirrors `cam5_replication/pang_cam5/` but rebuilt in matplotlib (project
plotting convention: no titles) instead of R/gggenomes.

## CARK8 / CARK9 (AT2G30740 / AT2G30730) — `plot_cark_combined.py`

`cark_combined_manhattan_hap.png/.pdf`, window **Chr2:13,090,000–13,137,000 (~47 kb)**,
climate axis **bio1**.

Three stacked panels (shared x):
1. **LFMM GEA Manhattan** — cohort `wza_in` variants (SNP/indel/SV), −log10 p,
   coloured by **Δp climate = mean over gen 1–3 of (hot − cold) allele freq**
   (red = allele higher in hot sites); per-class genome-wide Bonferroni lines;
   gene-model track on top (CARK8/9 highlighted).
2. **Founder panel** — 231 accessions ordered cold→hot by home bio1 (left strip =
   home bio1 origin), ALT ticks across the window (SNP grey / indel orange / SV red
   diamond). Over 47 kb every founder is a unique haplotype, so this is the
   per-founder view, not a collapse.
3. **Founder-LD triangle** — r² on common (MAF≥0.05) panel variants + the aggregated
   lead SV; connector drops from the lead SV.

### Key finding: the "CARK SV new-peak" is a block-merge artifact, not a CARK signal
- The block's Bonferroni-clearing signal is a **~5.1 kb insertion at Chr2:13,127,635**
  (SV, LFMM p=3×10⁻⁸ on bio1, MAF 0.05) — **~29 kb downstream of CARK8/CARK9**,
  intergenic, **191 bp from AT2G30810 = GASA12** (gibberellin-regulated) and inside
  the LD block of **AT2G30800 = HVT1** (RNA helicase).
- **Why it was attributed to CARK.** The `wza_in_clq09_tile` block column is
  `merge_small_blocks(assign_tiling(...))` (`blocks_tiling.py`). For the sparse SV
  class (~27k SVs genome-wide) the 14 tiling blocks between the CARK LD-island and
  the SV hold 0–1 SVs each, so `merge_small_blocks` (floor 2 records) **chains them
  backward** into the last SV-dense block — Chr2_5291 (CARK). The SV physically lives
  in tiling block Chr2_5306; the merged **SV-class WZA block spans 40 kb**
  (13,095,207–13,135,224) though its r²≥0.9 LD island is only ~1.4 kb.
- **Founder LD confirms no link:** r²(lead SV, any CARK-block common variant) ≤ 0.05
  (median 0.008); r²(lead SV, its own local block near AT2G30800/810) = 0.88.
- **Δp coloring caveat.** A naive net Δp (final − p0, pooled over sites) is ≈0 for
  this SV and renders white — it is antagonistically sorted (higher in hot at every
  generation: hot−cold = +0.027/+0.027/+0.012 for gen 1/2/3) while being **purged
  everywhere over time** (hot 0.062→0.027, cold 0.034→0.015), so the net cancels.
  The figure colours by the hot−cold **climate divergence** LFMM actually detects.
- **Not a warm-origin haplotype** (unlike cam5's CAM5 hap): founder-panel SV carriers
  are scattered across the home-bio1 range (carrier mean bio1 9.6 vs non-carrier 9.6,
  MWU p=0.64). The climate signal is in the *evolved cohort* AF trajectory, not in
  founder provenance.

### Across-garden dynamics of the lead SV — `plot_sv_garden_dynamics.py`
`cark_sv_garden_dynamics.png/.pdf` + `cark_sv_site_traj.csv`. Per-**site** (garden)
allele-frequency trajectory of the 5.1 kb insertion over the 3 experimental
generations (gen 0 = founding SEEDMIX p0 = 0.026, shared by all gardens), split
COLD vs WARM at the median garden bio1 and coloured by each garden's mean annual
temperature. Resolves the pooled hot/cold signal to all 31 gardens:
- **Warm gardens**: the SV climbs (slope +0.010; several gardens reach AF 0.1–0.25
  by gen 3).
- **Cold gardens**: flat/low (slope −0.002, p=0.48) — the allele stays rare.
This is the per-garden version of the hot−cold divergence that drives the bio1 GEA
hit; the effect is real but modest (not a strong sweep). Data: per-generation AF
matrices `gen_matrices/gen{1,2,3}_nonsnp_af.npy` (a fast column-slice per gen) +
site bio1 from `lib.load_climate()`.

## Functional-context screen of ALL significant SV/indel/non-SNP blocks — `screen_sig_blocks.py`
`sig_block_functional_screen.csv` (all 1,275 unique lead variants) +
`functional_candidates.csv` (the 515 CDS/UTR/promoter ones).

**Question (the CARK lesson generalised):** a block's attributed gene is
meaningless if the lead variant sits far away (tiling+merge artifact). So we take
the significance set — union over {sv, smallindel, nonsnp} × 20 climate axes of the
per-class genome-wide Bonferroni (0.05/n, MAF>0.05) **block-lead** variants (raw
pre-WZA LFMM p; a candidate net, NOT a calibrated hit list) — and classify each at
**its own position** against the full TAIR10 GFF: CDS / 5′-3′ UTR / intron /
promoter (≤1 kb upstream, strand-aware) / TE / proximal-intergenic / gene-desert,
with distance-to-nearest-gene and a `block_gene_mismatch` flag (variant not in the
block's span-genes = a CARK-type re-attribution).

**Result (1,275 unique lead variants; SV counts in parens):**
CDS 60 (18) · UTR 90 (7) · promoter 365 (71) · exon-noncoding 6 (2) · intron 184
(17) · TE 256 (72) · proximal-intergenic 293 (37) · gene-desert 21 (5). So **515
variants (98 SVs)** fall in coding/UTR/promoter space — the functional candidates.
41 are re-attributions (`block_gene_mismatch`); 51 recur across ≥3 axes.
Note TEs claim 72/227 SVs — SVs are heavily TE-associated, as expected.

**Strongest / most robust candidates (recurrent across axes = replicated):**
- `AT1G62630` disease-resistance protein — CDS (MNP), bio15, nlp 13.9, 4 axes.
- `EMB1241` GrpE co-chaperone — 5′/3′UTR indel, 7 axes, nlp 11.9.
- `AT3G22142` lipid-transfer/2S-albumin — **450 bp SV in CDS**, bio15, nlp 11.6.
- `CRK33` cysteine-rich RLK (defense) — **SV promoter**, bio6, nlp 10.9, 6 axes.
- `GPX6` glutathione peroxidase (oxidative stress) — **1.2 kb SV promoter**, 6 axes.
- `WRKY19` TF — **SV in CDS**, 5 axes; `CYP706A7` — **4.6 kb SV in CDS**.
- CARK lead SV is here too, correctly re-placed as **promoter of AT2G30800 (HVT1)**
  with `block_gene_mismatch=True` — the SV *is* ≤1 kb upstream of HVT1, so it does
  have a plausible functional handle (HVT1, not CARK).

The screen includes **all** significant SV/indel/non-SNP block-leads — it does NOT
apply the newpeak SNP-missed filter, so blocks that also carry SNPs are kept.
`augment_snp_cosig.py` then adds, per candidate, whether a genome-wide
Bonferroni SNP (MAF>0.05, any axis) sits within ±2 kb (`snp_cosig_2kb`,
`snp_cosig_dist`, `snp_cosig_axes`) — to separate SNP-shadowed from SV/indel-unique
signals without dropping either. Of the 515 functional candidates, **181 are
SNP-shadowed and 334 are SV/indel-unique** (SV-only: 21 shadowed, 75 unique).

Caveat: raw LFMM p (uncalibrated); rank/triage before wet-lab. Prioritise by tier
(CDS>UTR>promoter), recurrence (`n_axes`), MAF, and SV/indel-uniqueness
(`snp_cosig_2kb`=False); then LD-check the SV↔gene link as done for CARK.
Regenerate: `python screen_sig_blocks.py && python augment_snp_cosig.py`.

### Themed shortlist — `theme_filter.py`
`themed_candidates.csv`: the 123 SV/indel functional candidates (CDS/UTR/promoter,
vclass sv|smallindel) whose gene has a role in CLIMATE / FLOWERING / CIRCADIAN-LIGHT
/ STRESS (matched on annotated protein_name + symbol + TAIR/UniProt categories).
Best-supported (recurrent across axes and/or SV): EMB1241 (GrpE co-chaperone, UTR
indel, 7 axes), AIPP3/RVR1 (Repressor of Vernalization 1, promoter indel, 4 axes),
GPX6 (glutathione peroxidase, 1.2 kb SV promoter, 6 axes), AGL100 (MADS-box TF,
373 bp SV promoter), AFP2 (ABI5-binding/ABA, promoter indel), FAMA (stomatal,
promoter indel), HSP90-4 (heat shock, promoter indel), FRL2/FRL4B (FRIGIDA-like,
flowering), TPT (photosynthesis-acclimation, SV promoter, bio1).

### Per-candidate locus dissection — `plot_locus_combined.py` → `loci/`
Parametrized cam5/CARK-style combined figure (Manhattan Δp + founder-haplotype +
founder-LD + printed LD-confirm of lead↔gene) run on the shortlist. Verdicts and the
robust-vs-fragile triage in `loci/CANDIDATE_VERDICTS.md`. Top picks: **EMB1241, GPX6,
TPT, AGL100** (+ AIPP3/RVR1). Fragile (rare/isolated/low-LD): FAMA, HSP90-4, FRL2,
AT3G22142. Rule: robustness tracks recurrence + founder-frequency + gene-LD, NOT p.

## Data sources
- Founder panel (231 accessions, SV-inclusive, arch3 symbolic-ID merge):
  `panel/arch3/chr2/merged_231_chr2_final.vcf.gz` → region extract `cark_region.vcf.gz`.
- Founder home bio1 (per-accession WorldClim origin):
  `/global/scratch/users/tbellg/gea_grene-net/key_files/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv`
  (`ecotypeid`, `bio1`…). NB: personal-scratch path, not the shared project tree.
- GEA per-variant LFMM p: `phase1_replication/results/multiaxis/wza_in_clq09_tile/lfmm_{snp,smallindel,sv}_gen9_bio1.csv`.
- Δp (hot/cold group-mean AFs): `group_means.npz`.
- Gene models: `lib.load_genes()` (TAIR10). Gene annotation: `phase1_replication/annotate_genes_tair_uniprot.py`.

## Regenerate
```bash
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python   # on a compute node
BCF=/global/home/users/tbellg/miniforge3/envs/kmate/bin/bcftools
$BCF view -r Chr2:13090000-13137000 -e 'AC=0' -Oz -o cark_region.vcf.gz \
    /global/scratch/projects/fc_moilab/tbellg/kmate/panel/arch3/chr2/merged_231_chr2_final.vcf.gz
$BCF index -t cark_region.vcf.gz
MPLBACKEND=Agg $PY plot_cark_combined.py     # -> cark_combined_manhattan_hap.png/.pdf + cark_gea_variants_bio1.csv
```

## Outputs
- `cark_combined_manhattan_hap.png` / `.pdf` — the 3-panel figure.
- `cark_gea_variants_bio1.csv` — window GEA variants (pos, cls, MAF, nlp, Δp_climate).
- `cark_region.vcf.gz` — 231-founder region genotypes (2066 variants).

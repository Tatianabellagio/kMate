# Driver-vs-passenger: are SVs the causal/lead variant in climate-selected blocks?

2026-07-03. Fine-mapping test that beats the LD confound: within each climate-selected
LD block, is the SV the lead/causal variant beyond the best tagging SNP, or a passenger?

> **Consolidated 2026-07-06:** merges the former `HANDOFF.md` (re-run guide + inputs) into this
> file; both originals are in git (`5cefcfa`). Result also summarized in memory `gea-driver-passenger-svs`.

## Design (agreed this session)
- **Blocks = clq0.9** (BigLD LD islands, median ~223 bp / 7 var; `lib.assign_clq_blocks(0.9)`),
  finer/tighter-LD than the phase-1 hapFIRE map → ~single-haplotype units for fine-mapping.
- **All-class pooled WZA screen** (SNP+indel+SV in one block-p) — nominates blocks on
  TOTAL class-agnostic evidence (screening on SV-only WZA would be circular).
- **SuSiE-RSS** within each SV-containing selected block → PIPs / 95% credible sets.
  Driver = SV is the **lead (top-PIP) of a credible set**; passenger = SV never leads.

## Pipeline (this folder)
1. `build_allclass_wza.py` — re-stamp gen9 kendall (snp/sv/smallindel) with clq0.9, pool,
   run WZA. **deg-2 verified clean** (0 NaN, max 6,131 SNPs/window, RMSE 0.54 vs deg-7 0.37)
   → dropped the deg-7/cap2000 workaround. Also emits `block_composition_*.csv`.
2. `build_finemap_inputs.py --level {site,pool}` — per-block SuSiE inputs:
   - `site` (honest primary): 31-site z (flower-wgt plot→site, N=31) + founder-panel LD.
   - `pool` (full-power sensitivity): 355-pool z + pool-AF LD (n=355, internally consistent).
3. `susie_one.R` (susieR 0.14.2, r_env; needs `LD_LIBRARY_PATH=/usr/lib64`) + `run_susie.py`
   — susie_rss per block + driver readout. `--indir … --n …`.
4. `driver_passenger_FINAL.csv` — site-vs-pool comparison.

## RESULT — SVs are PASSENGERS, robustly

Screen: 58,308 clq0.9 blocks; **0 SV-exclusive** (SNP density defeats it → every SV
co-locates with SNPs → all fine-mappable). Top-block SV enrichment is a **block-size
artifact** (SV blocks median 73 rec vs 7; size-matched enrichment 1.0–1.3×, n.s.).
71 BH q<0.05 blocks, **22 contain an SV** = fine-map worklist.

| within-block fine-map | site N=31 (honest) | pool N=355 (full power) |
|---|---|---|
| blocks with any 95% CS | 3/22 | 17/22 |
| **SV leads a credible set (driver)** | **0/22** | **0/22** |
| SV out-PIPs the best SNP | 3/22 | 0/22 |
| max SV-PIP (any block) | 0.044 | 0.165 |
| fine-mapped lead class | SNP/indel | SNP/indel |

- **SVs never lead a credible set at either honest OR maximum power.** Where an SV enters a
  CS (6/22 pool) it is a low-PIP passenger member of a SNP/indel-led set (e.g. Chr5_4221:
  smallindel lead pip 0.91, SV pip 0.005; Chr4_6307: 410-variant diffuse CS, SV pip 3e-4).
- The **site↔pool contrast** shows the per-locus climate signal only *localizes* under the
  inflated N=355 (17 CS) and mostly evaporates at the honest N=31 (3 CS, diffuse) — i.e. the
  block signal is largely a pool-within-site pseudo-replication artifact. SuSiE flagged the
  site-level founder-LD mismatch directly (3 blocks: "prior variance unreasonably large").

## Caveats / open
- bio1 only. Block SELECTION uses 355-pool kendall (inflated); the driver readout is
  within-block *relative* so it is robust to that, but the 71 "BH-sig" set is not a
  site-permutation-calibrated selection.
- Founder-LD (site) vs pool-AF-LD (pool) are two reference choices; pool-AF LD is
  self-consistent with the z's (recommended). No individual genotypes → reference-LD SuSiE.

## Files / how to re-run (merged from HANDOFF)
Env: kmate python `/global/home/users/tbellg/miniforge3/envs/kmate/bin/python`; SuSiE in `r_env`,
needs `export LD_LIBRARY_PATH=/usr/lib64` (susieR 0.14.2). Compute node only. Outputs →
`analysis/grenenet_selection/extras/driver_passenger/results/`.
1. `build_allclass_wza.py --model kendall --gen 9 --climate bio1`
   → `allclass_*.records.csv`, `wza_*_deg{2,7}.csv`, `block_composition_*.csv`,
     `finemap_worklist_*.csv` (22 SV-containing BH-sig blocks). Runs WZA via `../wza_script.py`.
2. `build_finemap_inputs.py --level site` and `--level pool`
   → `finemap_inputs/` and `finemap_inputs_pool/` (per-block npz: z, R, is_sv, cls, pos…).
   site = 31-site flower-wgt z + founder LD (`panel/arch3/chrN/var_pa_*`); pool = 355-pool z +
   pool-AF LD (consistent). Both matched 100% of variants to the founder panel.
3. `run_susie.py [--indir …_pool] [--n 355] [--out …_pool.csv]` (default site, n=31; uses `susie_one.R`)
   → `finemap_susie_results{,_pool}.csv`, then `driver_passenger_FINAL.csv`.

### Inputs this depends on (already on disk)
- Per-class gen9 kendall: `analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/kendall/kendall_{snp,sv,smallindel}_gen9_bio1.csv`
  (chrom,pos,ref_len,alt_len,MAF,block,tau,pval; row-aligned to class_matrices records/af).
- `phase1_replication/class_matrices/{cls}_gen9_af.npy` (355 pools × records) + `gen9.pools.csv`
  (355 pools, 31 sites, bio1 per site, total_flowers for weighting).
- Founder panel `panel/arch3/chr{N}/var_pa_231_arch3_chr{N}.{var_pa,var_called,meta}.npz`;
  match key = (pos,ref_len,alt_len).

### Open / next steps
- **bio2–19 sweep** (only bio1 done): re-run per-class kendall (`../phase1_replication/run_kendall.py`)
  for other bioclim, then repeat steps 1–3 — confirms the passenger call isn't bio1-specific.
- Optional: site-permutation null on block-p before calling blocks "selected" (the 71-block BH set
  is not permutation-calibrated).
- Optional: gene-annotate the SNP/indel-led credible sets (the actual climate signals) via
  `../phase1_replication/build_significant_genes.py`.

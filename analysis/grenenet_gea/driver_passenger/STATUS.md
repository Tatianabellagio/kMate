# Driver-vs-passenger: are SVs the causal/lead variant in climate-selected blocks?

2026-07-03. Fine-mapping test that beats the LD confound: within each climate-selected
LD block, is the SV the lead/causal variant beyond the best tagging SNP, or a passenger?

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

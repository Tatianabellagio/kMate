# HANDOFF — driver-vs-passenger SV fine-mapping (2026-07-03)

## The question
Within climate-selected LD blocks, is the SV the **causal/lead** variant beyond the best
tagging SNP (driver), or always ≤ the best SNP (passenger)? This is the only framing that
beats the SV–SNP LD confound. Prior work leaned "passenger"; this session did it cleanly,
per selected block, with an LD-aware method.

## What was decided (design)
- **Blocks = clq0.9** BigLD islands (`lib.assign_clq_blocks(chrom,pos,r2=0.9)`; files
  `results/grenenet_gea/blocks_mcf90/chr{N}_clq0.9_blocks_clq0.9.tsv`, 58,376 blocks,
  median ~223 bp/7 var). Finer than the phase-1 hapFIRE map used by `build_class_matrices`.
- **All-class pooled WZA** for the block screen (SNP+indel+SV → one block-p). Rationale:
  nominate blocks on TOTAL class-agnostic evidence; SV-only WZA would be circular. Screen is
  the *less biased* choice (avoids winner's-curse toward SVs).
- **deg-2 WZA** (canonical Booker) — VERIFIED clean on clq0.9 (0 NaN, max 6,131 SNPs/window,
  SD-fit RMSE 0.54 vs deg-7 0.37). The deg-7/cap2000 workaround is NO LONGER NEEDED.
- **Fine-map = SuSiE-RSS** per SV-containing selected block. Driver = SV is the **lead
  (top-PIP) of a 95% credible set**. (NB: "SV merely a member of a CS" is NOT a driver —
  it's usually a passenger member of a SNP-led set.)
- Ran BOTH: **site-level (N=31, honest primary)** and **pool-level (N=355, full-power
  sensitivity, consistent LD)**. User explicitly wanted the pool sensitivity after the
  honest run showed no resolution.

## RESULT (headline) — SVs are PASSENGERS, robust at both power levels
`results/grenenet_gea/driver_passenger/driver_passenger_FINAL.csv`

| within-block fine-map | site N=31 | pool N=355 |
|---|---|---|
| blocks with any 95% CS | 3/22 | 17/22 |
| **SV leads a credible set** | **0/22** | **0/22** |
| SV out-PIPs best SNP | 3/22 | 0/22 |
| max SV-PIP (any block) | 0.044 | 0.165 |
| lead class of fine-mapped signals | SNP/indel | SNP/indel |

Also: 0 SV-exclusive blocks even at clq0.9 (SNP density); block-level SV "enrichment" among
top WZA hits is a **block-size artifact** (SV blocks median 73 rec vs 7; size-matched 1.0–1.3×,
n.s.). Site↔pool contrast shows the per-locus signal only localizes under the inflated N=355 →
block signal is largely pool-within-site pseudo-replication (SuSiE flagged the site-level
founder-LD mismatch on 3 blocks: "prior variance unreasonably large").

## Files / how to re-run (all in `analysis/grenenet_gea/driver_passenger/`)
Env: kmate python `/global/home/users/tbellg/miniforge3/envs/kmate/bin/python`;
SuSiE in r_env, needs `export LD_LIBRARY_PATH=/usr/lib64` (susieR 0.14.2 installed this session).
Compute node only (never login). Outputs → `results/grenenet_gea/driver_passenger/`.
1. `build_allclass_wza.py --model kendall --gen 9 --climate bio1`
   → `allclass_*.records.csv`, `wza_*_deg{2,7}.csv`, `block_composition_*.csv`,
     `finemap_worklist_*.csv` (22 SV-containing BH-sig blocks). Runs WZA via `../wza_script.py`.
2. `build_finemap_inputs.py --level site`  and  `--level pool`
   → `finemap_inputs/` and `finemap_inputs_pool/` (per-block npz: z, R, is_sv, cls, pos…).
   site = 31-site flower-wgt z + founder LD (`panel/arch3/chrN/var_pa_*`); pool = 355-pool z +
   pool-AF LD (consistent). Both matched 100% of variants to the founder panel.
3. `run_susie.py [--indir …_pool] [--n 355] [--out …_pool.csv]` (default site, n=31)
   → `finemap_susie_results{,_pool}.csv`. Uses `susie_one.R`.

## Inputs this depends on (already on disk)
- Per-class gen9 kendall: `results/grenenet_gea/phase1_replication/kendall/kendall_{snp,sv,smallindel}_gen9_bio1.csv`
  (cols chrom,pos,ref_len,alt_len,MAF,block,tau,pval). Row-aligned to class_matrices records/af.
- `phase1_replication/class_matrices/{cls}_gen9_af.npy` (355 pools × records) + `gen9.pools.csv`
  (355 pools, 31 sites, bio1 constant per site, total_flowers for weighting).
- Founder panel `panel/arch3/chr{N}/var_pa_231_arch3_chr{N}.{var_pa,var_called,meta}.npz`
  (231×R sparse alt-presence; meta has pos,ref,alt,ref_len,alt_len). Match key = (pos,ref_len,alt_len).

## OPEN / next steps
- **bio2–19 sweep** (only bio1 done): re-run kendall per-class for other bioclim, then repeat
  steps 1–3. Confirms the passenger call isn't bio1-specific. (kendall = `phase1_replication/run_kendall.py`.)
- Block SELECTION still uses 355-pool kendall (inflated); driver readout is within-block
  *relative* so robust, but the 71-block BH set is not site-permutation-calibrated. Optional:
  site-permutation null on the block-p before calling blocks "selected".
- Effective-N binomial screen (user's separate track: weight to 31-site level, p<0.05 84%→48%)
  not wired into this folder; the fine-map already uses site-level Pearson z, which is the same
  effective-N idea applied to the fine-map statistic.
- Could add: gene annotation of the SNP/indel-led credible sets (the actual climate signals),
  reusing `phase1_replication/build_significant_genes.py` machinery.
```

# ★ FINAL UNIT MAP + FINAL BENCHMARK (2026-06-19) — read this first

**THE units to benchmark** (k-mer-covered, LD-grown; = kMate h-window AND selection unit):
- genome-wide: `results/grenenet_gea/blocks_mcf90/final_units_dynld_K500.tsv`
- per-chrom (use for Chr1 benchmark): `results/grenenet_gea/blocks_mcf90/chr{N}_units_dynld_K500.tsv`
- cols: `chrom  start_pos  end_pos  n_variants  panel_kmers  covered`. **22,939 units; 72% COVERED**
  (panel_kmers>=500 -> local-fit; 16,403 units), 28% desert (global fallback). median 32 var /
  2.1 kb / 926 panel k-mers. Built by `dynamic_ld_blocks.py` (grow CLQ0.9 blocks along the LD
  gradient until >=500 panel k-mers; HapFM untouched).
- **CODE-AUDITED + BUG-FIXED + VALIDATED (2026-06-19):** 3-agent adversarial review. k-mer
  tagging byte-identical to `block_em.assign_kmers_to_blocks`; merge is version-stamped
  (highest-current-LD-first) and uses gap-inclusive interval k-mer counts (`km_of`, validated
  byte-identical to block_em per-block counts). This is the version to benchmark — earlier
  positional/kmeradaptive/ve07 maps were DELETED (superseded; don't look for them).
- characterization notebook: `analysis/grenenet_gea/notebooks/blocks_units_decision.ipynb`.

**HOW TO RUN the final benchmark** (window-kMate on the unit map, score vs truth):
- runner pattern: `benchmarks/ldblock_window_test/run_ldblock_window.sbatch <POOL> p231 <BLOCKS_TSV> 50`
  (calls `python -m kmate.per_sample_per_chrom --block-mode window --blocks-tsv <map>
  --min-kmers-per-block 50 ...`; Chr1). Pass BLOCKS_TSV = chr1_units_dynld_K500.tsv.
- multi-seed array example: `benchmarks/ldblock_window_test/run_coarse_sweep.sbatch` (6 seeds).
- **scorer (VALIDATED, reproduces base R²≈0.958):** `benchmarks/ldblock_window_test/score_ldblk.py`
  — row-aligned to `recomb_truth_RAW.tsv.gz` (NOT atomized), drops nonfinite + (truth>0|est>0).
  `--validate <out.tsv> <truth.gz> <hblocks.npz>` or `--sweep`. %local from h_blocks Chr1_status.
- sims: `benchmarks/p231/sims/cov10_n50_g3_s{42..47}_self97_hotspots_dom500_p231_chr1/` (6 seeds,
  selfing/low-recomb) AND `..._g3_s42_hotspots_p231_chr1` (outcross, HIGHER recomb).
- panel: kmer_pa `data/kmer_pa_231_arch3_filt2inv/kmer_pa`, var_pa `panel/arch3/chr1/var_pa_231_arch3_chr1.*`.

**WHAT THE FINAL BENCHMARK SHOULD SHOW:** (1) on the dynld unit map, ~69% local-fit + accuracy
holds (R²≈0.99); (2) KEY open test — run the **higher-recombination (outcross hotspots) arm**
and show local-fit accuracy BEATS pure-global there (the insurance the selfing sims couldn't
demonstrate — see memory gea-blockwindow-benchmark-result). Coarse-floor sweep results already
on disk: `benchmarks/ldblock_window_test/coarse_sweep_scores.csv`.

---

# LD-block / haploblock definitions — handoff for benchmarking

2026-06-19. These blocks are the candidate **windows for window-based kMate** (replacing the
old fixed 10 kb windows). Goal for the benchmarking agent: run the p231/p80 sims through
**block-based** window-kMate and measure accuracy on these windows.

## THE BLOCKS (what to benchmark on)

Final locked map — **CLQcut 0.9, panel-support filter `min-called-frac 0.9`** (variants kept
only if genotyped in ≥208/231 founders):

```
results/grenenet_gea/blocks_mcf90/chr{1..5}_clq0.9_blocks_clq0.9.tsv
```
- TSV cols: `chrom  start_pos  end_pos  n_variants` (1 row per block; coords are TAIR10 1-based).
- **58,376 blocks** genome-wide (Chr1 15789, Chr2 8512, Chr3 10674, Chr4 9339, Chr5 14062).
- Size: **median 223 bp / 7 variants**, mean ~1 kb (right-skewed, tail to 328 kb). 98% < 10 kb.
  Cover ~51% of the 119 Mb genome (variant-dense regions; sparse/repeat regions are gaps).
- vs 10 kb windows: there are 11,916 of those → these blocks are **4.9× more numerous and far
  thinner**. ⚠️ Many 2–7-variant blocks are likely under-powered for h-estimation — consider a
  **min-variant floor / merge tiny adjacent blocks** (or fall back to a fixed window where blocks
  are too thin) before the run. Open design decision.

## HOW THEY WERE BUILT (faithful HapFM partition on our all-class panel)

- Definition code: `analysis/grenenet_gea/recompute_blocks.py`
  (coarse `CompleteLDPartition` corr=0.2 + fine `BigLD` gpart density, copied verbatim from
  `/global/scratch/users/tbellg/HapFM/bin/` — `block_partition.py` + `BigLD.R`).
- Driver (the exact run that produced blocks_mcf90): `analysis/grenenet_gea/blocks_recompute_mcf90.sbatch`
  (`--corr 0.2 --window 50 --maf 0.05 --clqcut 0.9 --min-called-frac 0.9`, array=chrom).
- Panel input (founder genotypes, all-class SNP+indel+SV): `panel/arch3/chr{N}/var_pa_231_arch3_chr{N}.{var_pa,var_called,meta}.npz`
  (var_pa = 231×R 0/1 alt-presence; var_called = 231×R call mask; meta['pos'] = positions).
  Blocks are defined on MAF>0.05, call-rate≥0.9, unique-position records; founder missing imputed
  to per-variant mean for LD.

## HAPLOTYPE UNITS WITHIN BLOCKS (if benchmarking the unit layer too)

- Each n_eff≤2 block = 1 unit (block AF); each n_eff>2 block is split into HapFM-xmeans
  haplotype clusters (gate n_eff>2), keeping clusters with PC1-VE≥0.7.
- Clustering code: `analysis/grenenet_gea/block_cluster_pc1ve.py` (xmeans; **needs the numpy.warnings
  shim** — see below) and `block_unit_frontier.py` (per-block + per-cluster PC1-VE).
- Final unit registry: `results/grenenet_gea/blocks_mcf90/final_units_ve07.csv`
  (66,032 units; cols unit_id, chrom, block_start, block_end, n_eff, unit_type[block|hapcluster],
  cluster, cluster_freq, n_sig, pc1_ve). NB: does NOT yet store per-unit variant membership
  (cluster ids are from a non-seeded xmeans run).

## ENVIRONMENTS & GOTCHAS

- Python: `/global/home/users/tbellg/miniforge3/envs/kmate/bin/python`
- BigLD R: `/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript`; needs
  `export LD_LIBRARY_PATH=/usr/lib64` (system libxml2 for igraph). gpart::BigLD was installed
  from stripped source (see memory `gpart-bigld-install`).
- `pyclustering` (xmeans) crashes on modern numpy — shim `np.warnings=warnings` BEFORE import
  (already in `block_cluster_pc1ve.py`; see memory `pyclustering-numpy-warnings-shim`).
- Run heavy jobs via **sbatch** (template: `blocks_recompute_mcf90.sbatch`), not background.

## SIMS FOR BENCHMARKING

`benchmarks/p231/` and `benchmarks/p80/` (built sim outputs under `benchmarks/{p80,p231}/sims/`,
gitignored; truth = recomb_truth_raw + _atomized). See memory `sim-outputs-under-benchmarks`.
Window-mode kMate previously needed `threadpoolctl` (memory `threadpoolctl-window-blocker`).

## COARSENESS BENCHMARK (2026-06-19 — the live task)

Benchmark finding: on the base CLQ0.9 map only **~42% of blocks get a local h fit** (≥50
observed k-mers at 10×); ~30% coverage-limited (fallback→global), ~28% panel-limited (zero
k-mers). Median-7-variant blocks are TOO THIN. Looser CLQcut doesn't fix it (CLQ0.5 still
median 10 var, 43% ≥14-var). So we must MERGE to a size floor and benchmark "how coarse".

Coarsened candidate window maps (greedy-merge adjacent blocks to a min-variant floor;
`make_coarse_blocks.py`): `results/grenenet_gea/blocks_mcf90/coarse/floor{8,15,25,40}.tsv`
(+ per-chrom `chr{N}_floor{F}.tsv`). Cols: chrom start_pos end_pos n_variants n_merged.
Ladder: floor8=34404 win/med17var/1kb, floor15=25671/26var/1.8kb, floor25=19573/38var/2.9kb,
floor40=14766/55var/4.4kb (base=58376/7var/223bp).

BENCHMARK: run window-kMate on each floor map through p231/p80 sims; per floor report
(1) % windows locally fit, (2) accuracy on fit windows (MAE/R²/r), (3) overall w/ fallback.
Ideally also sweep a min-OBSERVED-KMER floor (you have the per-block k-mer counts) — variant
floors are the proxy starting point. Sweet spot = smallest floor where local-fit saturates +
accuracy plateaus. COHERENCE cost (median PC1-VE per window per floor) computed separately on
the GEA side (merging across LD boundaries lowers PC1-VE) — combine the two axes to pick.

## CONTEXT

Why blocks not 10 kb: blocks are LD-coherent (≈single-haplotype) so per-block h estimates a few
(median n_eff 2.46) haplotype freqs, vs a 10 kb window forcing 231 founder freqs over a mosaic.
Different estimand → the old 10 kb benchmark is NOT a clean accuracy bound; must re-benchmark.
Full session rationale in memory `gea-block-haplotype-units`.

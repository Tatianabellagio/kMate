# Session handoff — hapFIRE-projection for SV allele frequencies

**Project root:** `/home/tbellagio/scratch/hapfire_sv/`
**Last updated:** 2026-05-07

## TL;DR (start here)

The pipeline pivoted long ago from "hapFIRE on SNPs + project to SVs" to **`cactus_em` end-to-end**: a per-sample k-mer Poisson EM on the 231-founder simplex, projected through `cn_var_231_v2` to per-record AF for SNPs + INS + DEL in one pass.

**Production recipe (locked in 2026-05-07):**
```bash
python poolfreq/src/per_sample_per_chrom.py \
    --cn-kmer-prefix poolfreq/data/cn_full_231_v2/cn \
    --cn-var       poolfreq/data/cn_var_231_v2.cn_var.npz \
    --cn-var-meta  poolfreq/data/cn_var_231_v2.meta.npz \
    --reads <r1.fq> <r2.fq> --sample <name> --out <out.tsv> \
    --threads 8 --chroms Chr1 Chr2 Chr3 Chr4 Chr5 \
    --block-mode window --window-bp 10000 \
    --global-anchor-weight 0.3 \
    --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5
```

Wall: ~25 min on 1×8-core, ~40 GB peak per Chr1 cov50 sample. ~3× faster than xwu's hapFIRE end-to-end.

**Cross-regime SNP R² @ cov50× (poly Chr1)**:

| regime | winner ★ | global | hapfire_coarse | window_200kb |
|---|---|---|---|---|
| n200_g1 (easy) | **0.9960** | 0.9967 | 0.9946 | 0.9903 |
| n50_g1 (med) | **0.9888** | 0.9867 | 0.9827 | 0.9830 |
| n50_g3 (hard) | **0.9664** | 0.9489 | 0.9418 | 0.9726 |
| n50_g3_skewed (stress) | **0.8117** | 0.8305 | 0.7715 | 0.6557 |

Coverage-stable: cov50→cov10 SNP R² Δ ≤ 0.0014 (see `FINAL_RESULTS_cov10.ipynb`).

**Authoritative result notebooks (READ FIRST):**
- `FINAL_RESULTS.ipynb` — head-to-head 5 methods, 4 regimes, cov50×
- `FINAL_RESULTS_cov10.ipynb` — same at cov10×, plus coverage-robustness deltas
- `RECOMB_SWEEP_RESULTS_slim.ipynb` — full method × LD-tier × var_type sweep
- `METHODS_TRIED.md` — design-space map: every method/route we explored, status, and the published work each piece relates to (kallisto, Mora, KAGE 2, StrainFacts, microhaplotype pool, etc.)
- `RESULTS_LOG.md` 2026-05-06 / 2026-05-07 entries — sweep-by-sweep reasoning (Routes 1/2/3, λ-sweep, smoothing-sweep, kallisto-EM postmortem)

## What's pending

1. **Production scale-out on ~2,500 evolved GrENE-Net samples.** SLURM template already validated at `poolfreq/tests/run_site_array_perchrom.sh`. Estimated wall: ~1.5–5 days at cluster-wide concurrency. Per-chrom driver fits in 64 GB on memex (no longer needs the 200 GB bse-2021-001 bottleneck).
2. **bigld_haplotype mode** (in-progress as of 2026-05-03 evening) — full Chr1 cn build + 3-regime re-runs pending. Tests whether sequence-unique-per-block dedup beats fixed-window EM. **Not blocking** production scale-out; can proceed in parallel.
3. **Real-data validation against window_200kb** on the existing 57 site04 evolved BAMs. The 200 kb baseline is already computed; rerun the 10 kb anchor+smooth recipe on the same set and compare.
4. **Recomb-aware HMM transition tuning** (4 cM/Mb AT recomb rate) — current smoother uses a flat α; principled λ may close the remaining 0.6 pp on n50_g3.
5. **Empirically test the panel-EM coverage-saturation claim.** We currently *assert* that the per-k-mer Poisson EM saturates above ~5× because of the 16/32 k-mers-per-allele cap × 16M panel k-mers — and the cov10 vs cov50 deltas (≤0.0014 SNP R²) are consistent with that — but we have not actually run the curve. Plan: re-use `setup_cov10_sim.sh`'s structure to generate cov∈{1, 2, 3, 5, 10, 20, 50}× reads on the n50_g3 regime (most stress-testing of the recipe), run the 5 methods, plot R² and per-record std vs coverage. Expected: panel-EM curves flat above ~5×; freqk monotone in coverage; HARP-window between. Until we run this, the "above ~5×" claim in the cov10 wrap cell is a prediction from `MODEL_SPEC.md` arithmetic, not a measurement.

## Document map

| File | Purpose |
|---|---|
| `HANDOFF.md` (this file) | What to do next; chronologically newest entries on top |
| `FINAL_RESULTS.ipynb` | Authoritative cov50 head-to-head; production winner table |
| `FINAL_RESULTS_cov10.ipynb` | Same at cov10; coverage-robustness deltas |
| `RECOMB_SWEEP_RESULTS_slim.ipynb` | Full method × LD-tier × var_type table |
| `METHODS_TRIED.md` | Design-space map: methods tried (★/✓/◯/✗/⋯) + literature analogs |
| `RESULTS_LOG.md` | Chronological log of validation/benchmark numbers |
| `MODEL_SPEC.md` | Math: generative Poisson model, EM, identifiability |
| `ALGORITHM.md` | Implementation walkthrough of the EM solver |
| `CACTUS_EM_MATH.md` | Detailed math for the cactus_em window/global/anchor |
| `BACKGROUND.md` | Why we're doing this (project framing) |
| `SUMMARY.md` | **Historical** — original hapFIRE-projection writeup (preserved for context) |
| `CACTUS_PANGENOME_AUDIT.md` | cactus-vs-syri reconciliation |

---

## 2026-05-03 evening — bigld_haplotype mode, in-progress

New `--block-mode bigld_haplotype` for cactus_em. Unifies: BigLD-fine-block
resolution + sequence-unique-per-block dedup of founders + EM on the
class simplex with full-haplotype k-mer evidence.

**Why**: under heavy recombination (cov50_n50_g3 sim), neither current
bigld method beats `window_200kb` (R² 0.838/0.882 vs 0.946). bigld_haplotype
should close that gap by giving fine-block EM ~1000× more discriminative
k-mers than the current bigld_panel cn_full has.

**Done**:
- Builder: `poolfreq/scripts/build_block_haplotype_cn.py`
- EM core + per-sample driver: `poolfreq/src/block_haplotype_em.py`,
  `poolfreq/src/per_sample_bigld_haplotype.py`
- SLURM: `sims/visor_freqk/scripts/run_cem_bigld_haplotype.sh`
- Compare script wired in: `compare_recomb_stratified.py` includes
  `bigld_haplotype` in its methods set

**Pending** (in order):
1. Build full-Chr1 `block_haplotype_cn/chr1_full.npz` — long-running
   (~30-60 min wall; per-founder FASTA hashing per block is the
   bottleneck). Launch detached via `nohup setsid` so it survives shell
   teardown.
2. End-to-end smoke run on one recomb sim (e.g. cov50_n50_g3) to
   sanity-check TSV output.
3. Re-run all three Chr1 sims (n=50_g=1, n=50_g=3, n=200_g=1) with
   bigld_haplotype mode and re-stratify.

**Design notes** (full design + smoke metrics in
`memory/project_bigld_haplotype_design.md` and the 2026-05-03 evening
RESULTS_LOG entry):
- BigLD block index used **only** for `(chrom, start, end)` positions.
- Per-block dedup is by **full FASTA k-mer set** — captures SV variation
  among same-SNP-haplotype founders. Drops the SNP-only `hap_idx_d1`
  approach the first prototype used.
- Cross-block dedup drops k-mers shared across blocks (~35% loss in
  smoke test on 50 blocks; clean per-block EM independence).

---

---

## 2026-05-02 — per-chrom window-mode is the memory + speed win

`poolfreq/src/per_sample_per_chrom.py --block-mode window` is now validated as a **strict drop-in replacement** for the genome-wide window-mode driver. Tested on MLFH040120180306 (a real site04 evolved sample) against the existing genome-wide window-mode result:

| | Genome-wide window | Per-chrom window |
|---|---|---|
| Wall time | 1:08:34 | **40:12 (-41%)** |
| Peak RSS | 137 GB | **39 GB (-72%)** |
| Output AF vs reference | — | **bit-identical** (max abs diff 1e-5, R²=1.000000) |
| Cluster nodes that fit | bse-2021-001 only (200 GB) | any memex 64+ GB (15+ available) |

Why bit-identical: `block_em.define_windows` and `project_blocks_to_records` are already strictly per-chrom (within-chrom interpolation only), so loading one chrom at a time only changes *when* the data is in RAM, not the math. Validation script: `poolfreq/tests/compare_perchrom_window_vs_gw.py`.

Implication for production scale-out: **drop the 200 GB SLURM allocation in favor of 64 GB on memex**. The genome-wide concern from the 2026-05-01 entry (memory bottleneck) is resolved. Production SLURM template for evolved samples: `poolfreq/tests/run_site_array_perchrom.sh` (parametrized via env vars for MANIFEST/OUT_DIR/BLOCK_MODE).

Updated production scale-out estimate (~2,500 samples, window-mode):

| Setup | Per-sample wall | Concurrent | Total wall |
|---|---|---|---|
| Per-chrom on memex (8 nodes × 4 tasks) | 40-60 min | ~30 | **~3-5 days** |
| Push concurrency (cluster-wide, 100 slots) | 40-60 min | 100 | **~1-1.5 days** |

---

## 2026-05-02 — k-mer counting alternatives benchmarked, no meaningful win

Tested KMC 3 (drop-in for jellyfish) and rapidgzip (parallel gzip front-end) on SEEDMIX_S1 reads + Chr1 query (20.8M k-mers, 8 cores). Goal: reduce the k-mer counting stage which is ~50-75% of per-sample wall.

| Tool | Wall | Peak RSS | vs jellyfish | Output identity |
|---|---|---|---|---|
| zcat + jellyfish (current) | 9:10 | 23 GB | 1.00× | reference |
| KMC 3 (count + intersect) | 7:57 | 41 GB | 1.10× | 99.996% (saturates at 255 without `-cs10000`) |
| `cat \| rapidgzip 8t \| jellyfish 8t` | 8:34 | 24 GB | 1.07× | bit-identical |
| rapidgzip R1; rapidgzip R2 \| jellyfish 8t | 8:29 | 24 GB | 1.08× | bit-identical |

KMC's 80% memory regression and rapidgzip's marginal speedup mean **neither is worth integrating** at this scale. Earlier "1.45-1.88× rapidgzip speedup" was a CLI bug (`rapidgzip -d -c R1 R2` only processed one file). Once corrected to use `cat | rapidgzip` or sequential rapidgzip calls, the speedup collapses to ~7%.

Why so little gain: jellyfish at 8 threads is already saturating CPU at 618% with zcat as front-end; the gzip decompression bottleneck only bites at lower thread counts. The next real lever (per the literature) is **SSHash** (specialized k-mer dictionary, ~10× lookup speedup) or **LP-MPHF + custom Rust counter** (3-10×). Both are ~1-2 days of integration work and can be deferred until per-chrom window-mode + memex concurrency proves insufficient at full GrENE-Net scale.

Bench scripts: `poolfreq/bench_kmc/run_bench.sh` (KMC), `run_rapidgzip_bench2.sh`.

---

## 2026-05-02 — recombinant pool simulation built, first result counter-intuitive

`sims/visor_freqk/pool_sweep_82_recomb` closes HANDOFF.md item 1 (recombination simulation). Pipeline: build N mosaic haploids via Poisson crossover at 4 cM/Mb (`make_recomb_mosaics.py`), pool reads via VISOR SHORtS, run cactus_em window + global, derive truth from ancestry tracks. See `sims/visor_freqk/RECOMB_SIM.md`.

**First sim (cov 10×, N=50, G=1, seed=42)** — surprising result:

| Mode | n_polymorphic | R² | MAE | RMSE |
|---|---|---|---|---|
| cactus_em window | 3.51M | 0.934 | 0.030 | 0.057 |
| cactus_em global | 3.51M | **0.986** | **0.018** | **0.025** |

**Global beats window by 5pp R² across SNP / INS / DEL.** Interpretation: at N=50 with only 1 generation of recombination, the pool is globally well-mixed — per-window founder mixtures barely differ from the genome-wide mean, so window-mode just adds per-window EM noise without buying local structure. Window-mode is expected to win in regimes where ancestry blocks are sharper (G≥3, smaller N).

This contradicts our prior assumption that window-mode is always better for evolved samples. The site04 result (R²=0.97 vs hapFIRE) was on real ~3-gen evolved samples and used window-mode without a global comparison; we don't actually know if window or global is better there. **Next**: G=3 simulation, and a cheap global-vs-window comparison on one site04 sample to see which regime applies in production.

---

## 2026-05-03 — runtime benchmarks: cactus_em vs hapFIRE

Recorded from actual SLURM jobs on the 231-eco recomb sim and on real site04 BAMs:

| Method | Wall per sample | Peak RSS | Notes |
|---|---|---|---|
| **cactus_em** (per_sample_per_chrom.py, window + global, 231-cn, sim) | **~1:13** | 47 GB | jobs 59791, 59792 |
| **cactus_em** (per_sample_per_chrom.py, window-mode only, 231-cn, real site04 BAMs) | **~1:00-1:17 (median 1:08)** | 130-140 GB | site04 array 58866, n=57 |
| **hapFIRE** (full pipeline, `-haplotype True`, custom partition) | **~3:00-3:10** | 27-29 GB | sim jobs 59631, 59632 |
| BAM align (bwa+dedup+filter) | 0:08-0:15 | 13 GB | one-off per sample |

**Headline: cactus_em is ~3× faster than hapFIRE end-to-end.** For the production target of ~2,500 GrENE-Net samples, this is days vs weeks of cluster time at the same concurrency.

**Why hapFIRE is slow**: Phase 2 (HARP `like` + `freq`) runs serially across the ~13 coarse independent LD blocks in a Python for-loop. ~20 min per HARP call × 13 = ~4h serial. No internal parallelism. To bring hapFIRE down to ~40 min/sample we'd need to refactor `haplotype_frequency_estimation` to `multiprocessing.Pool` across chroms or blocks.

**For production scale-out** (~2,500 samples):

| Method | Per-sample | C=10 | C=30 |
|---|---|---|---|
| cactus_em window | 70 min | ~12 days | **~4 days** |
| hapFIRE full | 3 h | ~31 days | ~10 days |

**Production-relevant takeaway**: even if SV accuracy comparison is a wash between the two methods, cactus_em wins on **runtime** by 3-4×. This is a real production lever independent of accuracy.

---

## 2026-05-02 — measured per-sample compute + production scale-out estimate

After running cactus_em on all **57 site04 evolved samples** (window-mode, 231-cn v2 genome-wide driver), we have hard numbers:

**Per-sample compute (genome-wide window-mode driver, 4 CPUs/task, 200 GB SLURM alloc):**

| Metric | Value (n=57) |
|---|---|
| Wall time per sample (min/median/mean/max) | **53 / 62 / 63 / 77 min** |
| Peak RSS (MaxRSS) | **130-134 GB** (median 130.6) |
| CPU efficiency (AveCPU/Elapsed) | ~3× / 4 cores ≈ **75%** |
| Total compute for 57 samples | **~60 wall-h, ~240 core-h** |
| Cluster concurrency observed | **7-12 tasks** at once on bse-2021-001 (only ≥200 GB node) |

**Per-chrom driver alternative (`poolfreq/src/per_sample_per_chrom.py`):** validated R²=0.998 vs genome-wide on SEEDMIX_S1 (MAE 0.006). Tradeoff: peak RSS drops to ~44 GB (3× reduction) at the cost of ~1.5× wall time (k-mers re-counted per chrom). Fits on the 64+ GB memex nodes, unlocking ~3× more cluster concurrency.

**Production scale estimate (~2,500 GrENE-Net samples):**

| Setup | Per-sample wall | Concurrent slots | Total wall-clock |
|---|---|---|---|
| Genome-wide on bse-2021-001 only | 63 min | ~10 | **~10-15 days** |
| Per-chrom on memex (8 nodes × 4 tasks) | 95 min | ~30 | **~5 days** |
| Push concurrency (cluster-wide, 100 slots) | 63-95 min | 100 | **~1.5-2 days** |

Per-chrom mode is per-sample slower but **net-faster at scale** because it removes the 200 GB-node bottleneck. Recommended path for scale-out: switch SLURM template to `per_sample_per_chrom.py` + 64 GB allocation + memex partition.

---

## 2026-05-01 — cactus_em memory bottleneck (production scale-out concern)

**Current per-sample memory:** 100-150 GB peak. **96 GB SLURM allocation OOMs** (verified 2026-05-01 on site04 array — every task killed at OOM). **Production-safe SLURM spec: ≥200 GB per task.**

**Where the memory goes (231-cn pipeline):**
- `cn_full_231_v2` (founder × k-mer sparse matrix, 231 × ~80 M k-mers): ~30-50 GB held in RAM during EM
- Jellyfish k-mer hash during counting: ~10-30 GB
- EM solver float32 workspace: ~5-10 GB (was 26 GB pre-fix; see RESULTS_LOG entry "EM solver bottleneck fix")
- `cn_var_231_v2` projection matrix: ~3-5 GB

**Comparison.** freqk uses 16-32 GB (per-bubble streaming, no joint matrix). cactus_em is 3-5× heavier on memory and 2-3× longer per-sample runtime — the cost of joint k-mer EM that buys R²=0.996 vs hapFIRE on real SEEDMIX SNPs (vs freqk's 0.94). hapFIRE itself also needs ~96 GB on full Chr1 with 231 ecotypes (HARP stage), so it's not unique to us.

**Production scale.** 2,414 evolved samples × 200 GB × ~30 min ≈ 1,200 SLURM CPU-hours and ~240 K GB-hours. Concurrency limited by available 200+ GB nodes on the partition (typically 12-24 concurrent on memex). End-to-end wall: ~10-15 hours at 12-way concurrency.

**Mitigations (not implemented, ranked by ease):**
1. **Memory-map cn_full from disk** — sparse `.npz` can be `mmap`-loaded so the OS page cache shares it across array-task workers on the same node. Should cut per-task RAM from 50 GB → 10-20 GB at modest (~1.5×) speed cost. Easy.
2. **Subsample k-mers** — coverage is saturated below 5×, suggesting we don't need all 80 M k-mers. Keep the top ~10 M most-discriminating per founder pair. Cuts cn_full ~8×.
3. **Per-chromosome EM** — solve the simplex per chrom and average. Reduces peak RAM but introduces per-chrom convergence noise.
4. **PanGenie genotyping replacement** — bypass the EM by using PanGenie's per-individual genotyper (the other team's pipeline). Different method, separate validation needed; see `pangenie_genotyping/`.

**When to fix:** if we ever scale to 10× more samples, OR if the cluster's memex partition becomes saturated for these jobs. Until then, 96 GB array-tasks at 24-way concurrency clear the 2.4 K production set in ~50-100 wall-hours.

---

## 2026-04-29 — cactus_nocap rebuild finished, was a no-op

The 42-hour rebuild of cactus without `--maxLen` (job 56169) completed.
Compared to the original capped pang VCF:

| Metric | capped | nocap | Δ |
|---|---|---|---|
| vcfbub-filtered records | 4,450,411 | 4,450,117 | -294 (0.007%) |
| max SV size | 299,982 bp | 299,982 bp | 0 |
| SVs ≥1Mb | 0 | 0 | 0 |
| SVs 100kb-1Mb | 20 | 20 | 0 |

**The maxLen flag never gated any real SV in this panel.** It controls the
maximum length of input contigs used during alignment chaining — not the
output SV size. The capped run already produced SVs up to 300 kb, far above
the maxLen=10000 setting. The 0.007% record-count delta is bidirectional
across size bins and is consistent with cactus's documented run-to-run
stochasticity (minigraph chaining tiebreakers + multi-threaded job ordering).

**Action**: keep using the original capped pang at
`/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/`. No need to
rebuild downstream artifacts.

---

## 2026-04-29 — SEEDMIX recipe-truth headline numbers

Per-record alt_freq predictions vs recipe-projected truth, 8 SEEDMIX replicates:

| Panel | Mode | R² (8-rep mean) | Pearson r |
|---|---|---|---|
| 82-founder | global | 0.994 (S1 only) | 0.997 |
| 82-founder | window | 0.97 | 0.98 |
| 231-founder | window raw | 0.685 | 0.985 |
| **231-founder | window CALIBRATED** | **0.97** | 0.985 |

The 231-founder pipeline has a panel-intrinsic 1.43× scale bias from
imputation noise putting too much mass on the 151 imputed founders. A
single linear calibration (`predicted = 1.43 × truth + 0.003`) recovers
nearly all the R². Slope is rock-stable across replicates: **1.4325 ± 0.0055**.

**Production recipe**:
1. Compute calibration once from a SEEDMIX rep (saved at `poolfreq/data/calibration_231_seedmix_s1.json`)
2. Apply to all 2,414 evolved-sample TSVs via `poolfreq/src/calibrate_alt_freqs.py apply`

**For evolved samples**: use `--block-mode window` (captures mosaic ancestry)
with this calibration. Expected per-record R² ~ 0.97 against per-record truth.

Plots in `plots/` (generated by `notebooks/run_plots.py`).

---

---

## Major direction change (2026-04-27 afternoon)

The pipeline is pivoting from **freqk on syri panel** to **PanGenie on cactus pangenome** for the SV evidence step. Reasons:
- freqk drops 88% of panel SVs at the index step (catastrophic on real data); systematically over-calls on the 12% it keeps (off-target k-mer hits).
- The cactus pangenome (already built at `/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/`) covers ~93% of syri SVs within 5kb (after correcting an earlier metric error — see `CACTUS_PANGENOME_AUDIT.md`).
- PanGenie's graph-aware k-mer indexing retains 100% of cactus bubbles (vs freqk's 12% on syri).
- Use `raw.vcf` + `vcfbub -l 0 -r 100000` — recovers 19 percentage points of syri SVs vs the cactus pre-filtered output.

Pending jobs: 56169 (cactus rebuild without `--maxLen` cap, ~40h), 56172 (vcfbub on raw without -r filter), 56173 (PanGenie-index on the better input).

---

## Implementation progress (2026-04-27 evening)

`MODEL_SPEC.md` documents the math. Code lives under `poolfreq/`.

**Done:**
- `poolfreq/src/block_solver.py` — per-block CVXPY solver (WLS + IRLS) with contamination term ω. All 7 synthetic tests pass.
- `poolfreq/src/kmer_count.py` — Jellyfish-based k-mer counter for pool-seq BAMs. Tested on visor_freqk sim BAM.
- `poolfreq/src/build_kmer_cn.py` — builds founder × k-mer copy-number matrix from PanGenie kmers.tsv.gz + diploidized VCF + reference. Handles merged bubbles via per-founder haplotype reconstruction.
- `poolfreq/tests/test_end_to_end.py` — 200-bubble Chr1 region: real cn (82 × 18,523) + synthetic Poisson counts. IRLS recovers ||h - h_true|| = 0.005, all 5 carriers identified.

**Pending implementation:**
- Encoding matrix `enc[v, j, s]` for cactus-alt → syri-SV decomposition (Step 2 in MODEL_SPEC.md). ~2-3 days.
- Real BAM end-to-end: count k-mers from a pool-seq BAM, solve, validate vs known truth.
- Apply to SEEDMIX, then to 2,415 evolved samples.

**Pending jobs:**
- 56169 cactus rebuild without `--maxLen` (~40h, completion check). Verifies no other LV=0 records were dropped beyond what raw.vcf already gives us.

## Job results from 2026-04-27 afternoon

- **56172** (vcfbub on raw without -r): syri SV containment 77.2% (vs 71.2% with -r 100000). Going off the -r filter recovers a few percent more, but bubbles >100kb (max 3.47Mb!) become problematic for PanGenie. Stay with -r 100000.
- **56173** (PanGenie-index on vcfbub-from-raw): 796K bubbles with unique k-mers genome-wide (Chr1 206K, Chr2 127K, Chr3 159K, Chr4 118K, Chr5 185K). All chromosomes indexed cleanly.
- **56178** (PCA pang_69 vs GrENE-Net 231 on Chr1 SNPs): all 69 pang accessions within GrENE-Net P99 of centroid, but pang_69 median dist (137) is 2.3× GrENE median (59). Most distant: Col-0 (already in both panels), Qar-8a, Cant-1, Sah-0, Are-X (Asian/Caucasian/Argentine relicts). Verdict: most pang_69 cluster near GrENE; adding them to the rebuild is net-positive but the most-distant ~10 (Are-X, Sah-0, Cant-1) could be excluded to avoid imputation noise.

---

## TL;DR — what we showed

We adapted hapFIRE for SVs by post-processing: run hapFIRE on SNPs unmodified, then project its 231-vector ecotype-frequency output onto SVs via a precomputed founder × SV genotype matrix:

```
f_SV(v) = Σ_e G_SV(e, v) · f_ecotype(e)
```

Identical to how hapFIRE projects to per-SNP frequencies internally. ~50 lines of Python.

**Simulation result (ground truth available, 25 hapFIRE jobs across 4 reps × 3 coverages × 5 freqs):**

| Coverage | hapFIRE-proj MAE | freqk MAE | fold-improvement |
|---|---|---|---|
| 10× | 0.0013 | 0.0997 | **77×** |
| 20× | 0.0010 | 0.0756 | **76×** |
| 50× | 0.0006 | 0.0583 | **97×** |

The cov10× point is the headline: hapFIRE stays accurate at low coverage at the SV locus because it leans on ~800 K Chr1 SNPs; freqk's k-mer counts are starved.

**Real-data result (8 SEEDMIX replicates):** hapFIRE-proj vs recipe-based panel-truth r=0.99 across all 8 samples; freqk vs same r=0.46. With the caveat that the r=0.99 is partly tautological (both quantities are dominated by AC carrier-count) — see SUMMARY.md.

---

## Limitations we already know about (in priority order)

1. **No-recombination assumption** is hapFIRE's biggest weakness for our use case. The simulation does not test it. The seedmix is the founder pool itself (still no recombination). Evolved GrENE-Net samples have 1-3 generations of recombination and selection → expect degraded per-block accuracy. **This is the single most important thing to test before recommending hapFIRE-proj as the production method for evolved samples.**

2. **Panel ascertainment ceiling.** Our SV panel covers 80 of 231 GrENE-Net ecotypes (≈33% of seedmix mass, 17.7% by recipe). SVs polymorphic only among the un-genotyped 151 are invisible to projection. Bigger long-read panel would lift this ceiling, but only if the new accessions overlap GrENE-Net 231 (the user has ~200 assemblies but most don't overlap; only +4 vs the current 80 by a quick check of `ASSEMBLIES_Best_version_of_dataset.csv`).

3. **231-simplex identifiability.** With many similar founder haplotypes, hapFIRE's CVXPY projection can have multiple near-optimal solutions. Per-ecotype Pearson r between hapFIRE freqs and the seedmix recipe is ~0 — hapFIRE doesn't precisely recover individual ecotype freqs, but aggregates (sums over multiple carriers) are accurate.

---

## Open next steps (ranked by value-to-effort)

### 1. Apply to real evolved GrENE-Net samples (½ day, no new SLURM)
The existing per-sample hapFIRE outputs are at `/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/samples/ecotype_frequency/MLFH*_ecotype_frequency.txt` (2,415 samples). The founder × SV matrix is built (`data/founder_sv_matrix.parquet`). Just project all 2,415 samples and compare to the existing freqk outputs in `freqk_gr/results/samples_p_svs.parquet`. Output: per-sample × per-SV table that slots into the GEA pipeline. **This is where the real value is — it's the actual deliverable, and only requires running `project_seedmix_svs.py` over a different input directory.**

### 2. Recombinant-pool simulation (1 week, tests the assumption that matters)
Modify the visor_freqk pipeline to build mosaic ecotype FASTAs (pick 2 founders, switch chunks at crossover positions). Pool 100-1000 mosaics, simulate reads, run hapFIRE, project. Truth = which mosaic carries which SV at which positions. This is the test that probes the no-recombination assumption directly. Spec:
- New script `build_mosaic_fasta.py` that takes 2 input FASTAs + a recombination rate and produces a mosaic
- Loop over N mosaics in `02_run_hack_var.sh`-style
- Reuse VISOR SHORtS for read simulation
- Reuse hapFIRE jobs and the projection script unchanged

### 3. Per-block ecotype freqs (½ day, immediate accuracy bump for real data)
Currently hapFIRE saves only the genome-wide weighted average ecotype frequency (top-20% by haplotype-diversity). Per-block freqs are computed but discarded (hapFIRE.py:138-139 are commented out). Patch hapFIRE.py to also save them, then make the projection use the local-window founder freq (assigning each SV to its containing LD block). This mitigates the case where a recombinant in one block doesn't pollute estimates elsewhere. About 5 lines of Python in hapFIRE.py + a small change in `project_sv_freq.py`.

### 4. LD-aware k-mer reweighting (1 week, fixes freqk's failure mode)
For each SV, find its 50 nearest SNPs in tight LD (already in the founder VCF). Use hapFIRE's per-SNP frequency as an *expected* AF for each SV-flanking k-mer. If freqk's k-mer count says 50% but linked SNPs say 5%, that k-mer is most likely contaminated by an off-target hit — drop it from freqk's count. In the seedmix this would automatically filter the AC=1 SVs where freqk reports AF=1.0. Inverse of HAF-pipe.

### 5. Bigger founder panel (depends on assembly availability)
The current panel uses 82 assemblies (80 unique 1001G IDs after collapsing). The `ASSEMBLIES_Best_version_of_dataset.csv` lists 534 unique Accession_IDs total but only 84 overlap GrENE-Net 231. If the user can point at a bigger collection of long-read assemblies that includes more GrENE-Net ecotypes, rebuild the SV panel:
- pipeline is at `/home/tbellagio/scratch/pang/sv_panel/` (alignment → syri → filter → bcftools merge)
- `build_founder_sv_matrix.py` then re-runs unchanged on the new merged VCF
- projection downstream is automatic

---

## What the user wanted next (their exact phrasing)

> "ok in the meantime, can you think of a more realistic way to test this?"
> "i have long read assemblies from at least 200 ecotypes, but the only overlap with grenent is 81"
> "could we use hapfire approach of learning form ecotypes to empower the kmer based methods based on arabiodpsis LD?"

Short answers:
- More realistic test → **option 2 (recombinant pool)** above. The simulation we have is best-case for hapFIRE; we haven't probed the recombination assumption.
- Bigger panel → useful but limited by GrENE-Net overlap (84 max from current data, vs 80 we have). Worth it only if the new ~200 assemblies add non-GrENE-Net-overlap ecotypes that we'd never project against anyway.
- LD-empowered k-mer → **option 4 above** is the clean version. **Option 1 (posterior combine)** is the half-day version and is already implemented (`scripts/posterior_combine.py`).

---

## File map (current — 2026-05-07)

```
hapfire_sv/
├── HANDOFF.md                        ← this file (TL;DR + what's pending)
├── FINAL_RESULTS.ipynb               ← AUTHORITATIVE cov50 head-to-head
├── FINAL_RESULTS_cov10.ipynb         ← AUTHORITATIVE cov10 + coverage-robustness
├── RECOMB_SWEEP_RESULTS_slim.ipynb   ← full method × LD-tier × var_type sweep
├── RECOMB_SWEEP_RESULTS.ipynb        ← extended (Tier 1-9.6) sweep notebook
├── SIMULATION_RESULTS.ipynb          ← original no-recomb simulation (Chr1 hapFIRE-proj)
├── BIGLD_BLOCKS_DISTRIBUTION.ipynb   ← BigLD block-size distribution
├── METHODS_TRIED.md                  ← design-space map + literature analogs
├── RESULTS_LOG.md                    ← chronological reasoning log
├── MODEL_SPEC.md                     ← Poisson generative model + EM
├── ALGORITHM.md                      ← solver implementation walkthrough
├── CACTUS_EM_MATH.md                 ← cactus_em window/global/anchor math
├── CACTUS_PANGENOME_AUDIT.md         ← cactus-vs-syri panel reconciliation
├── BACKGROUND.md                     ← project framing
├── SUMMARY.md                        ← HISTORICAL — original hapFIRE-proj writeup
├── archive/                          ← superseded methods + lost benchmarks (see archive/README.md)
│   ├── README.md
│   ├── route3_kallisto_em/           ← kallisto-style EC EM (LOST: −6 pp R² @ 10kb)
│   ├── bench_kmc/                    ← KMC + rapidgzip bench (NO WIN; 36 GB)
│   └── hapfire_projection/           ← original Chr1 hapFIRE-proj deliverable (SUPERSEDED)
│       ├── SIMULATION_RESULTS.ipynb
│       ├── scripts/                  ← original projection scripts (12 files)
│       ├── data/                     ← founder_sv_*.parquet, vcf/, seedmix_recipe.tsv
│       └── results/                  ← rep*, simulation_*, seedmix_*, hapfire_*svs*
├── poolfreq/                         ← cactus_em production code
│   ├── src/
│   │   ├── per_sample_per_chrom.py   ← PRODUCTION DRIVER (per-chrom, 64GB-friendly)
│   │   ├── per_sample_driver.py      ← genome-wide driver (legacy, 200GB)
│   │   ├── per_sample_bigld_haplotype.py  ← bigld_haplotype mode (in-progress)
│   │   ├── per_sample_kallisto_em.py ← Route 3 prototype (failed, kept for postmortem)
│   │   ├── em_solver.py              ← Poisson EM with global-anchor KL prior
│   │   ├── block_em.py               ← per-window EM + project_blocks_to_records
│   │   ├── block_haplotype_em.py     ← bigld_haplotype EM + smooth_h_across_blocks
│   │   └── kmer_count.py             ← jellyfish wrapper
│   ├── scripts/                      ← cn-builders, panel utilities
│   ├── tests/                        ← validation + per-sample SLURM templates
│   ├── bench_kmc/                    ← KMC vs jellyfish bench (no win)
│   └── data/
│       ├── cn_full_231_v2/           ← founder × k-mer panel (per-chrom)
│       ├── cn_var_231_v2.cn_var.npz  ← founder × variant projection (Beagle-imputed)
│       └── cn_var_231_v2.meta.npz    ← variant metadata
├── sims/visor_freqk/                 ← VISOR HACk + SHORtS recomb sims
│   ├── pool_sweep_82_recomb/         ← 4 regimes × cov50/cov10
│   ├── scripts/                      ← sim builders, method runners, eval
│   └── RECOMB_SIM.md                 ← sim pipeline doc
├── pangenie_genotyping/              ← sister pipeline (PanGenie per-individual)
├── pangenome_comparison/             ← pang_82 vs pang_135
├── preprocess_qc/                    ← per-sample read QC + LD-comparison notebook
├── imputation/                       ← Beagle 5.5 (kept on disk; not used in production)
├── contamination_test/               ← planned 71-LR-assembly mix-in (Phase A/B)
├── notebooks/                        ← ad-hoc analysis (blocks_snp_sv_density, etc.)
├── plots/                            ← rendered figures
├── results/
│   ├── session_summary/              ← af_long.parquet + af_summary.parquet (notebook inputs)
│   ├── apples_decomposed_summary.tsv, blocks_snp_sv_density.tsv.gz,
│   │   diagnostic_disagreement/, hapfire_perblock_vs_cem_by_ldtag.tsv,
│   │   hapfire_proj_vs_cem_svs_summary.tsv, recomb_window_size_compare.tsv,
│   │   contamination_preview/        ← cactus_em-era ancillary outputs
│   └── ...                           ← (original hapFIRE no-recomb sim + SEEDMIX
│                                       outputs moved to archive/hapfire_projection/results/)
├── data/                             ← shared inputs (sv_panel_to_accession_id, vcf_samples_231,
│                                       seedmix_recipe_normalized — used by imputation/)
└── logs/                             ← SLURM .out/.err
```

The `poolfreq/`, `sims/`, `pangenie_genotyping/`, `imputation/`, `preprocess_qc/`, and `contamination_test/` trees are the active codebase. Historical / superseded / lost-benchmark artifacts live under `archive/` (see `archive/README.md`).

---

## Reproducibility notes

- **hapFIRE source:** `/carnegie/nobackup/scratch/xwu/haplotype_frequency_estimation/hapFIRE_sourcecode/` (Xing Wu's, unmodified).
- **HARP binary:** `…/hapFIRE_sourcecode/bin/harp`. Must be on PATH for hapFIRE to call it.
- **Conda env:** `hapfm` (`/home/tbellagio/miniforge3/envs/hapfm/`). Has cvxpy, pandas, scipy, sklearn, networkx, python-louvain (the last installed during this session, was missing).
- **VCF sample list:** the GrENE-Net VCF has 232 columns; the 232nd is a blank trailing column. Drop with `bcftools query -l … | grep -v "^$"` to get the actual 231.
- **Memory:** hapFIRE on full Chr1 with 231 ecotypes needs **96 GB**. 16 GB OOM'd at HARP stage.
- **Time:** ~1.5 h per BAM at cov50× on a single 4-CPU node. Two LD blocks on Chr1 (24 Mb + 6 Mb).
- **Chromosome naming:** the GrENE-Net VCF originally uses chrom "1"; the simulator uses "Chr1". We use a renamed VCF (`greneNet_final_v1.1.recode.Chr1.vcf.gz`) to match.
- **Posterior combine:** `posterior_combine.py` uses `SIGMA_ECO=1e-4` (calibrated from rep9 cov50 f30 per-ecotype error) and a placeholder `n_kmer=100` for freqk; on real data swap in `n_unique_kmers` from the freqk index file.

---

## Known issues / loose ends

- **`compare_methods.py` is legacy** (only handled rep9 cov50). Superseded by `plot_simulation.py` which does the full sweep. Keep or delete.
- **`run_all_analysis.sh` calls `seedmix_truth_analysis.py`** which prints the seedmix 3-way comparison. Currently the simulation is hardcoded to k=31; if you want k=21/41/etc. the freqk output paths in `plot_simulation.py` need updating.
- **VCF for hapFIRE is uncompressed (798 MB)** at `data/vcf/greneNet_Chr1_only.vcf`. Re-bgzip + tabix if disk is tight.
- **rep1 BAMs are partial** — only f10 and f50 BAMs exist for rep1 cov50. So the cov50 dataset has 15 results across {rep1: 2, rep9: 5, rep10: 5, rep11: 3}.
- **No simulation results for cov10/cov20 on reps other than rep9** — only rep9 has the full 5-freq sweep at all 3 coverages. Adding rep10 cov10/cov20 would make the low-coverage MAE estimates more robust.

---

## Quick commands to pick up

```bash
# refresh the simulation plot (after any new hapFIRE jobs land)
/home/tbellagio/miniforge3/envs/hapfm/bin/python \
    /home/tbellagio/scratch/hapfire_sv/scripts/plot_simulation.py

# rerun seedmix 3-way comparison
/home/tbellagio/miniforge3/envs/hapfm/bin/python \
    /home/tbellagio/scratch/hapfire_sv/scripts/seedmix_truth_analysis.py

# submit a new hapFIRE job
sbatch /home/tbellagio/scratch/hapfire_sv/scripts/run_hapfire_one.sh \
    /path/to/sim.srt.bam /path/to/output_prefix

# project an existing hapFIRE output onto SVs
/home/tbellagio/miniforge3/envs/hapfm/bin/python \
    /home/tbellagio/scratch/hapfire_sv/scripts/project_sv_freq.py \
    --ecotype-freq <path>/MLFH010120180409_ecotype_frequency.txt \
    --vcf-samples /home/tbellagio/scratch/hapfire_sv/data/vcf_samples_231.txt \
    --n-samples 231 --sv-freq 0.5 \
    --out /tmp/test.tsv
```

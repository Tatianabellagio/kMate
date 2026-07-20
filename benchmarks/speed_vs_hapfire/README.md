> ℹ️ **See [`PANEL_MISMATCH_BUG.md`](PANEL_MISMATCH_BUG.md).** An earlier
> poolsize×depth accuracy/speed comparison ran kMate on the arch3 panel but
> hapFIRE on the old greneNet SNP-only VCF, while reads came from arch3 — a
> panel mismatch that invalidated the pre-2026-07-15 `hapfire_vs_kmate_*`
> results (archived under `archive/panel_mismatch_prefit_2026-07-15/`). **Fixed
> 2026-07-15**: each tool now runs on its own native panel with reads simulated
> to match (kMate: arch3; hapFIRE: greneNet-derived FASTAs + the greneNet VCF),
> same founders-meta + seed so pool composition is identical across tools. See
> "Poolsize × depth (fair...)" below for current numbers.

# kMate vs hapFIRE — running-time benchmark

**Date:** 2026-06-03. **Question:** how much faster is kMate than hapFIRE (= a
wrapper around HARP) at estimating per-record allele frequencies from a pool-seq
sample? No prior head-to-head *runtime* comparison existed (only accuracy ones),
so this measures it directly on identical input.

## Setup (identical input, same node type, matched cores)

- **Pool (same simulated sample for both tools):**
  `visor_freqk/data/reads_var/del/rep1/cov50/var_del_1kb_n231_f30_err001/`
  — 231 founders pooled, ~38.6× on Chr1, 1 kb deletion at f=0.30.
  hapFIRE consumes `sim.srt.bam`; kMate consumes `r1.fq`/`r2.fq` (same pool).
- **Scope:** Chr1 only (the simulated reads are Chr1-only).
- **kMate panel:** production arch3 matrices `kmer_pa_231_arch3_filt2inv` (K=10.95M
  k-mers) + `var_pa_231_arch3_chr1` (2,154,423 records). Recipe: `--kmer-weight
  uniform --unit chrom` (scripts updated; the old `inv_mb --block-mode global`
  is superseded — `--block-mode global` is a deprecated alias for `--unit chrom`),
  `--threads 8`, BLAS/OMP pinned to 8.
- **hapFIRE panel:** `greneNet_final_v1.1` subset to Chr1 = **827,765 phased SNPs ×
  231 founders**. BAM contigs reheadered `Chr1→1` to match the VCF + Ensembl
  reference (TAIR10, identical sequence/coords). Phases 1–2 only (`-haplotype`
  omitted → BigLD/Phase-4 skipped); HARP `like`/`freq` over the python
  independent-LD partition produces `_snp_frequency.txt`.
- **Hardware:** both jobs ran on **savio3 node n0083** (32-core node, 8 cores
  requested), `/usr/bin/time -v`.
- **Fairness:** one-time precomputes are excluded for both (kMate's matrix build;
  hapFIRE's BigLD partition build — and BigLD never ran here at all). hapFIRE's
  parallelism is commented out in this build (`haplotype_generation.py:301-307`
  calls HARP serially), so it is single-threaded by design; kMate uses its 8 threads.

## Result

| Condition | kMate | hapFIRE | kMate advantage |
|---|---:|---:|---|
| **Single core (matched — the fair number)** | **8 min 53 s** (533 s) | **25 min 50 s** (1550 s) | **≈ 2.9× faster** |
| As-run wall-clock (kMate 8 threads vs hapFIRE 1) | 3 min 00 s (180 s) | 25 min 50 s (1550 s) | 8.6× |
| Total CPU time (user+sys) | 8 min 51 s (531 s) | 24 min 37 s (1477 s) | ≈ 2.8× |
| Peak RSS (process) | 24.2 GB | 41.0 GB | 1.7× less memory |
| Records output | 2,154,423 | 827,765 SNPs | — |

**Fairest claim (both single-threaded): kMate ≈ 2.9× faster.** hapFIRE's
parallelism is disabled in this build, so the larger 8.6× wall-clock gap reflects
kMate also using 8 threads. CPU-time (2.8×) is parallelization-agnostic and tracks
the single-core result. kMate single core: 319 s k-mer count + 157 s EM = 531 s.
(8-core kMate, BLAS pinned to 8: 3 min 00 s wall, 13 min 12 s CPU.)

kMate stage breakdown: 78 s k-mer count (jellyfish) + 45 s EM (200 iters) +
projection → 157–181 s total.

## Poolsize × depth (fair, each tool on its native panel)

Separate from the single-condition run above: a full N × depth × seed sweep
(N ∈ {2,5,20,50,150}, depth ∈ {1,10}×5 seeds, plus N=50 × depth ∈ {30,50}×5
seeds — 60 conditions), comparing h (founder-mixture) accuracy, founder
recovery, and speed/compute. kMate runs on arch3 (existing `benchmarks/p231`
results, no rerun); hapFIRE runs on greneNet-derived FASTAs + the greneNet VCF
(`sims_greneNet/`, `results/greneNet_fair/`). Both draw the identical pool
(same arch3 founders-meta + seed) so the comparison is apples-to-apples. No
shared AF-accuracy comparison — see `PANEL_MISMATCH_BUG.md` for why.

```bash
scripts/submit_greneNet_fair_grid.sh                       # sim + hapFIRE, full grid (idempotent)
/global/home/users/tbellg/miniforge3/envs/basic/bin/python scripts/score_hapfire_vs_kmate.py
/global/home/users/tbellg/miniforge3/envs/basic/bin/python scripts/plot_hapfire_vs_kmate.py
/global/home/users/tbellg/miniforge3/envs/basic/bin/python scripts/plot_speed_vs_coverage.py
```

Outputs: `results/hapfire_vs_kmate_table.tsv`,
`results/hapfire_vs_kmate_{h_accuracy,founders_recovered,speed,speed_vs_coverage}.png`.

### Result (60 conditions, 2026-07-15)

Averaged over depth ∈ {1,10}×5 seeds per N (the main grid):

| N | kMate R²(h) | hapFIRE R²(h) | kMate founders | hapFIRE founders | wall kMate | wall hapFIRE | speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 0.994 | 0.980 | 100 % | 100 % | 102 s | 979 s | 9.6× |
| 5 | 0.994 | 0.986 | 100 % | 100 % | 102 s | 979 s | 9.6× |
| 20 | 0.981 | 0.930 | 100 % | 100 % | 105 s | 1009 s | 9.6× |
| 50 | 0.950 | 0.866 | 99.0 % | 99.6 % | 106 s | 999 s | 9.5× |
| 150 | 0.580 | 0.723 | 96.7 % | 98.7 % | 93 s | 965 s | 10.4× |

- **Founder recovery is essentially tied** — both tools recover ~all pooled
  founders (top-K-vs-truth) across the whole grid; hapFIRE is marginally ahead
  at N=150 (98.7 % vs 96.7 %). The pre-fix "hapFIRE collapses to 36/50" result
  was the panel-mismatch artifact and is gone.
- **Proportion accuracy (R² of h):** kMate leads at N ≤ 50 (e.g. 0.95 vs 0.87 at
  N=50); at **N=150 the ordering flips** — hapFIRE 0.72 vs kMate 0.58. Both
  degrade as the pool approaches full-panel density (founder weights → near-
  uniform ~1/N, so R²'s denominator, sd(truth), shrinks and the metric gets
  fragile), but kMate degrades faster there. This N=150 reversal is real and
  worth flagging — not an artifact of the old bug.
- **Speed:** kMate ~**9.6× faster wall-clock** (median over all 60) and ~4.5×
  fewer CPU-seconds (the parallelization-agnostic number — kMate uses 8 threads,
  hapFIRE's HARP stage is single-threaded by build). At fixed N=50 the wall-clock
  gap holds ~9–10× as depth climbs 1×→50× (`..._speed_vs_coverage.png`).
- **Memory:** here kMate's peak RSS is *higher* (~32 GB vs ~18 GB) — the reverse
  of the single-condition 231-founder run above — because kMate carries the full
  arch3 k-mer + variant matrices while hapFIRE runs on the smaller greneNet
  SNP-only panel. Panel size, not method, drives this.

## Reproduce (single-condition speed run above)

```bash
sbatch benchmarks/speed_vs_hapfire/run_hapfire.sh   # hapFIRE (env: hapfire)
sbatch benchmarks/speed_vs_hapfire/run_kmate.sh     # kMate   (env: kmate)
# timings in results/{hapfire_time.txt,kmate_time_pinned8.txt}
```

- `hapfire` env: Python-only (cyvcf2/numpy/pandas/scipy/sklearn/networkx/
  python-louvain/cvxpy/pyclustering) + the statically-linked `harp` binary
  (`external/HapFIRE/src/harp`). No R needed because BigLD is not run.
- Input prep (one-time) lives in `work/`: Chr1-subset phased VCF, reheadered BAM,
  Chr1 partition subset.

## Caveats / scope

- Single representative point (rep1, cov50, f30, 1 kb del). For a paper figure,
  repeat across reps × coverages (`run_*.sh` parameterize trivially).
- hapFIRE timed through per-SNP frequency (Phase 2). Its optional fine-haplotype
  step (Phase 4) is excluded; it needs either the prebuilt partition via `-block`
  or R for BigLD, and is fast python post-processing relative to the HARP stage.
- Accuracy was validated separately (both produce valid AF); this run measures
  runtime only.

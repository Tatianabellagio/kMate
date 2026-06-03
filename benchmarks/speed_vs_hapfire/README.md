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
  inv_mb --block-mode global`, `--threads 8`, BLAS/OMP pinned to 8.
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

## Reproduce

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

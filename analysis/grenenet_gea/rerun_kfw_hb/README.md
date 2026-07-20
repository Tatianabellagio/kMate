# rerun_kfw_hb — GrENE-Net founder-frequency cohort (2026-07-07)

## What this is
Full re-estimation of the GrENE-Net cohort founder mixtures `h` under the
**full-panel Kf_w normalization fix** (commit `9669be7`) + the unified
`--unit` estimator (`a9bf1c0`). Supersedes `rerun_perfounder/` (built
2026-07-06 22:41, **before** the Kf_w-fullpanel fix landed 2026-07-07 10:22 —
i.e. stale).

- `seedmix/`  — 8 SEEDMIX reps (founding p0)
- `evolved/`  — 2168 evolved pool samples

## Estimator used: `--unit chrom` (whole-chromosome), NOT `--unit ld`

**The cohort was estimated with `--unit chrom` — one founder mixture `h` per
whole chromosome — for every sample.** This is the correct estimator for a
selfing pool (see below). It was launched via `grenenet/run_site_array_perchrom.sh`
with `BLOCK_MODE=global KMER_WEIGHT=uniform NORMALIZE=per_founder`; the driver
maps the (deprecated) `--block-mode global` → `--unit chrom`.

### Why chrom, not ld/r²=0.1 (the decisive evidence)
The r²=0.1 LD-block ("haploblock") path was explored and **rejected for the
selfing path**. On the equimolar `n231_g0` benchmark (true `h` uniform):

| estimator | AF-MAE | outliers >0.10 |
|---|---|---|
| `--unit chrom` | **0.0033** | 0.001% |
| `--unit ld` (r²=0.1) | 0.0080 | **0.625%** (625×) |

Mechanism: per-block EM re-exposes founder **non-identifiability in
low-diversity regions** (centromere/pericentromere). Those blocks have plenty
of k-mers but cannot distinguish 231 founders from *local* k-mers alone, so the
per-block EM drifts to a sparse vertex (eff_n 6–34 vs true ~231) and sprays
every record in the block. `--unit chrom` pools ALL of a chromosome's k-mers,
so diverse regions pin the founders down and low-diversity records inherit the
well-determined genome-wide `h`. Under ~97% selfing, ancestry is constant
genome-wide, so there is nothing to gain from per-block fitting — only the
centromere non-identifiability to lose. **ld is reserved for genuinely
recombinant pools.**

(Note: at r²=0.1 the within-block haploblock collapse is a no-op anyway — K_b≈231
genome-wide — so "r²=0.1 haploblocks" reduced to chrom in the frequencies, and
`ld` differs from `chrom` only via the harmful per-block *fitting*, not the
collapse.)

## Verification that the cohort ran chrom (3 independent checks, 2026-07-07)
Concern: was this launched on r²=0.1 blocks by mistake? No — verified:
1. **Structure:** each `*_ChrN.h_per_chrom.npz` holds a single `(231,)` vector
   per chrom (chrom fingerprint); **no** `*_h_blocks` per-block matrices (the ld
   fingerprint).
2. **Logs:** 10,828 task-log lines report `per-chrom unit=chrom`; **zero** report `ld`.
3. **File census:** 10,840 `.h_per_chrom.npz` (= 2168×5, chrom); **0**
   `.h_blocks_per_chrom.npz` (ld).

And `--block-mode global`→`--unit chrom` was independently shown byte-identical
(max|Δ|=0) to the pre-refactor full-panel-Kf_w global estimator (e2e job 35562432).

## Trap to be aware of
The **driver default is `--unit ld`**, which is WRONG for selfing pools. The
production runner is safe (it passes `--block-mode global`→chrom), but a bare
`per_sample_per_chrom` invocation on a selfing pool would silently use the bad
per-block path. Prefer `--unit chrom` explicitly for this project's pools.

## Config (production recipe)
- panel `data/kmer_pa_231_arch3_filt2inv` (in-house build), V_pa `panel/arch3` per-chrom triplet
- `--unit chrom` (via `--block-mode global`), `--kmer-weight uniform`, `--normalize per_founder`
- founder 9977 (previously ~4e-6, effectively zeroed pre-fix) is now defined
  everywhere; min founding p0 = 1.5e-4 → the old `P0_FLOOR=1e-4` is now a no-op.

## Status / next
Cohort complete (2168/2168, 8,489,646 records each). Integrity check =
`grenenet/validate_cohort.py`. Downstream still pending: repoint `lib.OUT`/
`lib.SEEDMIX` here, rebuild the AF store, re-run the selection chain (evolved
compositions are more concentrated under the Kf_w fix, so Δh/selection shift).
`rerun_perfounder/` retained until this is validated.

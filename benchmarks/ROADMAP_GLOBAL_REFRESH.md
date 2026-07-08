# Benchmark refresh roadmap — new corrected GLOBAL production pipeline (2026-07-07)

Scope: regenerate the **p231** and **p80** benchmarks (and the global-affected
auxiliary benchmarks) for the **new, fully-corrected GLOBAL estimator**.
**Window mode is explicitly OUT of scope** for this refresh — do not touch the
10 kb / star2 window recipe, its scripts, or its result sets.

---

## 0. What changed in production GLOBAL (why the benchmarks are stale)

The committed global estimator has changed on three axes since the current
benchmark numbers were produced. All three affect global:

1. **Per-founder M-step normalization (`normalize="per_founder"`, full-panel `Kf_w`)**
   — commits `a8ba02d` + `9669be7`. Removes the completeness bias; fixes the
   founder-`h` collapse. AF-MAE ~2.7× better in the ablation.
2. **Haploblock collapse (`eps=0`)** — commit `b4d6ce0`. Fit the `K_b ≤ F`
   distinct k-mer haplotypes per unit, not all 231 founders; equal-split back.
3. **"Global" is retired as a distinct mode; the estimator is unified under `--unit`**
   (this refresh, decided 2026-07-07). With haploblock collapse + per-founder
   normalization + local-only, `global` and `window` are the SAME algorithm at
   different unit sizes: `partition → per-unit (collapse→local EM) → project`.
   Production is now `--unit ld --ld-r2 0.1` (CompleteLDPartition r²=0.1 blocks from
   var_pa_231_arch3, ~18/Chr1, variable size; K_b≈231 inside each so collapse ≈ no-op
   but the unit is data-derived, not an arbitrary whole-chromosome). **This unification
   is Phase 0 and must be built first.**

Also switching the benchmark k-mer panel to the **real production in-house index**
`data/kmer_pa_231_arch3_filt2inv` (not the benchmark-local pang_135 rebuild), so the
refresh is a true end-to-end production benchmark.

**Production config the benchmarks must run:**
```
--unit ld --ld-r2 0.1        r²=0.1 CompleteLDPartition blocks (from var_pa_231_arch3)
--haploblock-eps 0           exact haploblock collapse per unit, then local EM
--normalize per_founder      full-panel Kf_w  (default)
--kmer-weight uniform        (ω=1/m_b retired)
--kmer-pa-prefix data/kmer_pa_231_arch3_filt2inv/kmer_pa   (production in-house index)
--var-pa panel/arch3/chr1/var_pa_231_arch3_chr1[_atomized]  (unchanged)
per unit: haploblock collapse (eps=0) → local EM → project the unit's records
NO anchor prior, NO cross-window smoothing, NO fallback
```

---

## What is REUSED (do NOT regenerate)

- **Simulated reads + truth tables** — estimator-independent. Reuse every existing
  sim regime under `benchmarks/{p231,p80}/sims/`. Truth stays fixed so old↔new
  numbers are directly comparable.
- **Competitor outputs (hapFIRE, vg giraffe)** in `accuracy_vs_competitors/` — these
  are competitor-side and independent of kMate's changes. **Do NOT re-run hapFIRE
  or vg.** Only re-run the kMate arm and re-score.
- **var_pa / var_called** (arch3) — unchanged, already production.
- **`fastas_231` / p80 fastas** — sim inputs, unchanged.

## What is RE-RUN

Estimator (kMate global) → score → notebook/plots, on the existing sims, for every
benchmark family below.

---

## Phase 0 — Unify the estimator under `--unit` + panel (PREREQUISITE, do first)

**0a. Unify global/window into one `--unit`-parameterized estimator** (decided
2026-07-07: full unify). There is one algorithm — `partition → per-unit
(haploblock-collapse → local EM) → project` — so collapse the two driver functions
into it and select the unit with a single flag:

```
--unit ld    [--ld-r2 0.1]     CompleteLDPartition LD blocks   (DEFAULT, r2=0.1)
--unit chrom                   one unit = whole chromosome     (former "global")
--unit bp    [--window-bp N]   fixed-bp windows                (former "window")
--unit tsv:<path>              explicit block TSV
```

Implementation:
  - **One code path.** `run_one_chrom_global` and `run_one_chrom_window` merge into a
    single `run_one_chrom(unit, ...)`. `chrom` is the degenerate partition `[whole
    chromosome]` fed to the same per-unit block-EM path that already does eps=0
    haploblock collapse + local fit. No anchor / no smoothing / no fallback is the
    only behavior (local-only becomes intrinsic, not a flag).
  - **`--unit ld`** computes CompleteLDPartition blocks from `var_pa` at `--ld-r2`
    (reuse `gen_ld_partitions.py`'s partitioner; cache the block TSV per panel/chrom).
  - **Port the two global-only features onto the unified path:** `--h-only` (skip the
    projection; trivial) and `--emit-af-se` (already collapse-aware — compute the
    K_b-class Fisher covariance per unit and equal-split-map to founders; lift the
    logic from the current global path). Both are per-unit now.
  - **Back-compat aliases (do not break production / the running cohort):** keep
    `--block-mode global` → `--unit chrom` and `--block-mode window` → `--unit bp`
    as deprecated aliases. The grenenet selfing cohort currently invokes
    `--block-mode global`; it must keep running byte-identically (unit=chrom is the
    exact whole-chromosome path, K_b==F no-op).
  - **Verify invariants after the merge:** `kmate selftest` PASS; `--unit chrom`
    byte-identical to the old `--block-mode global`; `--unit bp` reproduces the old
    window path; `--unit ld` runs and logs per-block K_b.
  - Update `cli.py`, driver docstring, ALGORITHM.md §5/§7, PIPELINE_STATE.md,
    README once merged. Window-mode *recipe* content (anchor/smoothing) is not
    reintroduced — it stays as the deprecated `--no-local-only` path only if kept.

**0b. Generate the r²=0.1 LD blocks per chromosome** (already done for Chr1–5):
`results/grenenet_gea/hap_blocks/ld_blocks_r2_0.10_<Chr>.tsv` (+ `_genome.tsv`),
from `gen_ld_partitions.py` (CompleteLDPartition on `var_pa_231_arch3`, MAF≥0.05,
callrate≥0.9). Benchmarks are Chr1-only → use `ld_blocks_r2_0.10_Chr1.tsv`
(18 blocks). Copy/symlink into each benchmark's `data/` for provenance.

**0c. Point benchmarks at the production in-house-index panel.**
Founder order verified identical across `data/kmer_pa_231_arch3_filt2inv`,
`var_pa_231_arch3`, and the benchmark-local build (231, same order) — safe to swap.
  - p231: set `--kmer-pa-prefix data/kmer_pa_231_arch3_filt2inv/kmer_pa` (replaces
    the benchmark-local `kmer_pa_p231_filt2inv`).
  - p80: production panel is the 231-founder arch3 index; the p80 control needs its
    OWN 80-cactus-founder filt2inv in-house index. **Decision:** either (i) build a
    p80 filt2inv in-house-index kmer_pa from the p80 VCF (matches production
    construction on the homogeneous panel), or (ii) keep p80 on its existing
    `kmer_pa_p80_filt2` and treat p80 as the balanced-panel *control* only. p80's
    role is isolating panel imbalance, not shipping a production panel — (ii) is
    acceptable if you just want the control; (i) if you want p80 fully
    production-faithful. **Confirm.**

**0d. Verification gates** (must pass before trusting any number):
  - founder order: kmer_pa == var_pa (verified for p231; re-check for p80 build).
  - est↔truth join == var_pa record count (100%), per existing G6 gate.
  - `K_b` printed per block ≈ 231 on the p231 panel (sanity: collapse ≈ no-op at
    r²=0.1); log the per-block `K_b` distribution.

---

## Phase 1 — p231 AF accuracy (headline)

Files: `benchmarks/p231/`. Regimes (reuse sims): `n50_g0 n231_g0 n50_g1 n231_g1
n50_g3 n50_g3_dom500` (+ any `_self97` variants already simulated). Both var_pa
arms: `raw` (SNP/indel/SV) and `atomized` (SNP-level).

1. **New runner.** Adapt `scripts/07f_run_kmate_filt2inv_p231.sh` to the Phase-0
   config: production in-house kmer_pa + `--unit ld --ld-r2 0.1 --kmer-weight uniform`
   (replacing `--block-mode global`). New output dir e.g.
   `results/kmate_ldr01_<cnvar>/<regime>/`.
2. **Score.** `scripts/score_p231.py` — MAE/RMSE/R²/outlier by regime × var-class
   (SNP/indel/SV). Point it at the new output dirs; keep the old dirs for the
   old↔new delta column.
3. **Notebook + plots.** Re-execute `results/FINAL_RESULTS_p231*.ipynb`
   (`scripts/exec_nb.sh`) and regenerate the scenario-grid PNGs via
   `plot_scenario_grids.py` + `_accuracy_panel.py`
   (`results/plots/p231_{scenario}_global[_miss50|_miss90].png`).
4. **Update** `submit_all_p231.sh` Phase C/D to call the new runner and **drop the
   legacy ω=1/m_b (`inv_mb`) A/B arm entirely** — it is stale: per-founder
   normalization now removes the panel-completeness imbalance at its source, so the
   m_b de-replication weight is redundant for global (confirmed 2026-07-07). Global
   benchmarks run `--kmer-weight uniform` only; the weight axis disappears from
   `score_p231.py`'s grid.

## Phase 2 — p80 AF accuracy (control)

Files: `benchmarks/p80/`. Same as Phase 1 with p80 scripts
(`07c_run_kmate_filt2_mb_p80.sh` → new global-ld runner, `WEIGHT=uniform`),
`score_filt2_mb_vs_uniform.py` / `score_filt2_mb_missingness.py`, and
`FINAL_RESULTS_cov10_p80.ipynb`. Panel per Phase-0c decision. p80 remains the
balanced-panel control (checks the cactus-vs-PG asymmetry is absent).

## Phase 3 — Founder-h recovery accuracy (most affected)

Files: `benchmarks/h_accuracy/`. The per_founder + haploblock changes affect the
founder **decomposition** more than AF, so this is the highest-signal refresh.
1. Run the new global estimator **h-only** (`--h-only`, writes `h_per_chrom.npz`)
   on `n231_g0` (and `n50_g0`, the skew/selection regime) — reuse sims.
2. `score_h_vs_truth.py --est <h_per_chrom.npz> --truth
   benchmarks/p231/sims/<sim>/pool_weights.tsv --label global_ldr01_<regime>`
   (run in `basic` env). Metrics: per-founder h RMSE/MAE, absorbed count,
   cactus-vs-PG mass balance, est-vs-true slope.
3. Compare against the retained `BASELINE_n231_g0.csv` and the existing
   `filt2mb`/`filt2u`/`filt2invu` arms already scored here. Expected: 0 absorbed,
   cactus mass ratio ≈ 1.0, slope ≈ 1.0.

## Phase 4 — Ecotype-count / resolution

Files: `benchmarks/ecotype_count/`. Uses `run_h_only.sbatch` + `score_ecotype.py`
(and the block-kmer-floor sweep). Re-run the kMate h-only arm under the new global
config; re-score resolvable-founder / ecotype counts. Ties directly to the
haploblock reframe (how many haplotypes the r²=0.1 blocks resolve). Reuse the
existing hapFIRE ecotype reference (`run_hapfire_ecotype.sbatch` output) — do NOT
re-run hapFIRE.

## Phase 5 — vs competitors (kMate arm only)

Files: `benchmarks/accuracy_vs_competitors/`. **Reuse existing hapFIRE + vg
outputs.** Re-run only the kMate arm under the new global config, then re-score
with `scripts/score_snp_fair.py` (SNP parity vs hapFIRE) and the SV-vs-vg path.
Read `SCORING_RULES.md` + `BENCHMARK_DESIGN.md` first (the join/MAR-loop pitfalls).
Update `benchmark_table_4tool.tsv` / `benchmark_4tool_RMSE_*.png` kMate rows only.

---

## Cross-cutting / bookkeeping

- **New output dirs, don't overwrite** the old ones — keep both so every notebook
  can show an old→new delta and nothing is lost if the refresh needs iterating.
- **Scripts to edit:** the `07*` runners (config swap), `submit_all_{p231,p80}.sh`
  (DAG), score scripts (new input dirs), notebooks (new dirs). Window-mode scripts
  (`07d_*_window*`) are untouched.
- **Envs:** estimator runs in `kmate`; scoring/notebooks/plots in `basic`
  (matplotlib hangs in `plotting`); score_p231 currently uses `hapfm` — migrate to
  `basic` if `hapfm` is being retired.
- **Compute:** sbatch everything (savio4_htc / co_moilab / savio_lowprio); the EM
  over a full Chr1 panel loads ~10 GB dense — do NOT run in an interactive shell.
- **Sanity vs production:** the p231 global-ld run on `n231_g0` should reproduce the
  seed-mix-style result (0 collapsed founders, cactus ratio ≈ 1.0). If not, stop.

## Open checkpoints to confirm before I start executing

1. **Phase 0c:** p80 panel — build a p80 filt2inv in-house index (production-faithful)
   vs keep `kmer_pa_p80_filt2` (control-only).

## Resolved

- **Estimator unified under `--unit`** (2026-07-07) — full unify: one code path,
  `--unit {ld,chrom,bp,tsv}`, default `ld:0.1`; `--block-mode global/window` kept as
  deprecated aliases so production is byte-identical. "global" retired as a distinct
  mode. See Phase 0a.
- **ω=1/m_b (`inv_mb`) is retired for global** — per-founder normalization removes the
  imbalance at its source, so m_b weighting is stale (2026-07-07). Benchmarks run
  `--kmer-weight uniform` only; no m_b comparison arm.

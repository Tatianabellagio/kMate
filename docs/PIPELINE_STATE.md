# Pipeline state — production kMate (2026-05-22, addendum 2026-05-27)

Supersedes `old_docs/SESSION_*` and the prior `PIPELINE_STATE_2026-05-19.md`. This file is the single source of truth on what's production-ready vs in evaluation vs deprecated. Updated 2026-05-27 to record the closed-out k-mer-filter / EM-weighting decision (see new §0 below).

## §0 — Update 2026-05-27 (k-mer filter + EM weighting resolved)

The "k-mer filter OPEN" item below is **resolved**. Production decision:

- **cn_full**: `cn_full_231_v3qc_v3_filt2` (drop ac=1 singletons; see `ALGORITHM.md` §2.1). ★
- **EM weighting**: `--kmer-weight inv_mb`, i.e. $\omega_k = 1/m_b$ per-bubble de-replication (see `ALGORITHM.md` §4.2). ★
- **Driver invocation**: `per_sample_per_chrom.py --block-mode {global|window} --kmer-weight inv_mb`.
- **Caveat (panel-conditional, not propagated into the paper)**: on the homogeneous p80 control panel without the cactus/PG imbalance, $\omega_k=1/m_b$ slightly under-performs unweighted EM (+2% to +41% AF MAE). The flag `--kmer-weight uniform` preserves the unweighted MLE byte-identically; balanced-panel users can opt out. See `METHODS_TRIED_AND_RESULTS.md` §3 and `ALGORITHM.md` §10 M6.
- **Stale results**: every TSV that predates the `--kmer-weight inv_mb` switch needs to be regenerated; this is the same TSV cohort flagged in §3 below for the MAR-projection re-run.

The bullets below are kept verbatim for historical context; treat §0 as authoritative wherever it conflicts.

## TL;DR

For each panel record, estimate per-record ALT allele frequency from pool-seq reads via a per-sample **weighted** Poisson EM on the 231-founder simplex (production $\omega_k = 1/m_b$), then project through `cn_var` to per-record AF. Production pipeline as of today:

- **Panel**: `founders_231_v3qc_v3` (78 cactus + 153 PG)
- **Decomposition**: **arch3** (annotate_vcf + convert-to-biallelic), NOT `bcftools norm -m -any`
- **Matrices**: `cn_var_231_arch3` (standard, for SV-level analyses) AND `cn_var_231_arch3_atomized` (per-base, for SNP-level GEA)
- **Chrom scope**: Chr1 built; **Chr2–5 extension is the open production task**
- **K-mer filter / cn_full build**: ~~**OPEN**~~ → **RESOLVED 2026-05-27 (see §0 above)**: `cn_full_231_v3qc_v3_filt2` + EM weighting $\omega_k=1/m_b$. The prior text in this bullet (subsampProtect1 leading, sweep in progress) is superseded.
- **Projection**: MAR — `(h @ cn_var) / (h @ cn_var_called)` in both global and window modes. Patched 2026-05-21.
- **Modes**: both `global` (one h per chrom) and `★★` (window 10kb + global anchor 0.3 + HMM smooth 5α0.5) are production; pick per-regime.
- **Naming**: the method is **kMate** (algorithm + math: `ALGORITHM.md`). Legacy code, result-dir paths (`*/cactus_em_*`), and the `sims/visor_freqk` sub-repo still use the prior name `cactus_em`; with the k-mer-filter decision now closed, the rename is unblocked but not yet executed.

## What changed since the last pipeline-state doc (2026-05-19)

| Change | Why |
|---|---|
| Arch decomposition replaces `norm -m -any` for cn_var | annotate_vcf + convert-to-biallelic produces symbolic-ID biallelic catalog with +17pp hapFIRE-SNP coverage gain (`project_arch3_chr1_validation` memory) |
| `cn_var_atomized` added as a production deliverable | Per-base SNP catalog with carriers UNIONed across overlapping records. −23% RMSE, −43% outliers for SNP-level GEA (`project_arch3_atomization_result`) |
| MAR projection in window mode | Previously window mode silently treated `.` as REF (different from global mode's called-mask normalization). Inconsistent semantics fixed 2026-05-21. Star2 numbers prior to this patch are stale. |
| K-mer filter declared **undecided** | Prior doc incorrectly said mixed-loose was production. mixed-loose was evaluated, didn't win unanimously. Active testing: `filt2`. **(Superseded 2026-05-27: filt2 + ω_k=1/m_b chosen; see §0.)** |

## 1. Three persistent matrices feed the pipeline

| Matrix | Path | Shape (Chr1) | Purpose | Build script |
|---|---|---|---|---|
| `cn_full` | **`data/cn_full_231_v3qc_v3_filt2/cn_Chr1.{cn,meta}.npz`** (filt2, drop ac=1) | (231, ~22.7M raw → filt2 subset) | founder × k-mer (k=31, from PanGenie-index) | `build_cn_full_filt2_v3qc_v3_chr1.sh` |
| `cn_var` | `panel/arch3/chr1/cn_var_231_arch3_chr1.{cn_var,cn_var_called,meta}.npz` | (231, ~6.3M) | founder × biallelic-record carriers + called-mask | `panel/arch3/chr1/jobA5_build_cnvar.sh` |
| `cn_var_atomized` | `panel/arch3/chr1/cn_var_231_arch3_chr1_atomized.{cn_var,cn_var_called,meta}.npz` | (231, ~7.5M) | founder × per-base SNP (atomized) | `panel/arch3/chr1/jobD1_atomize_cnvar.sh` |

The k-mer dictionary is the **pang_135 PanGenie-index** (built once on the full 135-assembly pangenome). Not rebuilt per panel; reused across v3, v3qc, v3qc_v3.

## 2. Driver — `src/per_sample_per_chrom.py`

Reads FASTQs, runs k-mer Poisson EM on the simplex, projects through cn_var. Output TSV schema (post 2026-05-21 patch):

```
chrom  pos  ref_len  alt_len  alt_freq  info  n_called  se
```

- `alt_freq` = `(h @ cn_var)[r] / (h @ cn_var_called)[r]` — MAR projection
- `info` = h-weighted called mass at record r (h-dependent, varies with method)
- `n_called` = integer count of called founders at record r (h-independent panel QC)
- `se` = Wald SE using `n_called` as effective sample size

Both `--block-mode global` and `--block-mode window` (the `★★` recipe) use the same projection semantics. **After the 2026-05-26 cleanup the window-mode defaults *are* the ★★ recipe** (window-bp 10000, global-anchor-weight 0.3, hmm-smooth-passes 5, hmm-smooth-alpha 0.5), so `--block-mode window` alone reproduces it. **Pass `--kmer-weight inv_mb` in both modes** (production EM weighting; see §0 and `ALGORITHM.md` §4.2). The LD-block modes, overlapping windows, and the older k-mer-budget rebalancing / carrier-weighting / contamination-ω variants were archived to `src/archive/`; the authoritative file list + recipes are in `src/README.md`.

## 3. Outstanding production work

| Task | Why |
|---|---|
| Build `cn_var_231_arch3` for Chr2-5 | Current Chr1-only build is enough to validate; whole-genome is needed for downstream GEA |
| ~~Choose production k-mer filter~~ | **DONE 2026-05-27**: `filt2` + `--kmer-weight inv_mb`. See §0 and `METHODS_TRIED_AND_RESULTS.md` §0/§3. |
| Re-run SEEDMIX baselines under MAR star2 + arch cn_var + `--kmer-weight inv_mb` | Prior numbers used `norm -m -any` cn_var, the buggy "star2 treats `.` as REF" projection, AND unweighted EM; need fresh validation under the full production recipe. |
| Production scale-out on 2,415 evolved GrENE-Net samples | SLURM template in `tests/`; ~1.5–5 days at cluster-wide concurrency |

## 4. What's deprecated (do not use)

- `cn_full_231_v3qc_v3_mixedloose/` and `cn_var_231_v3qc_v3.{cn_var,meta}.npz` for new production runs — superseded by arch-decomposed cn_var.
- "Star2 treats `.` as REF" projection — buggy, patched 2026-05-21. Any star2 result TSV without `info`/`n_called`/`se` columns is from the buggy code.
- mixed-loose / mixed-conserv / mixed-strict filters — evaluated, not production.
- Beagle imputation for SVs — hard-rejected (−25 to −30pp concordance loss).

## 5. Companion docs (still current)

- `BACKGROUND.md` — project framing
- `ALGORITHM.md` — kMate algorithm, math & production wiring (code-verified single source of truth; supersedes the old prose doc and `CACTUS_EM_MATH.md`)
- `old_docs/CACTUS_EM_MATH.md` — superseded formal-math doc (folded into `ALGORITHM.md`)
- `SIMULATIONS_METHODS.md` — methods-ready description of the pool-seq simulation framework
- `INVESTIGATION_CN_VAR_DECOMPOSITION.md` — context for the arch decomposition switch
- `MISSINGNESS_231PANEL.md` — F_MISSING characterization on the production panel
- `PIPELINE_FASTQ_PREPROCESSING.md` — read-side preprocessing pipeline (trim, dedup)
- `archive/exploration/panel_overlap_135_vs_82/RESULTS.md` — panel composition analysis
- `RESULTS_LOG.md` — chronological numerical record (large file; useful historical reference, not authoritative for current state)

## 6. Subprojects (own subfolder READMEs)

- `benchmarks/p80/` — homogeneous 80-cactus-founder control experiment. See `benchmarks/p80/README.md`.
- `sims/visor_freqk/` — pool-seq simulation framework. See `sims/visor_freqk/README.md` and `RECOMB_SIM.md`.

## 7. Historical record

All docs about prior pipeline states / superseded decisions are preserved in `old_docs/`. Useful for archaeology; do not treat as authoritative.

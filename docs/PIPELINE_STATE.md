# Pipeline state — production kMate

**This file is THE single source of truth for what kMate consumes in production.**
ALGORITHM.md, HANDOFF.md, and README.md all defer to the block below. If anything
anywhere disagrees with §0, §0 wins — fix the other place.

---

## §0 — PRODUCTION INPUTS (CANONICAL) — last set 2026-05-30

> **There is exactly ONE production pipeline: `arch3`.** arch3 (A1–A5) *decomposes*
> the merged 231-founder pangenome panel (symbolic-ID, not `bcftools norm`) into
> the per-chrom canonical biallelic VCF **`merged_231_chr{N}_final.vcf.gz`**. That
> VCF — and nothing else — is what every downstream matrix is built from. Both
> production matrices derive from it.
>
> **The final VCF is SEGREGATING-ONLY (set 2026-06-02).** A4 drops monomorphic records
> (`bcftools view -e 'INFO/AC=0 || INFO/AC=INFO/AN'`) — ALT carried by no founder *or* by
> all called founders. Genome-wide **8,489,646** records (was 10,325,364 before the filter;
> 17.8% were monomorphic, almost all AC=0). V_pa is built 1:1 from this VCF, so it is
> segregating-only too. K_pa is unaffected (monomorphic variants contribute no k-mers;
> filt2inv already bounds founders-per-k-mer to [2,230]). Stats + the exact commands:
> `results/panel_stats/PANEL_STATS.md`.

| What | Production artifact (the ONLY thing consumed) |
|---|---|
| **Panel decomposition** | `arch3` A1–A5 → `panel/arch3/chr{N}/merged_231_chr{N}_final.vcf.gz` |
| **K_pa** (`kmer_pa`, founder×k-mer) | `build_kmer_pa.py` ← in-house index `panel/pangenie_index/pang_135_haploid/ours_Chr{N}` + **`merged_231_chr{N}_final.vcf.gz`** + `TAIR10.chr.iupacN.fa`, `--treat-missing-as-n --filter-production` → **`data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr{N}.{kmer_pa,meta}.npz`** |
| **V_pa** (`var_pa`/`var_called`, founder×variant) | `build_var_pa.py` ← **`merged_231_chr{N}_final.vcf.gz`** (segregating-only) → **`panel/arch3/chr{N}/var_pa_231_arch3_chr{N}.{var_pa,var_called,meta}.npz`** — 8,489,646 records genome-wide. (`_atomized` for SNP-level is **stale**: Chr1-only, predates the segregating filter — rebuild before use.) |
| **k-mer index** | in-house `panel/pangenie_index` builder (`ours_Chr{N}`), Level-A+B validated equivalent to PanGenie's. PG index kept only as a reference comparator. |
| **K_pa filter** | `filt2inv` (drop ac=1 singletons **and** ac=F invariants), applied inline by `--filter-production` |
| **EM normalization** | **`per_founder` (default, set 2026-07-06)** — M-step divides each founder by its own effective k-mer content $K^w_f$ (RSEM effective-length correction). Fixes global-mode founder-$h$ collapse (~19–49/231 founders driven to ~0). Legacy multinomial kept as `--normalize global` (deprecated). See `docs/FOUNDER_NORMALIZATION_FIX.md`, `ALGORITHM.md` §4.3. |
| **EM weighting** | **GLOBAL mode (the GrENE-Net production estimator — heavily-selfing cohort, so global not window; emits BOTH per-sample founder $\mathbf h$ and per-record SNP/SV AF):** `--kmer-weight uniform` (drop $\omega=1/m_b$) + keep filt2inv — **supersedes** the old "$\omega=1/m_b$ production" note (per_founder + uniform wins both AF and $\mathbf h$: AF-MAE 0.0036 vs 0.0096, 0 vs 33 absorbed). **WINDOW mode** (recombinant pools only, not the selfing production path): `--kmer-weight inv_mb` as before. |
| **Driver** | `per_sample_per_chrom.py --unit {chrom|ld|bp|tsv} --normalize per_founder` — `--unit chrom` (the default; production selfing path) adds `--kmer-weight uniform`; the recombinant `--unit bp` path keeps `--kmer-weight inv_mb`. (`per_founder` is the default; `--normalize global` reproduces legacy. `--block-mode global|window` kept as deprecated aliases.) |
| **Conda env** | **`kmate`** (mamba; `/global/home/users/tbellg/miniforge3/envs/kmate`) — the only env. `hapfm`/`gwas`/`sequencing_pipeline`/`pang`/`pangenie`/`basic` are gone (cluster migration). |

**ARCHIVED — DO NOT CONSUME (these are the recurring confusion; they live under `archive/`):**
- `founders_231_v3qc*.vcf.gz` (the `bcftools norm` naive merge, any v3qc/v3qc_v2/v3qc_v3) — superseded by the arch3 decomposition. arch3 reproduces the same v3qc_v3 genotype QC (per-cell het→`.`, no V4 record-drop) via its A2 step; it just decomposes by symbolic-ID instead of `norm`, which avoids the carrier-loss bug.
- `data/var_pa_231_v3qc_v3.*` and `data/kmer_pa_231_v3qc_v3[_filt2]/` — old matrices off the naive merge.
- The simple `build_var_pa_v3qc_v3.sh` path — V_pa comes from arch3 A5, never directly off the haploid VCF.

**Why arch3 only:** `bcftools norm -m -any` scatters/drops carriers at co-located
multi-allelic sites (up to ~99% carrier loss at one measured SNP). arch3's
symbolic-ID propagation handles convergence/position-shift by construction. See
`docs/INVESTIGATION_CN_VAR_DECOMPOSITION.md`.

**Panel-conditional caveat (not for the paper):** on the homogeneous p80 control,
$\omega_k=1/m_b$ slightly under-performs unweighted EM; `--kmer-weight uniform`
reproduces the MLE byte-identically for balanced-panel users. See `ALGORITHM.md` §10 M6.

---

## §0.1 — Run recipe (per sample)

`--unit chrom` (default) is the production selfing path; `--unit bp` is the
recombinant path (the legacy anchored+smoothed `★★`/star2 recipe — window-bp
10000, anchor 0.3, hmm passes 5 α 0.5 — is now opt-in behind `--no-local-only`).
Env: `kmate`. Replace `Chr1` / chrom paths as needed (loop Chr1..Chr5).

```bash
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
CHR=Chr1
$PY src/per_sample_per_chrom.py \
    --kmer-pa-prefix data/kmer_pa_231_arch3_filt2inv/kmer_pa \
    --var-pa     panel/arch3/${CHR,,}/var_pa_231_arch3_${CHR,,}.var_pa.npz \
    --var-called panel/arch3/${CHR,,}/var_pa_231_arch3_${CHR,,}.var_called.npz \
    --var-meta   panel/arch3/${CHR,,}/var_pa_231_arch3_${CHR,,}.meta.npz \
    --reads <r1.fq> <r2.fq> --sample <name> --out <name>.tsv \
    --threads 8 --chroms $CHR --normalize per_founder --kmer-weight uniform \
    --unit chrom              # DEFAULT (production selfing): per_founder + uniform.
                              # Recombinant pools: --unit bp --window-bp 10000 --kmer-weight inv_mb
```

Output TSV (8 cols): `chrom  pos  ref_len  alt_len  alt_freq  info  n_called  se`.
Also writes `*.h_per_chrom.npz` (global) or `*.h_blocks_per_chrom.npz` (window).

## §0.2 — Reproducibility / environment

- **Conda env**: `kmate` (mamba; `/global/home/users/tbellg/miniforge3/envs/kmate`). The
  only env. Created on the current cluster after the migration removed `hapfm`/`gwas`/
  `sequencing_pipeline`/`pang`/`pangenie`. Contains numpy/scipy/pysam/pandas +
  samtools/bcftools/htslib/jellyfish/kmer-jellyfish + bbmap/minimap2/wgsim.
- **Reference FASTA**: `TAIR10.chr.iupacN.fa` (cactus normalizes IUPAC→N; `build_kmer_pa`
  uses it for bubble flanks only).
- **Chrom naming**: GrENE-Net VCF uses `1..5`; cactus / kMate / arch3 use `Chr1..Chr5`.
- **Sample list**: the 231-founder VCF header has a blank trailing 232nd column; drop with
  `bcftools query -l … | grep -v '^$'`.
- **HARP binary** (hapFIRE comparator only): `external/HapFIRE/.../bin/harp`.

## TL;DR (historical framing; §0 is authoritative)

For each panel record, estimate per-record ALT allele frequency from pool-seq reads via a per-sample Poisson EM on the 231-founder simplex (M-step normalization `per_founder`, the default; window-mode weighting $\omega_k = 1/m_b$, global-mode weighting uniform — superseded 2026-07-06 for global mode, see `docs/FOUNDER_NORMALIZATION_FIX.md`), then project through `var_pa` to per-record AF. Production pipeline as of today:

- **Panel**: 231 founders (78 cactus + 153 PG), decomposed by **arch3** → `merged_231_chr{N}_final.vcf.gz` (the only panel VCF; §0)
- **Decomposition**: **arch3 ONLY** (annotate_vcf + convert-to-biallelic), NOT `bcftools norm -m -any`
- **Matrices**: K_pa `kmer_pa_231_arch3_filt2inv` + V_pa `var_pa_231_arch3` (SV-level) AND `var_pa_231_arch3_atomized` (per-base, SNP-level GEA) — all from merged_231 (§0/§1)
- **K-mer filter / kmer_pa build**: `filt2inv` (drop ac=1 + ac=F, inline `--filter-production`). EM weighting: window-mode $\omega_k=1/m_b$ (`--kmer-weight inv_mb`); global-mode (production) `--kmer-weight uniform` + `--normalize per_founder` (superseded 2026-07-06 for global mode — see §0 and `docs/FOUNDER_NORMALIZATION_FIX.md`).
- **Projection**: MAR — `(h @ var_pa) / (h @ var_called)` in both global and window modes. Patched 2026-05-21.
- **Haploblock collapse (2026-07-07, ON by default both modes)**: each unit (window, or whole chromosome in global) is fit `unit → distinct haplotypes (K_b ≤ F) → EM → equal-split back to members`, instead of assuming all F=231 founders are separately identifiable (`haploblock_collapse=True`; `em_solver.haploblock_collapse_indices`). `--haploblock-eps` sets the merge tolerance (default 0 = exact k-mer-identical → **exact byte-identical no-op when K_b==F**, which is essentially always true on the arch3 panel at whole-chromosome / r²=0.1+ units, so the GLOBAL production path is numerically inert; it matters for window mode where a 10kb window may carry only a few haplotypes). `--emit-af-se` and the bootstrap compute uncertainty on the K_b classes (mapped to founders by the equal-split), not the singular full-F Fisher info. See `ALGORITHM.md` §4.4.
- **Estimation unit (`--unit`, 2026-07-07)**: one estimator; the mode is the unit, each fit locally (haploblock collapse → EM → project; no anchor / smoothing / fallback). **`--unit chrom`** (the CLI default) = one h per chrom (former `--block-mode global`; byte-identical alias retained; supports `--h-only` + `--emit-af-se`) is the **production estimator** for the heavy-selfing GrENE-Net cohort — robust on uniform and sparse panels. `--unit ld --ld-r2 0.1` (r²-LD CompleteLDPartition blocks from var_pa; `src/kmate/ld_partition.py`) is a per-block option that **collapses in low-diversity blocks** (centromere) → **wrong for selfing pools**; it is NOT the corrected production estimator. `--unit bp` (fixed windows) / `--unit tsv` (explicit blocks) for recombinant pools. `--block-mode global|window` kept as deprecated aliases. The old fixed-bp *window* cohort outputs (`grenenet_kmate_window*`) are **STALE** (old multinomial EM; decision 2026-07-07) — archived (~1.6 TB) to `results/archive/`; see `docs/RERUN_AFTER_FIX.md` Group D. Anchored+smoothed "star2" recipe is now behind `--no-local-only`.
- **Naming**: the method is **kMate** (algorithm + math: `ALGORITHM.md`). Legacy code and result-dir paths (`*/cactus_em_*`) still use the prior name `cactus_em`; with the k-mer-filter decision now closed, the rename is unblocked but not yet executed. (The former `sims/visor_freqk` sub-repo was absorbed into `sims/` 2026-05-30; see `sims/README.md`.)

## What changed since the last pipeline-state doc (2026-05-19)

| Change | Why |
|---|---|
| **Panel filtered to segregating-only (2026-06-02)** | A4 Step-3 filter widened from `AN=0`-only to `AC=0 \|\| AC=AN`. Dropped 1,835,718 monomorphic records (17.8%); VCF + V_pa overwritten in place, 10,325,364 → 8,489,646. Estimates on kept records unchanged (per-record projection). See `results/panel_stats/PANEL_STATS.md`, `scripts/apply_monomorphic_filter.py`. |
| Arch decomposition replaces `norm -m -any` for var_pa | annotate_vcf + convert-to-biallelic produces symbolic-ID biallelic catalog with +17pp hapFIRE-SNP coverage gain (`project_arch3_chr1_validation` memory) |
| `var_pa_atomized` added as a production deliverable | Per-base SNP catalog with carriers UNIONed across overlapping records. −23% RMSE, −43% outliers for SNP-level GEA (`project_arch3_atomization_result`) |
| MAR projection in window mode | Previously window mode silently treated `.` as REF (different from global mode's called-mask normalization). Inconsistent semantics fixed 2026-05-21. Star2 numbers prior to this patch are stale. |
| K-mer filter declared **undecided** | Prior doc incorrectly said mixed-loose was production. mixed-loose was evaluated, didn't win unanimously. Active testing: `filt2`. **(Superseded 2026-05-27: filt2 + ω_k=1/m_b chosen; see §0.)** |

## 1. Three persistent matrices feed the pipeline

All three are built from the arch3 canonical VCF `merged_231_chr{N}_final.vcf.gz` (§0).

| Matrix | Path | Purpose | Build script |
|---|---|---|---|
| `kmer_pa` (K_pa) | **`data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr{N}.{kmer_pa,meta}.npz`** (filt2inv) | founder × k-mer (k=31, in-house `ours_` index + merged_231 VCF) | `scripts/build_kmer_pa_production_arch3.sh` |
| `var_pa` (V_pa) | `panel/arch3/chr{N}/var_pa_231_arch3_chr{N}.{var_pa,var_called,meta}.npz` | founder × biallelic-record carriers + called-mask | `panel/arch3/chr{N}/jobA5_build_cnvar.sh` |
| `var_pa_atomized` | `panel/arch3/chr{N}/var_pa_231_arch3_chr{N}_atomized.{var_pa,var_called,meta}.npz` | founder × per-base SNP (atomized) | `panel/arch3/chr{N}/jobD1_atomize_cnvar.sh` |

The k-mer index is the **in-house `pang_135_haploid/ours_Chr{N}` builder** (`panel/pangenie_index`), Level-A+B validated equivalent to PanGenie's index. Built once on the 135-assembly pangenome graph; reused (it supplies bubble windows + candidate k-mers, the merged_231 VCF supplies founder carriers).

## 2. Driver — `src/per_sample_per_chrom.py`

Reads FASTQs, runs k-mer Poisson EM on the simplex, projects through var_pa. Output TSV schema (post 2026-05-21 patch):

```
chrom  pos  ref_len  alt_len  alt_freq  info  n_called  se
```

- `alt_freq` = `(h @ var_pa)[r] / (h @ var_called)[r]` — MAR projection
- `info` = h-weighted called mass at record r (h-dependent, varies with method)
- `n_called` = integer count of called founders at record r (h-independent panel QC)
- `se` = Wald SE using `n_called` as effective sample size

All units (`--unit chrom` default; `ld`/`bp`/`tsv`) use the same projection semantics. The legacy anchored+smoothed `★★`/star2 window recipe (window-bp 10000, global-anchor-weight 0.3, hmm-smooth-passes 5, hmm-smooth-alpha 0.5) is now opt-in behind `--no-local-only`; the default `--unit bp` window is fit local-only. **EM weighting is unit-specific:** the recombinant `--unit bp` path uses `--kmer-weight inv_mb`, the production `--unit chrom` (selfing) path uses `--kmer-weight uniform` + `--normalize per_founder` (superseded 2026-07-06 — see §0 and `ALGORITHM.md` §4.3). The LD-block modes, overlapping windows, and the older k-mer-budget rebalancing / carrier-weighting / contamination-ω variants were archived to `src/archive/`; the authoritative file list + recipes are in `src/README.md`.

## 3. Outstanding production work

| Task | Why |
|---|---|
| **Validate arch3 Chr2–5 → merged_231_chr{N}_final.vcf.gz** | All 5 chroms now built on disk (2026-05-30/31). Confirm the Chr2–5 decompositions match the Chr1-validated quality before relying on them for downstream GEA. |
| **Verify K_pa `kmer_pa_231_arch3_filt2inv` + V_pa `var_pa_231_arch3` for all 5 chroms** | All 5 chroms of both matrices are present on disk (off the merged_231 VCFs above), replacing every v3qc-named matrix. Sanity-check counts/coverage before the cohort run. |
| Re-run SEEDMIX / evolved baselines under the full production recipe | Any TSV predating arch3 var_pa + MAR projection is stale — and (2026-07-06) any predating the global-mode switch to `--normalize per_founder` + `--kmer-weight uniform` must be re-run under it (see §0). |
| Production scale-out on ~2,415 evolved GrENE-Net samples | SLURM template at `grenenet/run_site_array_perchrom.sh`; ~1.5–5 days at cluster-wide concurrency |

## 4. What's deprecated / archived (DO NOT USE — see §0)

- **`founders_231_v3qc*.vcf.gz` (any v3qc / v3qc_v2 / v3qc_v3 naive `norm` merge)** — superseded by the arch3 decomposition. Moved to `archive/`.
- **`data/var_pa_231_v3qc_v3.*` and `data/kmer_pa_231_v3qc_v3[_filt2]/`** — old matrices off the naive merge. Moved to `archive/`.
- `kmer_pa_231_v3qc_v3_mixedloose/`, mixed-loose / mixed-conserv / mixed-strict filters — evaluated, not production.
- "Star2 treats `.` as REF" projection — buggy, patched 2026-05-21. Any star2 result TSV without `info`/`n_called`/`se` columns is from the buggy code.
- Beagle imputation for SVs — hard-rejected (−25 to −30pp concordance loss).

## 5. Companion docs (still current)

- `BACKGROUND.md` — project framing
- `ALGORITHM.md` — kMate algorithm, math & production wiring (code-verified single source of truth; supersedes the old prose doc and `CACTUS_EM_MATH.md`)
- `old_docs/CACTUS_EM_MATH.md` — superseded formal-math doc (folded into `ALGORITHM.md`)
- `SIMULATIONS_METHODS.md` — methods-ready description of the pool-seq simulation framework
- `INVESTIGATION_CN_VAR_DECOMPOSITION.md` — context for the arch decomposition switch
- `NOCAP_INDEX_INVESTIGATION.md` — why PanGenie caps per-bubble k-mers, and why kMate's `inv_mb` weight makes them unnecessary in window mode (verdict: keep caps; global-mode production now uses `uniform` + per_founder, where the caps do matter — see that doc's 2026-07-06 note; no-caps untested on accuracy)
- `MISSINGNESS_231PANEL.md` — F_MISSING characterization on the production panel
- `PIPELINE_FASTQ_PREPROCESSING.md` — read-side preprocessing pipeline (trim, dedup)
- `archive/exploration/panel_overlap_135_vs_82/RESULTS.md` — panel composition analysis
- `RESULTS_LOG.md` — chronological numerical record (large file; useful historical reference, not authoritative for current state)

## 6. Subprojects (own subfolder READMEs)

- `grenenet/` — **cohort scale-out** of kMate across the ~2,415 evolved GrENE-Net pool-seq samples (SLURM). The runner, count-once DB reuse, and the shared-filesystem I/O fix (`JF_DIR` → node-local/RAM for the transient per-sample DB) are documented in `grenenet/README.md`. Production submit = single-phase `run_site_array_perchrom.sh` (`BLOCK_MODE=global`, `--mem=32G`, manifest `data/sample_manifest_usesample.tsv`).
- `benchmarks/p80/` — homogeneous 80-cactus-founder control experiment. See `benchmarks/p80/README.md`.
- `sims/` — pool-seq simulation framework (self-contained as of 2026-05-30; the former `visor_freqk` sub-repo was absorbed here). See `sims/README.md` and `docs/SIMULATIONS_METHODS.md`.

## 7. Historical record

All docs about prior pipeline states / superseded decisions are preserved in `old_docs/`. Useful for archaeology; do not treat as authoritative.

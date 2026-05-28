# Session handoff — kMate

**Last updated:** 2026-05-27
**Project root:** `/global/scratch/users/tbellg/kmate/`

## TL;DR

`kMate` end-to-end: per-sample **weighted** k-mer Poisson EM on the 231-founder simplex (production weight $\omega_k = 1/m_b$, per-bubble de-replication; see `ALGORITHM.md` §4.2), projected through `cn_var` to per-record AF (SNPs + indels + SVs in one pass). Production panel uses the **arch decomposition** (annotate_vcf + convert-to-biallelic). MAR-aware projection is the recipe in both `global` and `★★` window modes. **K-mer filter resolved 2026-05-27: `cn_full_231_v3qc_v3_filt2` (drop ac=1 singletons) + EM weighting $\omega_k = 1/m_b$.** See `docs/METHODS_TRIED_AND_RESULTS.md` §0/§3 for the full sweep history and the panel-conditional caveat (1/m_b is opt-in via `--kmer-weight {uniform,inv_mb}` for users on balanced panels).

**Naming:** the method is **kMate** (see `ALGORITHM.md`). Legacy code, result-dir paths (`*/cactus_em_*`), and the `sims/visor_freqk` sub-repo still carry the prior name `cactus_em`; with the k-mer-filter decision now closed (2026-05-27), the path/code rename is unblocked but not yet executed.

## Authoritative source on current state

**`docs/PIPELINE_STATE.md`** is the single source of truth for what's production vs in-evaluation vs deprecated. Read it before making decisions about the panel, projection, or filters.

## Production recipe (2026-05-27)

| component | choice |
|---|---|
| Panel | `founders_231_v3qc_v3` (78 cactus + 153 PG) |
| Decomposition | **arch3** (annotate_vcf + convert-to-biallelic), NOT `bcftools norm -m -any` |
| cn_var (SV-level) | `panel/arch3/chr1/cn_var_231_arch3_chr1.{cn_var,cn_var_called,meta}.npz` |
| cn_var (SNP-level) | `panel/arch3/chr1/cn_var_231_arch3_chr1_atomized.*` (per-base atomized) |
| cn_full | **`data/cn_full_231_v3qc_v3_filt2/cn_Chr1.{cn,meta}.npz`** (filt2: drop ac=1 singletons) |
| EM weighting | **`--kmer-weight inv_mb`** (ω_k = 1/m_b per-bubble de-replication, `ALGORITHM.md` §4.2) |
| Projection | MAR: `(h @ cn_var) / (h @ cn_var_called)`, both `global` and window modes |
| Chrom scope | **Chr1 only currently — Chr2–5 build is the open production task** |

Two output modes, both production-supported:

```bash
# global — default for SEEDMIX / F0 pools
python src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3qc_v3_filt2/cn \
    --cn-var       panel/arch3/chr1/cn_var_231_arch3_chr1.cn_var.npz \
    --cn-var-called panel/arch3/chr1/cn_var_231_arch3_chr1.cn_var_called.npz \
    --cn-var-meta  panel/arch3/chr1/cn_var_231_arch3_chr1.meta.npz \
    --reads <r1.fq> <r2.fq> --sample <name> --out <out.tsv> \
    --threads 8 --chroms Chr1 \
    --block-mode global \
    --kmer-weight inv_mb

# window ("★★") — for high-recomb regimes (evolved pools with multi-gen mosaic
# ancestry). The window-mode DEFAULTS are this recipe (window-bp 10000,
# global-anchor-weight 0.3, hmm-smooth-passes 5, hmm-smooth-alpha 0.5), so the
# flag alone reproduces it — pass those flags only to override.
python src/per_sample_per_chrom.py \
    [same inputs as global] \
    --block-mode window \
    --kmer-weight inv_mb
```

**Estimator code (2026-05-26 cleanup):** two modes only — `global` and `window`.
LD-block modes, overlapping windows, the older k-mer-budget rebalancing
(`--row-normalize-cn`), carrier-weighting and contamination-ω were archived to
`src/archive/`. The current production EM weighting is the cleaner
$\omega_k = 1/m_b$ composite-likelihood form, still wired into the active
driver (`--kmer-weight inv_mb`); see `ALGORITHM.md` §4.2. The authoritative
file list + invocation recipes are in **`src/README.md`**.

Output TSV (post-2026-05-21 patch) has 8 columns:

```
chrom  pos  ref_len  alt_len  alt_freq  info  n_called  se
```

`info` (h-weighted observed mass), `n_called` (h-independent panel count), `se` (Wald SE) are new per-record uncertainty metrics — see `docs/PIPELINE_STATE.md` §2.

## What's pending

1. **Arch 3 Chr2–5 panel build** — run A1→A5 for remaining chroms. Chr1 is validated; whole-genome needed for downstream GEA.
2. ~~**Choose production k-mer filter**~~ — **RESOLVED 2026-05-27**: `filt2` (drop ac=1) + EM weighting $\omega_k=1/m_b$ (`--kmer-weight inv_mb`). See `docs/METHODS_TRIED_AND_RESULTS.md` §0/§3.
3. **Re-validate SEEDMIX baselines under MAR + arch cn_var + production weighting** — prior numbers used `bcftools norm -m -any` cn_var, the (now-patched) "star2 treats `.` as REF" projection, AND unweighted EM. All star2 result TSVs without `info`/`n_called`/`se` columns are stale, as are all results that predate the `--kmer-weight inv_mb` switch.
4. **Production scale-out on ~2,500 evolved GrENE-Net samples** — SLURM template at `tests/run_site_array_perchrom.sh`. Blocked on (1).
5. **Subprojects**: `control_p80/` (homogeneous 80-cactus-founder control, all 6 regimes done; established the $\omega_k=1/m_b$ panel-conditional caveat — see `control_p80/results/FINAL_RESULTS_cov10_p80.ipynb`).

## Companion docs (still current)

| File | Purpose |
|---|---|
| `docs/PIPELINE_STATE.md` | Production-state SoT |
| `BACKGROUND.md` | Project framing |
| `ALGORITHM.md` | kMate algorithm, math & wiring (code-verified single source of truth) |
| `src/README.md` | Estimator source inventory — active files, two recipes, what was archived (2026-05-26) |
| `old_docs/CACTUS_EM_MATH.md` | Formal math (superseded; folded into `ALGORITHM.md`) |
| `docs/INVESTIGATION_CN_VAR_DECOMPOSITION.md` | Why we switched to arch decomposition |
| `docs/MISSINGNESS_231PANEL.md` | F_MISSING characterization on the production panel |
| `docs/PIPELINE_FASTQ_PREPROCESSING.md` | Read-side preprocessing pipeline |
| `docs/SIMULATIONS_METHODS.md` | Methods-ready description of the pool-seq simulation framework (regime matrix, parameters, citations) |
| `panel_overlap_135_vs_82/RESULTS.md` | Panel composition analysis |
| `docs/RESULTS_LOG.md` | Chronological numerical record (large file; historical reference, not authoritative) |
| `data/exclude_list.txt` | Assembly_IDs dropped from cactus panel (101003 + 100852) |
| `data/flag_list.tsv` | Per-Assembly_ID flag status |
| `sims/visor_freqk/README.md`, `RECOMB_SIM.md` | Pool-seq sim framework |
| `control_p80/README.md` | 80-founder homogeneous control |

Historical / superseded docs are preserved under `old_docs/`. Useful for archaeology; do not treat as authoritative.

## Reproducibility

- **Conda env**: `hapfm` (`/home/tbellagio/miniforge3/envs/hapfm/`)
- **HARP binary**: `/carnegie/nobackup/scratch/xwu/haplotype_frequency_estimation/hapFIRE_sourcecode/bin/harp` (for hapFIRE comparator runs only)
- **GrENE-Net VCF sample list**: VCF header has 232 cols; 232nd is blank trailing. Drop with `bcftools query -l … | grep -v "^$"` to get the 231.
- **Chrom naming**: GrENE-Net VCF uses `1..5`; cactus / kMate use `Chr1..Chr5`. Conversion handled per-pipeline.
- **Reference FASTA**: cactus VCFs normalize IUPAC codes to N. Use `TAIR10.chr.iupacN.fa` for `bcftools consensus` and downstream FASTA work.

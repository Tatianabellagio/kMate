# Session handoff — hapFIRE-SV

**Last updated:** 2026-05-22
**Project root:** `/carnegie/nobackup/scratch/tbellagio/hapfire_sv/`

## TL;DR

`cactus_em` end-to-end: per-sample k-mer Poisson EM on the 231-founder simplex, projected through `cn_var` to per-record AF (SNPs + indels + SVs in one pass). Production panel uses the **arch decomposition** (annotate_vcf + convert-to-biallelic). MAR-aware projection is now the recipe in both `global` and `★★` window modes. K-mer filter is **undecided** (testing `filt2`; `mixed-loose` is NOT production despite earlier docs claiming otherwise).

## Authoritative source on current state

**`PIPELINE_STATE_2026-05-22.md`** is the single source of truth for what's production vs in-evaluation vs deprecated. Read it before making decisions about the panel, projection, or filters.

## Production recipe (2026-05-22)

| component | choice |
|---|---|
| Panel | `founders_231_v3qc_v3` (78 cactus + 153 PG) |
| Decomposition | **arch3** (annotate_vcf + convert-to-biallelic), NOT `bcftools norm -m -any` |
| cn_var (SV-level) | `arch3/chr1/cn_var_231_arch3_chr1.{cn_var,cn_var_called,meta}.npz` |
| cn_var (SNP-level) | `arch3/chr1/cn_var_231_arch3_chr1_atomized.*` (per-base atomized) |
| cn_full | `poolfreq/data/cn_full_231_v3qc_v3/cn_Chr1.{cn,meta}.npz` (k-mer filter still TBD) |
| Projection | MAR: `(h @ cn_var) / (h @ cn_var_called)`, both `global` and window modes |
| Chrom scope | **Chr1 only currently — Chr2–5 build is the open production task** |

Two output modes, both production-supported:

```bash
# global — default for SEEDMIX / F0 pools
python poolfreq/src/per_sample_per_chrom.py \
    --cn-kmer-prefix poolfreq/data/cn_full_231_v3qc_v3/cn \
    --cn-var       arch3/chr1/cn_var_231_arch3_chr1.cn_var.npz \
    --cn-var-called arch3/chr1/cn_var_231_arch3_chr1.cn_var_called.npz \
    --cn-var-meta  arch3/chr1/cn_var_231_arch3_chr1.meta.npz \
    --reads <r1.fq> <r2.fq> --sample <name> --out <out.tsv> \
    --threads 8 --chroms Chr1 \
    --block-mode global

# ★★ — for high-recomb regimes (evolved pools with multi-gen mosaic ancestry)
python poolfreq/src/per_sample_per_chrom.py \
    [same as above] \
    --block-mode window --window-bp 10000 \
    --global-anchor-weight 0.3 \
    --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5
```

Output TSV (post-2026-05-21 patch) has 8 columns:

```
chrom  pos  ref_len  alt_len  alt_freq  info  n_called  se
```

`info` (h-weighted observed mass), `n_called` (h-independent panel count), `se` (Wald SE) are new per-record uncertainty metrics — see `PIPELINE_STATE_2026-05-22.md` §2.

## What's pending

1. **Arch 3 Chr2–5 panel build** — run A1→A5 for remaining chroms. Chr1 is validated; whole-genome needed for downstream GEA.
2. **Choose production k-mer filter** — testing `filt2` (AC≥2). `mixed-loose` is NOT production despite earlier docs. Other candidates: raw cn_full (no filter), `rownorm` (rejected per memory).
3. **Re-validate SEEDMIX baselines under MAR + arch cn_var** — prior numbers used `bcftools norm -m -any` cn_var AND the (now-patched) "star2 treats `.` as REF" projection. All star2 result TSVs without `info`/`n_called`/`se` columns are stale.
4. **Production scale-out on ~2,500 evolved GrENE-Net samples** — SLURM template at `poolfreq/tests/run_site_array_perchrom.sh`. Blocked on (1).
5. **Subprojects**: `control_p80/` (homogeneous 80-cactus-founder control, all 6 regimes done — see `control_p80/results/FINAL_RESULTS_cov10_p80.ipynb`).

## Companion docs (still current)

| File | Purpose |
|---|---|
| `PIPELINE_STATE_2026-05-22.md` | Production-state SoT |
| `BACKGROUND.md` | Project framing |
| `ALGORITHM.md` | cactus_em prose walkthrough |
| `CACTUS_EM_MATH.md` | Formal math |
| `INVESTIGATION_2026-05-19_CN_VAR_DECOMPOSITION.md` | Why we switched to arch decomposition |
| `MISSINGNESS_231PANEL.md` | F_MISSING characterization on the production panel |
| `PIPELINE_FASTQ_PREPROCESSING.md` | Read-side preprocessing pipeline |
| `SIMULATIONS_METHODS.md` | Methods-ready description of the pool-seq simulation framework (regime matrix, parameters, citations) |
| `panel_overlap_135_vs_82/RESULTS.md` | Panel composition analysis |
| `RESULTS_LOG.md` | Chronological numerical record (large file; historical reference, not authoritative) |
| `data/exclude_list.txt` | Assembly_IDs dropped from cactus panel (101003 + 100852) |
| `data/flag_list.tsv` | Per-Assembly_ID flag status |
| `sims/visor_freqk/README.md`, `RECOMB_SIM.md` | Pool-seq sim framework |
| `control_p80/README.md` | 80-founder homogeneous control |

Historical / superseded docs are preserved under `old_docs/`. Useful for archaeology; do not treat as authoritative.

## Reproducibility

- **Conda env**: `hapfm` (`/home/tbellagio/miniforge3/envs/hapfm/`)
- **HARP binary**: `/carnegie/nobackup/scratch/xwu/haplotype_frequency_estimation/hapFIRE_sourcecode/bin/harp` (for hapFIRE comparator runs only)
- **GrENE-Net VCF sample list**: VCF header has 232 cols; 232nd is blank trailing. Drop with `bcftools query -l … | grep -v "^$"` to get the 231.
- **Chrom naming**: GrENE-Net VCF uses `1..5`; cactus / cactus_em use `Chr1..Chr5`. Conversion handled per-pipeline.
- **Reference FASTA**: cactus VCFs normalize IUPAC codes to N. Use `TAIR10.chr.iupacN.fa` for `bcftools consensus` and downstream FASTA work.

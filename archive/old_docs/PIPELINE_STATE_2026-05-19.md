# Pipeline state — global cactus_em + 231-founder panel (2026-05-19)

Snapshot after the v3qc_v3 work (het-mask + bcftools merge bug fix + cn_var rebuild). This documents the cactus_em pipeline state as of 2026-05-19, which forms the foundation for **Arch 3** (production as of 2026-05-21). The cn_full (v3qc_v3 mixedloose) and EM layer are unchanged in Arch 3; what changed is the cn_var panel build (graph-annotated, symbolic-ID decomposition — see `arch3/chr1/jobA1…A7` scripts). Companion docs: `BACKGROUND.md` (project context), `ALGORITHM.md` (math), `HANDOFF.md` (current production pointer), `RESULTS_LOG.md` (numbers over time).

---

## 1. Global cactus_em — what it computes

For one pool-seq sample (pair of FASTQs), produce per-variant ALT allele frequency:

```
reads (FASTQ pair)
    │
    │  Jellyfish k-mer counting
    ▼
counts[k] : K-vector (one per k-mer in pangenome dict)
    │
    │  Multiplicative-update EM (Poisson mixture)
    │   solve_em(counts, cn_full, coverage, max_iter=200, tol=1e-7)
    ▼
h : F-vector  (founder mass fractions, simplex)
    │
    │  AF projection through cn_var with the called-mask
    │   af[k] = (h @ cn_var[:, k]) / (h @ cn_var_called[:, k])
    ▼
alt_freq : per-record AF (TSV chrom, pos, ref_len, alt_len, alt_freq)
```

**Three persistent matrices feed the pipeline**, all built once per panel build:

| Matrix | Shape | What it encodes | Built by |
|---|---|---|---|
| **cn_full** | `(F, K_kmers)` binary CSR | does founder f carry k-mer k? (k-mer = 31bp from PanGenie's bubble dictionary) | `build_kmer_cn.py` using PanGenie-index k-mers + haplotype FASTAs |
| **cn_var** | `(F, N_records)` binary CSR | does founder f carry ALT at biallelic-decomposed record k? | `build_cn_var.py` reading the haploid VCF |
| **cn_var_called** | `(F, N_records)` binary CSR | is founder f's GT *non-missing* at record k? | `build_cn_var.py`, NEW in v3qc_v3 — fixes the AC/AN denominator bug |

`F = 231` (78 cactus + 153 PG). The k-mer dictionary is the **pang_135** PanGenie-index (built once on the 135-assembly minigraph-cactus pangenome graph, not re-built per panel — see `SESSION_2026-05-17.md`).

**Driver**: `poolfreq/src/per_sample_per_chrom.py` (one chrom at a time for memory). Key flags:
- `--row-normalize-cn` + `--kf-correction-alpha=N` — the α-sweep exit path (now empirically inferior to mixed-loose; kept for diagnostics)
- `--cn-var-called` — auto-detected next to `--cn-var`; enables the AC/AN projection
- `--block-mode global` is the production mode for SEEDMIX / pool-seq AF; window/LD modes are alternative routes parked for high-recomb regimes

---

## 2. Panel evolution: v3 → v3qc → v3qc_v2 → v3qc_v3

We've gone through three named "v3qc" rebuilds, each fixing a specific bug discovered downstream. The history matters because **each rebuild was triggered by an observed AF artifact** that traced to a panel-construction step.

### v3 (baseline, 2026-05-12)

- 231 founders = 80 cactus + 151 PG
- No PG GT filter beyond PanGenie's defaults; no GQ mask, no V4 filter
- Multi-allelic records left intact in cactus VCF; PG decomposed inconsistently
- Result: production cactus_em-v3 with ~190× cactus-PG h-bias, ~5% AF MAE vs hapFIRE
- Locked in via memory `hapFIRE-SV project context`; baseline for all comparisons

### v3qc (2026-05-16) — first QC layer

Triggered by xwu's GrENE-Net comparison showing systematic AF inflation on cactus side.

Changes vs v3:
- Drop 2 cactus assemblies (101003 Set-1, 100852 Ped-0) — see `exclude_list.txt`
- Add 7 PG-substituted founders genotyped from 1001G short reads (5772, 9947 = the two excluded cactus, plus 5 with low LOO concordance kept as cactus but flagged)
- Apply **GQ ≥ 20** mask + **V4 < 0.01** (xwu's excess-heterozygosity filter) on PG side
- Decompose multi-allelics with `bcftools norm -m -any` BEFORE filtering (multi-allelic OR-semantics bug — see `SESSION_2026-05-17.md`)
- AC=0 cleanup on merged 231 panel

Result: ~30× h-bias, ~3% AF MAE. Still cactus-heavy.

### v3qc_v2 (2026-05-17/18) — cn_var/cn_full bug fixes

Triggered by k-mer composition analysis showing cactus founders carried ~2× more total k-mers than PG founders (asymmetry A).

Changes vs v3qc:
- **`--treat-missing-as-n` flag in `build_kmer_cn.py`**: at `./.` positions, write `N` into the reconstructed haplotype FASTA (not REF). This stops `./.` cells from being silently treated as REF carriers in cn_full.
- **cn_var meta gains REF/ALT base strings** (uses `dtype=object` to handle 100kb+ SV alleles; otherwise numpy tries to allocate 2.3 TB for a fixed-width string array)
- **cn_var_called sparse matrix added** — tracks which (founder, record) cells are non-missing. AF projection becomes `(h @ cn_var) / (h @ cn_var_called)` instead of `h @ cn_var / 1`. Without this, missing PG founders silently inflate cactus AF at cactus-only records.
- Multi-allelic alignment fix at the comparison layer: join hapFIRE by `(chrom, pos, REF, ALT)` instead of position only.

Result: cactus_em raw vs hapFIRE MAE = 0.023 with the FIXED cn_var, but residual 27.5% of records still have F_MISSING > 0.5 (the high-missingness chunk).

### v3qc_v3 (2026-05-18/19) — foundation for Arch 3 (current production)

Triggered by observation that the F_MISSING > 0.5 chunk was structural — caused by V4 dropping entire PG records (~2M genome-wide) and leaving cactus-only after merge.

Changes vs v3qc_v2:
- **Het-mask replaces V4 record drop** (`pangenie_genotyping/scripts/build_v3qc_v3_phase_a.sh`):
  - `bcftools +setGT -- -t q -i 'GT="het"' -n .` converts every diploid `0/1` cell to `./.`, then `fill-tags` recomputes AC/AN. The record is kept (V4 was dropping the whole row).
  - Justification: 1001G is inbred; `0/1` calls are short-read genotyping artifacts. Keeping the `0/0` and `1/1` cells preserves real carrier information.
- **Pre-haploidize PG before merge** (`haploidize_pg_hetmasked.sh`):
  - `bcftools merge` of haploid cactus + diploid PG has a context-dependent bug that drops PG GTs at AC=0 records on large datasets (~7M records). The isolated small-region merge works correctly; the full merge fails. Pre-haploidizing PG (so both inputs are haploid) eliminates the bug.
  - Verified empirically on 11 test positions in the broken merge — all preserved with the fix.

Result: cactus_em mixed-loose vs hapFIRE MAE = 0.015, 0.88% outliers, F_MISSING > 0.5 chunk collapsed from 27.5% → 0.49%.

---

## 3. The mixed-loose data-level rebalancing (the key methodological win)

Independent of the panel bug fixes, we tried multiple algorithmic mitigations for the residual cactus-PG h-bias:

| Approach | Idea | Result |
|---|---|---|
| `--row-normalize-cn` | Each cn_full row divided by K_f (per-founder mass=1) | 8.5× h-bias, 0.039 AF MAE |
| `--kf-correction-alpha α` | Post-EM rescale `h ∝ h_em / K_f^α` | α=3 gives h≈balanced, MAE 0.022; but per-founder R² with hapFIRE drops to 0 (rescaling averages, not founder identity) |
| `filt2` (drop ac<2 k-mer columns) | Eliminate per-founder uniqueness asymmetry | Worse on every metric — removes the signal we need |
| **mixed-loose**: filter cn_full to k-mers carried by ≥1 cactus AND ≥1 PG | Use only k-mers reliably represented on both panel sides — class-agnostic, conceptually parallel to hapFIRE's SNP-only strategy | **Best on every metric**: AF MAE 0.015, h-bias near 1, per-founder R² 0.23 with hapFIRE |

mixed-strict (`ac ≥ 5 each side`) and mixed-conserv (`ac ≥ 10 each side`) over-correct (cactus slightly under-weighted) and lose AF accuracy. **mixed-loose is the production cn_full filter** (used in Arch 3).

This was the conceptual breakthrough: rather than penalize the EM with post-hoc corrections, **drop the k-mers that aren't usable evidence in both panels.** The remaining ~9.8M k-mers per Chr1 (43% of cn_full) carry full discriminating power because both panels confirm them.

---

## 4. v3qc_v3 baseline numbers (SEEDMIX_S1 Chr1) — superseded by Arch 3

> **Arch 3 (current production)** numbers: MAE 0.0131, RMSE 0.0225, |Δ|>0.10 = 0.27%, hapFIRE SNP coverage 79.8% — see `AF_TRUTH_VS_ESTIMATE_arch3_chr1.ipynb`.

Using `cn_full_v3qc_v3_mixedloose` + `cn_var_v3qc_v3` with called mask (pre-Arch 3 baseline):

| Metric | Value |
|---|---|
| AF MAE vs hapFIRE | **0.0150** |
| AF MAE vs FIXED recipe (AC/AN, uniform h=1/231) | **0.0118** |
| R² vs hapFIRE | 0.968 |
| R² vs recipe | 0.991 |
| Outlier rate (\|diff\|>0.1) vs hapFIRE | 0.88% |
| Outlier rate vs recipe | 0.41% |
| h-bias (cactus per-founder / PG per-founder) | 2.29× (hapFIRE = 0.88× as reference) |
| Per-founder h R² with hapFIRE | 0.19 |

For comparison, **hapFIRE itself vs FIXED recipe: MAE 0.0148, R² 0.967, outliers 1.01%** — so our cactus_em mixed-loose matches the recipe more tightly than hapFIRE does.

---

## 5. Residual error sources (in priority order)

Stratified by F_MISSING (which is now the dominant signal):

| F_MISSING bin | n records (Chr1) | %outliers vs hapFIRE | %outliers vs recipe |
|---|---:|---:|---:|
| < 0.01 | 312,325 | **0.09%** | 0.02% |
| 0.01–0.05 | 130,541 | 0.78% | 0.13% |
| 0.05–0.1 | 48,286 | 2.60% | 0.32% |
| 0.1–0.2 | 24,210 | 6.44% | 0.52% |
| 0.2–0.5 | 1,834 | **25.46%** | **8.71%** |
| ≥ 0.5 | 21 | 52% | 17% |

**87% of records are at F_MISSING ≤ 0.05 and have <0.2% outlier rate.** The remaining error is concentrated in the 6,000 records with F_MISSING > 0.1 — exactly where Beagle SNP imputation could help (task #38). Of the 4,589 outliers vs hapFIRE, ~3,300 are in F_MISSING > 0.05 records; the rest split between centromere alignment dead zone (~500) and genuine cactus-vs-1001G GT disagreement at specific records (~1,000).

---

## 6. Known QC tags on the panel

### Cactus founders (78)

Documented in `data/flag_list.tsv`:
- **2 excluded** entirely (replaced with PG-substituted versions): 101003 Set-1 (5772 mislabel of T980; k-mer Jaccard evidence), 100852 Ped-0 (9947, 10× ONT + 30% LOO disagreement)
- **5 flagged but kept** (`flag-low-cov-but-concordance-OK` / `flag-het-but-concordance-OK` / `flag-ambiguous-pair`): Fernando's source flags didn't predict assembly quality issues empirically
- **6 flagged untriaged** (`flag-untriaged-tail`): in the LOO disagreement tail (rank 1-6 of 75 worst). 7164 (Hau-0) is the WORST at 52.31% — independently identified as an extreme k-mer outlier in the v3qc_v3 composition diagnostic. Possible 5772-like mislabel.

### PG founders (153)

No formal flag list, but the v3qc_v3 k-mer composition diagnostic identified 8 extreme founders. Three of them (9977 Vezzano-2.1, 9985 Slavi-1, 10013 Lerik1-3) all have **<8× short-read coverage on old Illumina GAII** (panel median is 22.6×). These are the natural "low-data" PG founders.

One PG anomaly worth investigating: **9507 IP-Coa-0** has 23.5× HiSeq coverage (normal) but only 260 private k-mers (lowest of any founder). Possible sample mislabel or extreme similarity to another founder.

---

## 7. Outstanding work

| ID | Task | Why |
|---|---|---|
| **#10** | Re-run 4 cov10 sim regimes + SEEDMIX S2–S8 on v3qc_v3 | Validate generalization; SEEDMIX_S1 alone could be a lucky run |
| **#38** | Recompute SNP missingness on v3qc_v3 and decide Beagle | DONE: cell-level PG SNP missing = 3.00%, 18.6% of records have F_MISSING > 5%. Beagle on the biallelic-SNP subset could meaningfully reduce the high-F_MISSING outlier tail |
| **#11** | Update `RESULTS_LOG.md` and memory entries with v3qc_v3 outcomes | |
| Future | Investigate **9507 IP-Coa-0** k-mer Jaccard against panel (5772-like signature) | |
| Future | Investigate **7164 Hau-0** k-mer Jaccard (52% LOO worst, flagged) | |
| Future | Full-genome (Chr2–5) cn_full builds for v3qc_v3 | Currently only Chr1 built and EM-tested |

---

## 8. Key files and paths

### Arch 3 production cn_var — Chr1
```
arch3/chr1/
├── merged_231_chr1_final.vcf.gz                     # 231-panel VCF (biallelic, haploid, AN>0)
├── cn_var_231_arch3_chr1.cn_var.npz                 # carrier matrix (231 × 2.6M)
├── cn_var_231_arch3_chr1.cn_var_called.npz          # called mask
├── cn_var_231_arch3_chr1.meta.npz                   # chrom, pos, ref, alt, founders
├── cn_var_231_arch3_chr1_atomized.cn_var.npz        # per-base atomized (7.46M rows) — use for SNP GEA
├── cn_var_231_arch3_chr1_atomized.cn_var_called.npz
└── cn_var_231_arch3_chr1_atomized.meta.npz
```

### Production cn_full (unchanged from v3qc_v3)
```
poolfreq/data/
├── cn_full_231_v3qc_v3_mixedloose/cn_Chr1.cn.npz   # production k-mer filter, 9.78M k-mers
└── cn_full_231_v3qc_v3/cn_Chr1.cn.npz              # pre-filter (archive)
```

### v3qc_v3 panel VCFs — inputs to Arch 3 side pipelines
```
pangenie_genotyping/data/v3qc_v3/
├── pangenie_153_hetmasked_filled_bi.vcf.gz  # PG input to Arch 3 A2
├── cactus_78_bi.vcf.gz                      # cactus input to Arch 3 A3
├── founders_231_v3qc_v3.vcf.gz             # *(archive)* pre-Arch 3 merged panel
├── founders_231_v3qc_v3.haploid.vcf.gz     # *(archive)*
├── pangenie_153_hetmasked_haploid.vcf.gz   # *(archive)* — Arch 3 re-haploidizes in A2
└── pangenie_153_qc_v3.vcf.gz              # *(deprecated)* stale by-product of AC=0 pre-drop
```

> `pangenie_153_qc_v3.vcf.gz` was a stale by-product of Phase A's AC=0 pre-drop (mistakenly dropped 3.27M AC=0 records). Arch 3 A2 uses `pangenie_153_hetmasked_filled_bi.vcf.gz` directly.

---

## 9. Reproducibility

### Arch 3 (current production) — Chr1

```bash
sbatch arch3/chr1/jobA1_annotate_chr1.sh       # annotate_vcf + biallelic catalog (feeds A2+A3)
sbatch arch3/chr1/jobA2_pg_chr1.sh             # PG side: transfer_id → convert-to-biallelic → haploidize
sbatch arch3/chr1/jobA3_cactus_chr1.sh         # cactus side: transfer_id → convert-to-biallelic
sbatch arch3/chr1/jobA4_merge_chr1.sh          # merge 231 + fill-tags + AN=0 filter
sbatch arch3/chr1/jobA5_build_cnvar.sh         # build cn_var + cn_var_called + meta
sbatch arch3/chr1/jobA6_project_h_seedmix.sh   # project h through new cn_var → AF TSV
sbatch arch3/chr1/jobA7_compare_vs_hapfire.sh  # validate vs hapFIRE (4-tuple join)
# Optional: atomized cn_var for SNP GEA
sbatch arch3/chr1/jobD1_atomize_cnvar.sh
```

### *(archive)* v3qc_v3 build scripts (cn_full only — still used in Arch 3)

```bash
sbatch poolfreq/scripts/build_cn_full_v3qc_v3_chr1.sh
sbatch poolfreq/scripts/build_cn_full_mixedloose_v3qc_v3_chr1.sh
```

The v3qc_v3 panel VCF build scripts (`build_v3qc_v3_phase_a.sh`, `phase_b.sh`, `haploidize_*.sh`, `build_cn_var_v3qc_v3.sh`) are superseded by Arch 3's A1–A5. They remain on disk as archive.

---

## 10. Bugs caught in this work (in case any recur elsewhere)

1. **Multi-allelic OR-semantics in `bcftools view -e 'INFO/X...'`**: filters on per-ALT INFO fields drop the entire record if any ALT triggers the expression. Fix: `bcftools norm -m -any | sort` before filtering. (`SESSION_2026-05-17.md`)
2. **`./.` treated as REF in cn_var**: AC/AN denominator was always 231, leading to AF under-estimation at high-F_MISSING records. Fix: cn_var_called sparse matrix + AC/AN projection.
3. **`np.array(ref_strings)` overflow on long SV alleles**: numpy infers fixed-width string dtype, tries to allocate `2.3 TiB` for a 100kb REF. Fix: `dtype=object`.
4. **`bcftools merge` mixed-ploidy + AC=0 + scale bug**: full-dataset merge of haploid cactus + diploid PG drops PG GTs at AC=0 records. Not reproducible on small subsets. Fix: pre-haploidize PG.
5. **bcftools `+setGT -- -i 'GT="het"' -n .` with an apostrophe in a comment**: the bash single-quoted string closes prematurely. Fix: avoid apostrophes in awk/bcftools-inline scripts.
6. **`alpha` post-EM correction creates per-founder noise**: increasing the K_f exponent reduces side-averaged h-bias but blows up per-founder std on the PG side. Fix: don't optimize for h-bias as a proxy — use mixed-loose at the data layer instead.

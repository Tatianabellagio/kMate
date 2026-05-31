# 231-panel missingness by variant class

**Generated**: 2026-05-21
**Source data**: `data/var_pa_231_v3qc_v3.var_called.npz` + `var_pa_231_v3qc_v3.meta.npz` — the
v3qc_v3 `norm`-merge matrix, now **superseded/archived** (`docs/PIPELINE_STATE.md` §0). The
production matrix is `panel/arch3/chr{N}/var_pa_231_arch3_chr{N}`. arch3 reproduces the same
per-cell genotype QC (het→`.`, GQ≥20), so the class-level missingness characterization below
carries over; the absolute record counts will differ under the arch3 biallelic catalog.
**Plot**: `benchmarks/p80/results/missingness_by_var_class_231panel.png`
**SLURM log**: `benchmarks/p80/logs/miss_231_63224.out`

## TL;DR

PanGenie SV genotyping is the limiting factor on the 231-panel. PG side has **25% missing on SVs** vs 3% on SNPs. Cactus side is ~95% called for every variant class. The "treat-`.`-as-REF" projection convention used by `cactus_em` window/star2 mode would silently under-call SV AF at production scale — this is not a defensible default for SVs.

## Numbers

### Overall (231 founders × 6,293,085 records, density 93.86%)

| class | n records | % panel | mean #miss / 231 | % fully called | % > 50 missing |
|---|---:|---:|---:|---:|---:|
| SNP   | 4,128,028 | 65.60% |  9.36 | **34.67%** |  4.43% |
| indel | 1,547,744 | 24.59% | 16.06 |     19.72% |  6.18% |
| **SV** |  **617,313** |  **9.81%** | **41.82** |  **1.54%** | **23.11%** |
| ALL   | 6,293,085 |  —     | 14.19 |     27.75% |  6.69% |

### Per-side decomposition

| class | cactus #miss (of 78) | % cactus called | PG #miss (of 153) | % PG called |
|---|---:|---:|---:|---:|
| SNP   | 4.12 | 94.72% |  5.25 | 96.57% |
| indel | 2.60 | 96.67% | 13.46 | 91.20% |
| **SV** | **3.21** | **95.89%** | **38.61** | **74.76%** |

Cactus side: clean and class-uniform (~95% called everywhere).
PG side: SVs are catastrophically worse than SNPs (25% missing vs 3%).

### Records retained by max-#missing threshold

| max_miss | SNP %kept | indel %kept | SV %kept | ALL %kept |
|---:|---:|---:|---:|---:|
|   0 |  34.67% |  19.72% |   **1.54%** |  27.75% |
|   1 |  47.99% |  28.19% |   2.35% |  38.64% |
|   5 |  65.61% |  42.97% |   4.38% |  54.03% |
|  10 |  74.72% |  53.29% |   6.55% |  62.76% |
|  25 |  88.53% |  75.30% |  18.77% |  78.43% |
|  50 |  95.57% |  93.82% |  76.89% |  93.31% |
| 100 |  99.61% |  99.54% |  96.69% |  99.31% |
| 150 |  99.92% |  99.87% |  97.80% |  99.70% |

A strict `F_MISSING ≤ 0.05` filter (~max 11 missing out of 231) **would discard ~93% of SV records.** Non-starter for SVs.

## Implications

1. **PanGenie + short reads is the SV genotyping bottleneck.** Not k-mer rebalancing, not the cactus side, not the panel mixture math. The PG side has structural failure on SVs and downstream methods inherit it.

2. **The MAR (called-mask) projection is the only defensible default for SVs.** At a typical SV record on the 231-panel, ~75% of PG founders are `.`. "Treat-`.`-as-REF" silently treats ~115 missing PG founders as non-carriers → systematic SV AF under-call.

3. **Star2 (window/block) projection on production must be patched to use the called mask.** The current code path in `block_em.py:project_blocks_to_records` skips the called-mask normalization, which is fine for cactus-only p80 but wrong for production SVs. Either both modes use the mask or neither does — letting `--block-mode` silently toggle the missingness model is a footgun (`HIGH` severity for SV-touching analyses).

4. **p80 doesn't probe this failure mode.** p80 is cactus-only, where missingness is ~4% across all classes. The interesting pathology lives entirely on the PG-SV intersection that p80 excludes by design. Numbers from p80 about SV inference quality (especially the "SVs better than SNPs" result) **do not transfer** to production.

5. **~23% of production SV records sit at F_MISSING > 50/231.** At very high missingness, AF estimator variance grows independently of which convention you use. Either propagate per-record uncertainty downstream, or acknowledge that the high-F_MISSING SV slice is unreliable. See `PIPELINE_STATE.md` for production-state context.

## Methodology notes

- `var_called[f, r] = 1` iff founder f's GT at record r is not `.`. Built by `build_var_pa.py` alongside `var_pa`.
- Cactus side identified as Accession_IDs in `data/sv_panel_to_accession_id.tsv` minus `data/exclude_list.txt` (78 of 80 cactus assemblies; 5772 and 9947 PG-substituted per v3qc decisions).
- Variant class: SNP = `ref_len==1 AND alt_len==1`; SV = `max(ref_len, alt_len) ≥ 50`; indel otherwise.
- This is panel-side missingness (founder GT availability). Pool-side sampling noise is a separate, orthogonal source of AF variance.

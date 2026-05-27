# Panel overlap — v3 vs GrENE-Net 231 SNP catalog (results)

**Date**: 2026-05-14
**Inputs**:
- v3 SNP positions (biallelic): `data/v3_snp_positions.nochr.tsv` (4,426,208 positions on Chr1-5)
- GrENE-Net 231 SNP positions: `data/grenenet231_snp_positions.tsv` (3,235,480 positions on Chr1-5)

Both panels cover the **same 231 founders** (80 cactus + 151 PG = exactly the GrENE-Net 231). Direct (chrom, pos) intersection — no founder-subset complication.

## Headline counts

| comparison | n |
|---|---:|
| v3 biallelic SNP positions (Chr1-5) | 4,426,208 |
| GrENE-Net 231 SNP positions (Chr1-5) | 3,235,480 |
| shared (in both) | 1,781,226 |
| **% of GrENE-Net SNPs in v3** | **55.05%** |
| % of v3 SNPs in GrENE-Net | 40.24% |
| only in v3 (v3 calls, GrENE-Net misses) | 2,644,982 |
| only in GrENE-Net (GrENE-Net calls, v3 misses) | 1,454,254 |

Compare with the OLD test (syri 82-accession panel vs GrENE-Net): syri covered 66.7% of GrENE-Net SNPs. v3 covers 55% — LESS, despite having the full 231 founders. Different pipelines produce different SNP catalogs (see § "Why" below).

## Three-way decomposition of the GrENE-Net catalog

Every GrENE-Net SNP falls into exactly one of three categories. Stratified by whether v3 has the SNP and (if not) whether cactus has a bubble at that position:

| chrom | total GN SNPs | cat1: shared | cat2: in-bubble dropped | cat3: no bubble |
|---|---:|---:|---:|---:|
| Chr1 | 827,765 | 468,661 (56.6%) | 308,507 (37.3%) | 50,597 (6.1%) |
| Chr2 | 536,387 | 277,051 (51.6%) | 224,419 (41.8%) | 34,917 (6.5%) |
| Chr3 | 617,020 | 334,290 (54.2%) | 230,640 (37.4%) | 52,090 (8.4%) |
| Chr4 | 531,361 | 286,182 (53.9%) | 198,029 (37.3%) | 47,150 (8.9%) |
| Chr5 | 722,947 | 415,042 (57.4%) | 263,571 (36.5%) | 44,334 (6.1%) |
| **TOTAL** | **3,235,480** | **1,781,226 (55.05%)** | **1,225,166 (37.87%)** | **229,088 (7.08%)** |

- **cat1 (55.05%)** — SAME SNPs: in both v3 cn_var and GrENE-Net. Currently usable; both pipelines call these. Pred-target for hapFIRE-projected-through-cn_var_v3.
- **cat2 (37.87%)** — SNPs INSIDE a cactus pangenome bubble but absent from `cn_var_v3`. Cactus's graph knows about that position (creates a bubble there) but the SNP record was dropped downstream. **Probable cause: `vcfbub -l 0 -r 100000` filter** strips nested SNPs/indels, keeping only top-level snarls. Recoverable by rebuilding cn_var_v3 with `vcfbub -l max_level` (or skipping vcfbub entirely).
- **cat3 (7.08%)** — SNPs OUTSIDE all cactus bubbles. Cactus has NO bubble at that position. Two flavours per the Mb-hotspot analysis on Chr1:
  - ~25% concentrated in known cactus alignment dead zones (centromere Mb 14-17, knob Mb 21-23.5; Chr4 NOR at Mb 0-5; pericentromeric regions of all chroms)
  - ~75% scattered across "easy" arm regions where all 80 cactus founders happened to agree with TAIR10 — so cactus's deconstruct didn't create a bubble — but the 151 PG founders DO have variation there per GrENE-Net's pipeline.

## Spatial pattern of cat3 (no-bubble) SNPs

Top Mb hotspots for the 229K outside-bubble SNPs:

| chrom | dominant Mb hotspots | n in cat3 |
|---|---|---:|
| Chr1 | Mb 13–16 (peri/centromeric) | 50,597 |
| Chr2 | Mb 5–8 (peri/centromeric) | 34,917 |
| Chr3 | Mb 14–16 + Mb 20 (peri/centromeric + knob3) | 52,090 |
| Chr4 | **Mb 2 alone has 14,290** (NOR / heterochromatic knob region) | 47,150 |
| Chr5 | Mb 10–14 (peri/centromeric) | 44,334 |

For Chr1: 19.6% of cat3 in centromere (Mb 14–17), 5.8% in knob (Mb 21–23.5), **74.6% scattered across the arms**.

## Why v3 covers LESS of GrENE-Net than syri did

| panel | pipeline | % of GrENE-Net SNPs covered |
|---|---|---:|
| syri 82-accession (old test) | per-assembly long-read SNP calls vs TAIR10 (dense, per-position) | 66.7% |
| **v3 cactus + PG (this study)** | cactus pangenome bubbles + PanGenie short-read genotyping (graph-bubble-restricted) | **55.05%** |

syri reports every position where any of the 82 long-read assemblies differs from TAIR10. v3's cn_var_v3 only retains SNPs that survive: cactus deconstruct → vcfbub -l 0 → bcftools norm -m -any. Each of those steps drops some calls. Net effect: v3 has ~4.4M SNP positions vs syri's ~5M, but they're a DIFFERENT subset.

## Implications for the hapFIRE-projected-through-cn_var_v3 pipeline

When `hapfire_v2panel_v3proj` / `hapfire_v3panel` / `hapfire_coarse_v3proj` projects hapFIRE's per-block h through cn_var_v3:

- **55.05% of GrENE-Net SNPs** are kept (cat1). These are the positions where both methods estimate.
- **37.87%** (cat2) are silently dropped because the SNP-level record doesn't exist in cn_var_v3 — even though the pangenome graph has a bubble there.
- **7.08%** (cat3) are also dropped because the pangenome has no bubble there.

The n_pred discrepancy we kept seeing in `FINAL_RESULTS_cov10_v3.ipynb` (hapfire n_pred ≈ 635-666K vs cactus_em n_pred ≈ 671K) is the n=441,280 SNP intersection × 1 chrom — the cn_var_v3 record set, not the GrENE-Net set.

## Production recovery path

The biggest single fix would be **rebuilding cn_var_v3 with `vcfbub -l max_level`** (or skipping vcfbub entirely):
- v3 SNP coverage of GrENE-Net would jump from **55% → 93%** (cat1 + cat2)
- Trade-off: cn_var_v3 records grow from 7.69M → ~15M+ (atomization keeps every nested SNP/indel), with downstream memory/compute cost
- PanGenie genotyping at nested-bubble positions may have correctness issues at multi-allelic complex sites
- Not tested yet on Chr1

The remaining 7.08% (cat3) needs either:
- A different cactus pipeline that bubbles every variant position regardless of cactus-founder uniformity (changes whole pangenome topology)
- A separate "PG-only variant" track that augments cn_var_v3 with positions PanGenie called outside cactus bubbles

## Reproduce

```bash
# 1. Extract v3 position lists (~5 min SLURM)
sbatch extract_v3_positions.sh

# 2. Open analysis notebook (3 Venns + per-chrom + cat1/cat2/cat3 split)
overlap_v3_vs_grenenet.ipynb
```

## Files

| File | Purpose |
|---|---|
| `extract_v3_positions.sh` | SLURM job to extract v3 position lists |
| `data/v3_*_positions.tsv` (+ `.nochr.tsv`) | v3 (chrom, pos) lists (3 versions × 2 chrom-naming) |
| `data/grenenet231_snp_positions.tsv` | 3.24M GrENE-Net 231 SNP positions, from old test |
| `overlap_v3_vs_grenenet.ipynb` | Analysis notebook (Venns, per-chrom, cat1/2/3 split) |
| `README.md` | Methodology + chrom-prefix gotcha + file map |
| `results_summary.md` | This file |

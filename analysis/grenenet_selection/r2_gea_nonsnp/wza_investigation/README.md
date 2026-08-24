# WZA investigation — was WZA the right choice for phase-1 block aggregation?

Owner: T. Bellagio · started 2026-06-16 · scope: **gen1, bio1** (kept narrow on purpose).

Re-examines the GrENE-Net phase-1 use of the **WZA** (Booker et al. 2024, *Mol.
Ecol. Resour.*) to pool per-SNP GEA p-values into haplotype-block windows, after
the NaN-p-value problem flagged in the Booker email thread.

## TL;DR
WZA is not the problem; **our window definition is**. WZA's SNP-number correction
assumes a reasonably bounded, densely-sampled SNP-per-window distribution. Our
hapFIRE/BigLD haplotype blocks span **1 → ~3,000 SNPs/block** (1.05M SNPs in
16,674 blocks; centromeric LD blocks are huge). The deg-2 correction polynomial
turns down in that sparse heavy tail → predicts negative SD → NaN. The fix used by
both the WZA author (email) and the RepAdapt application is to **cap/downsample
SNPs per window** — a switch built into WZA (`--sample_snps`) that phase-1 never
set. Phase-1 instead dodged the NaN with a **deg-7** polynomial, which over-fits
the tail and inflates significance (CAM5 deg-7 p=3e-8 vs deg-2 p≈1e-3).

---

## Worry 1 — what LD threshold defines the 16,674 blocks? **ANSWERED**

The blocks are **hapFIRE** haplotype blocks (`scratch/HapFIRE`), built with a
**two-level LD partition** (`haplotype_generation.py`, `BigLD.R`):

| level | method | threshold | output |
|---|---|---|---|
| **independent blocks** | `CompleteLDPartition` (r² to up/down-stream neighbours) | **r² = 0.1** (`-r`), window **100** (`-w`) | `*.independent_genomewide_partition.txt` |
| **fine blocks** (the 16,674 we use) | **BigLD** (`gpart`, `CLQmode="density"`, `hrstType="fast"`) | **`CLQcut` r² = 0.5** default (`-c`), `MAFcut` (`-m`) | `*.fine_genomewide_partition.txt` |

So: independent breakpoints at **r² 0.1**, then BigLD sub-partitions each at
**r² 0.5** (density mode). Our 16,674-block map (`gea_grene-net/ARCHIVE/
linages_wza_picmin/kendall_0_w_id_n_blocks.csv`, used by `lib.assign_ld_blocks`)
matches the **fine** partition: 1,054,574 SNPs, median **9** SNPs/block, span
median **666 bp**, but a heavy tail — p99 1,066 SNPs / 117 kb, max **3,028 SNPs /
586 kb** (centromeric LD). 255 singleton blocks; 14 blocks > 2,000 SNPs.

> Caveat: the exact CLQcut / MAFcut / window Tati passed lived on the old Carnegie
> cluster (`/home/tbellagio/HapFM/output_fine_genomewide_partition.txt`, gone). The
> block-size profile is consistent with the hapFIRE **fine BigLD r²=0.5** default.
> The frozen SNP→block map is what we use, so the partition itself is reproducible
> from that file regardless.

**This is the root cause.** Centromeric/pericentromeric LD makes BigLD emit a few
enormous blocks, producing the 1→3,000-SNP spread that breaks WZA's correction.

## Worry 2 — does RepAdapt use genes as windows? **YES**

RepAdapt (Whiting et al.; `github.com/JimWhiting91/RepAdapt`, branch `main`,
`R/01_GEA_pipeline/collate_GEA_to_gene_WZA_empirical_pvals.R`):
- **Windows = genes + 500 bp flank** (`gene_flank_size`); SNPs assigned by overlap
  (`foverlaps(snp_bed, gene_bed)`).
- **They cap SNPs per gene**: `snp_per_gene <- quantile(table(gene_id), probs=...)`
  (floor 10), then average the WZA over **100 random down-samples** — i.e. they
  USE WZA's SNP-cap mechanism by default.

Booker's own `general_WZA_script.py` has the same cap (`--sample_snps -1` → 75th
pct, `--resamples 100`). His email to Tati recommended the same: cap ~2000 + deg-7.
**Phase-1 set neither cap** (`--sample_snps` unset → max 1e6), only deg-7.

## Worry 3 — does WZA over-weight low-SNP windows (e.g. CAM5 ~16 SNPs)?

Tested empirically — see `analysis/grenenet_selection/r2_gea_nonsnp/wza_investigation/results/`
(`overweighting_by_snpbin.csv`, `fig2`, `fig3`) and the summary below.

---

## Code
- `wza_core.py` — Booker WZA, parameterized (degree, roller/minEntries, SNP cap).
  Raw weighted-Z byte-identical to canonical; only the correction is exposed.
- `run_investigation.py` — gen1 SNP×bio1: sweeps degree {2,7} × cap {none, 2000,
  q95, q75}; reports NaN rate, Bonferroni hits, CAM5 (block 2_1265 / 2_1264); plus
  the low-SNP-overweighting test. Outputs → `analysis/grenenet_selection/r2_gea_nonsnp/wza_investigation/results/`.

## Results
(filled in after the run — see `RESULTS.md`.)

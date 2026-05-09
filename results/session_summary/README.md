# Session results — methods compared

| Method | What it does | Evidence | Resolution |
|--------|--------------|----------|-------------|
| **freqk** | Per-record k-mer-presence ratio. No panel-EM, no block model. Indexes the same Beagle-imputed VCF (merged_v2) the recomb truth uses, computes per-bubble allele frequency directly from k-mer counts. Pos offset of +2 applied at join time (freqk_pos = VCF_pos − 2). | k-mers (alignment-free, panel-free) | Per-bubble (≥31 bp due to k-mer flanks) |
| **hapfire_coarse** | Existing hapFIRE — partition Chr into ~30 coarse independent LD blocks, run HARP per coarse block + analytical projection to fine-block AFs. | Read-level SNPs (alignment-based) | Coarse LD block |
| **hapfire_fine** | NEW. Same HARP but per BigLD fine block (~7 kb). Group-size-corrected projection. | Read-level SNPs (alignment-based) | BigLD fine block |
| **bigld_panel** | cactus_em on the 231-founder simplex per BigLD block, using PanGenie-style bubble-allele k-mers (strict allele-uniqueness + 16/32 cap). | Bubble-allele k-mers (panel-derived) | BigLD fine block |
| **clean_cc5** | cactus_em block-haplotype EM. Reconstruct each panel haplotype's full block sequence from VCF, k-merize, dedup at class level, filter to k-mers carried by ≥5 founders. | Block-haplotype k-mers (panel-derived) | BigLD fine block (class-level h then scattered to founders) |
| **clean_smooth** | clean_cc5 + Li-Stephens HMM smoothing across blocks (passes=10, α=0.2). | same | same |
| **dose_aware** | clean reconstruction with k-mer multiplicity (`cn[c, k]` = count, not 0/1). | same | same |
| **global** | cactus_em --block-mode recomb_global. Single chrom-wide h. | Block-haplotype k-mers | Whole chrom |
| **window** | cactus_em --block-mode recomb_window. Fixed 200-kb windows. | Block-haplotype k-mers | 200-kb window |

## Datasets
- `af_long.parquet`: long-form per-record (chrom, pos, ref_len, alt_len, regime, method, est_af, truth_af, var_type, sv_size, max_r2_xwu, ld_tier).
- `af_summary.parquet`: pre-computed R²/MAE/RMSE/slope per (regime, method, var_type, ld_tier).
- `identity_metrics.parquet`: per-block h-vector identity (L1, JSD, cosine, top-5 recall) at founder + class level.

## Regimes (3, all chr1, cov=50)
- `n200_g1`: 200 founders, 1 generation of recomb (high diversity, weak block structure)
- `n50_g1`: 50 founders, 1 gen of recomb (moderate)
- `n50_g3`: 50 founders, 3 gens of recomb (low diversity, strong fine-block structure — hardest regime)

LD tiers (max r² between SV and any panel SNP within 50 kb): `weak` < 0.2, `medium` 0.2-0.5, `strong` 0.5-0.8, `perfect` ≥ 0.8.

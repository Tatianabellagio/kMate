# `benchmarks/p231/` — 231-founder benchmark (paper headline)

The main kMate accuracy benchmark on the **full 231-founder heterogeneous panel**
(78 cactus long-read + 153 PanGenie short-read founders). Companion to
`benchmarks/p80/` (the all-long-read homogeneous control that isolates the panel-
imbalance effect). Method under test: **filt2 + `--kmer-weight uniform` + `--normalize per_founder`, `--unit chrom`** (the corrected front-runner).
`inv_mb` (ω_k=1/m_b) is retained only as a clearly-LABELED legacy A/B column; it is
SUPERSEDED for global by per-founder normalization + uniform weighting.

## Provenance / audit (verified 2026-05-27)

Everything traces to the current **v3qc + arch3** lineage; nothing stale is forced in.

- **Canonical VCF** `panel/arch3/chr1/merged_231_chr1_final.vcf.gz`: v3qc QC (cactus side
  `view -s ^5772,9947` drops the 2 flagged assemblies; PG side `setGT GQ>=20`),
  arch biallelic decomposition + merge + AN=0 filter, **unimputed**, het→missing
  masked. 1.90M SNP / 637K indel / 84K SV. This is the latest QC'd arch panel.
  > **NOTE (2026-06-02):** this benchmark predates the segregating-only filter. The
  > production panel was later filtered to drop monomorphic records (`AC=0 || AC=AN`),
  > so the current Chr1 `merged_231_chr1_final.vcf.gz` has 2,154,423 records
  > (1.56M SNP / 534K indel / 59K SV), not the pre-filter counts above. The benchmark
  > numbers below were computed on the pre-filter panel; segregating-site accuracy is
  > unchanged by the filter (per-record projection). See `docs/PIPELINE_STATE.md` §0.
- **FASTAs** `fastas_231/`: rebuilt by `bcftools consensus -H 1` from that VCF
  (haploid; founder order identical to var_pa). NOT the stale v3 fasta dirs.
- **kmer_pa** `data/kmer_pa_p231[_filt2]/`: rebuilt from the SAME VCF via
  `build_kmer_pa.py` + the pang_135 k-mer dictionary (135-asm graph, matches
  merged_231's annotation topology), `--treat-missing-as-n`. Validated against the
  then-current `kmer_pa_231_v3qc_v3_filt2` by `03c_compare_kmer_pa.py` (kmer_pa is
  consensus-derived ⇒ representation-invariant ⇒ expected ~identical). *(Note: that
  v3qc_v3 matrix is now archived; current production is `kmer_pa_231_arch3_filt2inv`,
  in-house index + filt2inv — `docs/PIPELINE_STATE.md` §0. This benchmark intentionally
  uses plain `filt2` and the pang_135 PG dictionary, not the production filt2inv/in-house index.)*
- **var_pa** (projection + truth target): REUSED `panel/arch3/chr1/var_pa_231_arch3_chr1`
  in TWO arms — `_atomized` (7.46M per-base records, SNP-level GEA benchmark) and
  raw (2.62M records, SNP/indel/SV classes). Both built from the canonical VCF.
- **het handling**: both kmer_pa and var_pa mask het→missing (Arouisse-2020
  precedent); consistent across the k-mer and variant sides.

## Key design choices vs benchmarks/p80

1. **RANDOM crossovers** at the A. thaliana rate (4 cM/Mb, `DEFAULT_RECOMB_RATE`),
   NOT forced at hapFIRE BigLD block boundaries — forcing crossovers at the LD
   blocks a block-based method keys on would be circular for a benchmark. The
   LD-block crossover option was removed entirely (2026-05-30); `sample_crossovers()`
   always places Poisson(L·rate) crossovers at uniform-random positions.
2. **Truth on both var_pas** so each projection arm joins its truth 100% on
   (chrom,pos,ref_len,alt_len) — the est↔truth consistency gate. (Pre-arch3 sim
   truth files joined only 16–21% to atomized and must not be reused; the old
   `visor_freqk` sim data they lived in was removed 2026-05-30.)
3. **kmer_pa rebuilt from arch3** for single-VCF provenance (vs production reuse).

## Regimes (match benchmarks/p80, cov10, seed 42, Chr1)

`n50_g0`, `n231_g0`, `n50_g1`, `n231_g1`, `n50_g3`, `n50_g3_dom500` (dom = one
individual at 50% of pool reads). g0 = perfect founder mix (no recombination).

## Pipeline

```
05_build_fastas_p231.sh         231 consensus FASTAs from merged_231 (array 1-231)
03_build_kmer_pa_p231.sh        kmer_pa from merged_231 + pang_135 dict
03b_build_kmer_pa_filt2_p231.sh filt2 (drop ac<2 singletons) -> front-runner kmer_pa
03c_compare_kmer_pa.py          validate rebuilt vs production kmer_pa
06_run_sim_p231.sh N G          mosaics (RANDOM crossovers) + VISOR reads + truth x2
06b_run_sim_p231_skewed.sh      dom500 variant
07c_run_kmate_filt2_mb_p231.sh REGIME CNVAR [WEIGHT]   EM (--unit chrom, uniform + --normalize per_founder; inv_mb = legacy A/B only) + project
score_p231.py                   MAE/RMSE/R2/outlier by regime x cnvar x weight x class
submit_all_p231.sh              the full DAG
```

## Verification gates

- G1/G2 (FASTAs): 231 built, every one applied >0 variants, founder IDs == kmer_pa. ✅
- G4 (truth): `recomb_truth_atomized` n rows == 7,464,709; `recomb_truth_raw` == 2,617,370.
- G6 (join): est↔truth inner join == var_pa record count (100%).
- kmer_pa validation: rebuilt vs production overlap ~100%, carrier agreement ~100%.

## Note

Reads/kmer_pa use a haploid panel; the founder-list source for the mosaic builder
is the production kmer_pa meta (founders identical to merged_231, order-verified) —
immaterial to provenance (only names are used).

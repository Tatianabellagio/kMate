# panel/

Construction of the **231-founder reference panel** that kMate estimates against
(78 cactus + 153 PanGenie-genotyped founders → `merged_231_chr1_final.vcf.gz`).

The panel is built from two sides and then decomposed into the per-record copy-number
matrices (`cn_var`) the estimator projects through:

| dir | role |
|---|---|
| `pangenie_index/` | In-house k-mer index tables for the panel (was `kmer_index/`). Per-founder k-mer fingerprints used by the k-mer count step. Heavy `*.tsv.gz` tables are gitignored (regeneratable from the PanGenie index). |
| `pangenie_genotyping/` | Short-read side: PanGenie genotyping of the 153 PG founders (download → preprocess → genotype → QC → merge). Produces the `v3qc_v3` haploidized panel VCFs under `data/v3qc*/`. Leave-one-out (LOO) concordance tooling lives here too. |
| `arch3/` | **arch decomposition** (annotate_vcf + convert-to-biallelic, NOT `bcftools norm -m -any`): turns the merged multi-allelic panel VCF into per-record / per-base (atomized) `cn_var` matrices. Chr1 is production-validated; Chr2–5 is the open build (`jobA1..A5`, `jobD1` atomize). See `docs/INVESTIGATION_CN_VAR_DECOMPOSITION.md` for why arch replaced norm. |
| `imputation/` | **Deprecated** (Beagle-based founder imputation; rejected — see `docs/RESULTS_LOG.md` 2026-05-06 literature lock-in). Tooling kept for archaeology; content untracked. |

## Pipeline (Chr1, production)

1. **Long-read / cactus side** → 78 cactus founders (pang_135 → pang_82 → QC-filtered).
2. **Short-read / PanGenie side** (`pangenie_genotyping/`) → 153 PG founders, haploidized → `founders_231_v3qc_v3.haploid.vcf.gz`.
3. **Merge** → `arch3/chr1/merged_231_chr1_final.vcf.gz` (231 founders).
4. **arch decomposition** (`arch3/`) → `cn_var_231_arch3_chr1.{cn_var,cn_var_called,meta}.npz`
   and the atomized (per-base) `*_atomized.*` variants.

The resulting `cn_var` (this folder) plus `cn_full` (k-mer copy number, in `data/`) are the
two inputs to `src/per_sample_per_chrom.py`. See root `HANDOFF.md` for the production recipe.

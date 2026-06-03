# panel/

Construction of the **231-founder reference panel** that kMate estimates against
(78 cactus + 153 PanGenie-genotyped founders → `merged_231_chr1_final.vcf.gz`).

The panel is built from two sides and then decomposed into the per-record presence/absence
matrices (`var_pa`) the estimator projects through:

| dir | role |
|---|---|
| `pangenie_index/` | In-house k-mer index tables for the panel (was `kmer_index/`). The production `kmer_pa` index (`ours_Chr{N}`); the PanGenie genotyping index is a comparator only. Heavy `*.tsv.gz` tables are gitignored (regeneratable via `scripts/build_kmers_tsv.py`). |
| `pangenie_genotyping/` | Short-read side: PanGenie genotyping of the 153 PG founders (download → preprocess → genotype → QC). Produces the per-side VCFs under `data/v3qc_tmp/` that `arch3/` ingests (its final `founders_231_v3qc_v3.*` naive-norm merge is superseded — see that dir's README). Leave-one-out (LOO) concordance tooling lives here too. |
| `arch3/` | **arch decomposition** (annotate_vcf + convert-to-biallelic, NOT `bcftools norm -m -any`): re-merges the haploid biallelic side VCFs by symbolic-ID into `merged_231_chr{N}_final.vcf.gz` and builds the per-record / per-base (atomized) `var_pa` matrices. All 5 chroms built on disk (Chr1 validated; Chr2–5 present as of 2026-05-30/31, validation pending) — `jobA1..A5`, `jobD1` atomize. See `docs/INVESTIGATION_CN_VAR_DECOMPOSITION.md` for why arch replaced norm. |
| `imputation/` | **Deprecated** (Beagle-based founder imputation; rejected — see `docs/RESULTS_LOG.md` 2026-05-06 literature lock-in). Tooling kept for archaeology; content untracked. |

## Pipeline (production)

1. **Long-read / cactus side** → 78 cactus founders (pang_135 → pang_82 → QC-filtered), haploid biallelic side VCF.
2. **Short-read / PanGenie side** (`pangenie_genotyping/`) → 153 PG founders, het-masked + haploidized side VCF (`data/v3qc_tmp/`).
3. **arch3 symbolic-ID merge** (`arch3/`, A2–A4, `bcftools merge --merge none`) → `arch3/chr{N}/merged_231_chr{N}_final.vcf.gz` (231 founders; the canonical production panel VCF). **Segregating-only** since 2026-06-02 — A4 drops monomorphic records (`AC=0 || AC=AN`); 8,489,646 records genome-wide. See `docs/PIPELINE_STATE.md` §0.
4. **arch decomposition** (`arch3/`, A5/D1) → `var_pa_231_arch3_chr{N}.{var_pa,var_called,meta}.npz`
   and the atomized (per-base) `*_atomized.*` variants (the on-disk atomized is Chr1-only and **stale** vs the segregating filter — rebuild before use).

The resulting `var_pa` (this folder) plus `kmer_pa` (k-mer presence/absence, in `data/`) are the
two inputs to `src/per_sample_per_chrom.py`. See `docs/PIPELINE_STATE.md` §0/§0.1 for the production recipe.

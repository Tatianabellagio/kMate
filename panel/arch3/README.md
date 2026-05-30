# arch3 — multi-allelic → biallelic decomposition (the "arch" step)

As-built documentation of how the merged 231-founder **multi-allelic** pangenome VCF
is decomposed into the **biallelic** per-record matrices (`var_pa`, `var_called`) that
the kMate estimator projects through, plus the per-base **atomized** variant.

> Design rationale and the bug this avoids: `docs/INVESTIGATION_CN_VAR_DECOMPOSITION.md`.
> This README is the *implemented* pipeline; that doc is the *why*.

## Why not `bcftools norm -m -any`

Naïve decomposition (`bcftools norm -m -any`, and also `--atomize` / `vcfwave`) splits
each multi-allelic record by **pairwise REF↔ALT alignment**. On a pangenome graph this
**scatters carriers across shifted positions** and **silently drops carriers where ALT
paths converge** to the same atomic variant. Measured effect: up to ~99% carrier loss at
SNPs co-located with multi-allelic indels (e.g. Chr1:10421645 read AC=2/231 vs truth 229).
Because our simulation truth was computed *from the same var_pa*, the bug was invisible to
closed-loop validation for months. See the investigation doc §2–§3.

## The approach: symbolic-ID propagation ("arch", borrowed from HPRC)

Instead of re-aligning, we use the **HPRC human-pangenome** decomposition strategy
(`prepare-vcf-MC` + eblerjana `convert-to-biallelic.py`): every atomic (biallelic) variant
nested inside a graph bubble carries a **symbolic ID**, and per-sample genotypes are
propagated **multi-allelic → biallelic by matching those IDs** — no alignment, so
convergence and position-shift are handled by construction.

External tools (expected under `external/`, which is **gitignored** — not shipped; obtain
them per "Prerequisites" below):
- `genotyping-pipelines/prepare-vcf-MC/.../annotate_vcf.py` — HPRC; emits the annotated
  multi-allelic VCF (each ALT tagged with `INFO/ID` = the symbolic IDs of the atomic
  variants on that ALT's graph path) **and** the biallelic catalog (one record per atomic
  variant, keyed by the same IDs).
- `pangenie-tools/.../convert-to-biallelic.py` — eblerjana/HPRC; for each sample's GT on a
  multi-allelic record, sets that sample's GT on every biallelic catalog record whose
  symbolic ID appears in the called ALT path.
- `panel/arch3/transfer_id_annotation.py` — **ours**; copies `INFO/ID` from the
  annotated catalog onto the PanGenie/cactus genotyped VCFs (which were genotyped on the
  *same* graph but lack the IDs), so `convert-to-biallelic` has IDs to match on. Custom
  because `bcftools annotate -c INFO/ID` corrupts the angle-bracketed graph-node IDs.

**Inherited limitation:** the atomic catalog from `annotate_vcf` is itself built with
`vcfwave` internally, so ~0.8% of atomic variants are missing (this is the Arch-3 fallback
caveat from the investigation doc, accepted for production). The Arch-1 graph-native path
(`vg deconstruct -a -e` + `resolve-nested-genotypes`) would close that gap and is the
documented upgrade path.

## Pipeline DAG (Chr1 production; scripts in `chr1/`)

```
                          A1  annotate_vcf  (135-sample pangenome VCF + GFA)
                          │   → annotated multi-allelic VCF + biallelic catalog
                          │     (both carry the symbolic INFO/ID)
              ┌───────────┴───────────┐
              ▼                       ▼
   A2  PG side (153)         A3  cactus side (78)
   transfer_id →             transfer_id →
   convert-to-biallelic →    convert-to-biallelic →
   fill-tags →               (already haploid)
   haploidize (het→.)        sort/bgzip
              └───────────┬───────────┘
                          ▼
                          A4  bcftools merge --merge none  (→ 231, stays biallelic)
                          │   fill-tags ; drop AN=0
                          ▼
                 merged_231_chr1_final.vcf.gz   ← canonical biallelic panel
                          │
              ┌───────────┴────────────┐
              ▼                         ▼
   A5  build_var_pa.py        D1  build_var_pa_atomized.py
   path-aware records         per-base SNP-level atomization
   (SNP / indel / SV)         (carriers UNIONed over the aligned
   → var_pa_231_arch3_chr1     overlap of every source record)
                              → ..._arch3_chr1_atomized
```

## Per-step

| step | script | what it does | key decisions |
|---|---|---|---|
| A1 | `jobA1_annotate_chr1.sh` | `annotate_vcf` on the 135-sample graph VCF → annotated multi-allelic + biallelic catalog with symbolic `INFO/ID` | the big/slow step (~30 min, ~20 GB); outputs feed both A2 and A3 |
| A2 | `jobA2_pg_chr1.sh` | PG-153 side: subset → `transfer_id` → `convert-to-biallelic` → fill-tags → **haploidize (het→`.`)** | het→missing follows the Arouisse 2020 *A. thaliana* precedent; no V4 filter (per-cell het→missing covers it) |
| A3 | `jobA3_cactus_chr1.sh` | cactus-78 side: subset → `transfer_id` → `convert-to-biallelic` → sort/bgzip | already haploid (long-read assemblies) — no haploidize step |
| A4 | `jobA4_merge_chr1.sh` | merge the two haploid biallelic sides → 231; fill-tags; drop `AN=0` | **`--merge none`**: keep records biallelic — a multi-allelic merge would make `build_var_pa` store `ALT[0]` only but count any non-zero allele as its carrier (mis-attribution) |
| A5 | `jobA5_build_cnvar.sh` | `build_var_pa.py` on the merged panel | → `var_pa_231_arch3_chr1.{var_pa,var_called,meta}.npz` (raw, path-aware; SNP+indel+SV) |
| D1 | `jobD1_atomize_cnvar.sh` | `build_var_pa_atomized.py` on the same merged panel | per-base SNP catalog; carriers unioned over the aligned overlap of every source record; pure INS/DEL beyond the overlap do **not** atomize (they stay in the raw var_pa) |

## Two outputs, two uses
- **Raw** `var_pa_231_arch3_chr1` (A5): one column per biallelic record — SNPs, indels, and
  SVs. Path-aware; use for SV-aware AF and the general projection.
- **Atomized** `var_pa_231_arch3_chr1_atomized` (D1): one column per single-base
  substitution, carriers unioned across every record implying it. Use for SNP-level GEA /
  comparison against linear-reference SNP callers (closed MNP-vs-SNP encoding outliers).

## Prerequisites

These build scripts are SLURM jobs written for the Berkeley Savio cluster; the `#SBATCH`
headers (account/partition/qos) are site-specific. To run elsewhere, adjust those headers
and provide the following:

- **Tools** (the scripts reference absolute conda-env paths — `bcftools`, `bgzip`, `tabix`,
  and a Python with `pysam`, `numpy`, `scipy`): edit the `BCF`/`BGZIP`/`TABIX`/`PY` vars at
  the top of each job, or put equivalents on `$PATH`.
- **External decomposition tools** under `external/` (gitignored, not shipped):
  - `annotate_vcf.py` from HPRC [`eblerjana/genotyping-pipelines`](https://github.com/eblerjana/genotyping-pipelines) (`prepare-vcf-MC`).
  - `convert-to-biallelic.py` from [`eblerjana/pangenie`](https://github.com/eblerjana/pangenie) (`pipelines/run-from-callset/scripts`).
- **Inputs** (live outside the repo; set the paths near the top of jobA1–A3):
  the 135-sample minigraph-cactus pangenome VCF + GFA (A1), and the haploid PG/cactus
  side VCFs from `panel/pangenie_genotyping/` (A2/A3).
- Override `$ARCH3_CHR1_DIR` if launching from an sbatch spool copy outside the source tree.

## Validation

Use **open-loop** truth only (FASTA-lookup / independent estimators). Do **not** validate
against `compute_recomb_truth.py`, which reads the same `var_pa` and so cannot see a
decomposition bug (investigation doc §3, §7).

`jobA7_compare_vs_hapfire.sh` is the kept open-loop check: it joins kMate per-SNP AF
against an **independent** estimator (xwu's hapFIRE on the 1001G SNP catalog) on the
4-tuple `(chrom, pos, ref, alt)` — joining on `pos` alone manufactures off-diagonal
scatter at multi-allelic split records. Its `NEW` input is a per-sample AF TSV (from the
production driver `src/per_sample_per_chrom.py`, or the archived `jobA6`); the hapFIRE
inputs are external. Headline result on SEEDMIX_S1 / Chr1 (518,570 shared SNPs):

| panel | MAE vs hapFIRE | \|Δ\|>0.10 |
|---|---|---|
| v3qc_v3 (`norm -m -any`) | 0.0150 | 0.88% |
| arch3 raw | 0.0140 | 0.64% |
| arch3 atomized | 0.0131 | 0.37% |

Other one-off projection/outlier diagnostics from the Chr1 build were retired to
`chr1/archive/` (see its README).

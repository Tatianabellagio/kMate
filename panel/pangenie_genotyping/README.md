# pangenie_genotyping — building the 231-founder haploid panel VCF

Genotypes the **153 short-read founders** (no long-read assembly) with PanGenie,
QCs and het-masks them. This is the **short-read genotyping side** of the panel:
its per-side VCFs (under `data/v3qc_tmp/`) are the inputs that the `panel/arch3/`
decomposition ingests (A2/A3) and re-merges by symbolic-ID into the production
panel VCF `panel/arch3/chr{N}/merged_231_chr{N}_final.vcf.gz`.

> **Production panel VCF = `panel/arch3/chr{N}/merged_231_chr{N}_final.vcf.gz`**
> (`docs/PIPELINE_STATE.md` §0), NOT the v3qc_v3 outputs below. arch3 does **not**
> consume `founders_231_v3qc_v3.haploid.vcf.gz` — it re-merges the haploid
> biallelic *side* VCFs (cactus-78 + het-masked-haploid PG-153) from
> `data/v3qc_tmp/` via `bcftools merge --merge none`. The old
> `founders_231_v3qc*.vcf.gz` naive-`norm` merge is **archived / superseded** by
> arch3 (§0). The k-mer index used here for *genotyping*
> (`data/pang_135_pangenie_index`) is also distinct from the in-house production
> `kmer_pa` index in `panel/pangenie_index/`.

## Production chain (`scripts/`)

```
build_ena_manifest.py                         data/ena_manifest.tsv (151 missing ecotypes)
   │
launch_downloads.sh → download_one.sh         raw ENA fastqs  (+ concat_lanes.sh for multi-lane)
   │
launch_preprocess.sh → preprocess_one.sh      Trimmomatic + Clumpify → data/preprocessed/
   │
build_pangenie_index.sh                       data/pang_135_pangenie_index  (genotyping graph index)
   │
launch_pangenie.sh → pangenie_one.sh          PanGenie genotype each founder → data/genotyped/*.vcf.gz
   │
merge_vcfs.sh                                 → data/merged/pangenie_151.vcf.gz
   │
qc_pg_v4_filter.sh                            Step 2 merge (151 + 2 LOO) → v3qc_tmp/pangenie_153_raw.vcf.gz
   │
build_v3qc_v2_phase_a.sh                      → v3qc_v2/{cactus_78_bi, pangenie_153_filled_bi}.vcf.gz
   │                                            (v2 Phase B is superseded; v3 reuses ONLY these two Phase-A files)
build_v3qc_v3_phase_a.sh                      het-mask → v3qc_v3/pangenie_153_hetmasked_filled_bi.vcf.gz
   │
haploidize_pg_hetmasked.sh                    → pangenie_153_hetmasked_haploid.vcf.gz
   │
build_v3qc_v3_phase_b.sh                      merge cactus_78_bi + hetmasked_haploid; re-decompose; drop AC=0
   │                                            → v3qc_v3/founders_231_v3qc_v3.vcf.gz
haploidize_v3qc_v3_vcf.sh                     → founders_231_v3qc_v3.haploid.vcf.gz   ← SUPERSEDED naive-norm merge (archived; arch3 re-merges the v3qc_tmp side VCFs instead)
```

The **production** 231-founder panel VCF is built downstream by `panel/arch3/`
(symbolic-ID merge of the `v3qc_tmp/` side VCFs) → `merged_231_chr{N}_final.vcf.gz`.
The `founders_231_v3qc_v3.haploid.vcf.gz` produced by the last step above is the
old single-VCF `norm -m -any` merge, retained for provenance only.

### Cactus side (78 founders) and the 151 → 153 count

The chain above is the **short-read (PanGenie) side**. The **cactus side**
(`data/v3qc_tmp/cactus_78.vcf.gz`, consumed by `build_v3qc_v2_phase_a.sh`) is the
80-assembly minigraph-cactus pangenome VCF subset to drop the two duplicate
assemblies (5772, 9947) → 78 founders. That subset step is preserved in
`scripts/archive/build_v3qc_merged.sh` (Step 1); it is not re-run by the current
chain (the `cactus_78.vcf.gz` it produced is reused as an input).

The PanGenie side genotypes **151** founders (the `pangenie_151` set, which
already includes 100001/100002 via the xwu-BAM path in `preprocess_one.sh`);
`qc_pg_v4_filter.sh` then GQ≥20-masks those and merges in **2** leave-one-out
genotyped founders — **5772 and 9947** (`data/loo_genotyped/{5772,9947}_genotyping.vcf.gz`)
— at its Step-2 merge → the **153** short-read founders. 78 cactus + 153 PG =
the 231-founder panel.

## Validation (kept, off the build chain)

The PanGenie leave-one-out (LOO) concordance harness checks PG-genotyping
quality against cactus assembly truth for the founders that have both. It does
**not** feed the production panel; kept in place as reusable validation:
`download_loo_one.sh`, `launch_loo_downloads.sh`, `preprocess_loo_one.sh`,
`launch_loo_preprocess.sh`, `launch_loo_pangenie.sh`, `pangenie_loo_one.sh`, and
`loo_concordance.py` (with `data/smoke/smoke_concordance.sh` as a tiny smoke test
for the concordance scorer). Superseded LOO variants (`validate_loo.sh`, the
Beagle-imputed `loo_concord_imputed_one.sh`) are in `scripts/archive/`.

## Prerequisites

These are SLURM jobs written for the Berkeley Savio cluster; `#SBATCH` headers are
site-specific. Each script roots its paths at `$BASE`, which defaults to this
directory (resolved from the script location) — **override `$PANGENIE_GT`** when
launching from an sbatch spool copy outside the source tree. Tools
(`bcftools`/`bgzip`/`tabix`, Trimmomatic, Clumpify/BBMap, PanGenie, and a Python
with `pysam`/`numpy`/`scipy`) and the external inputs are referenced by absolute
conda-env / scratch paths near the top of each script; edit those or place
equivalents on `$PATH`. External inputs to supply:
- the cactus pangenome VCF + GFA, and the TAIR10 reference;
- ENA fastqs for the 151 short-read founders;
- the **cactus assembly-ID → 1001G-ID rename map** (`$SAMPLE_RENAME`, required by
  `merge_vcfs.sh` and `pangenie_loo_one.sh` to subset/reheader the 80 cactus
  founders; two columns, `assembly_id`↦`ecotype_id`).

## `scripts/archive/`

Superseded earlier panel generations and concluded one-offs (history preserved):
the v1 (`build_v3qc_merged.sh`) and v2 (`build_v3qc_v2_phase_b.sh`,
`haploidize_v3qc_v2_vcf.sh`) finals, older haploidizers
(`haploidize_v3qc_vcf.sh`, `haploidize_merged_vcf.sh`), the rejected
missingness/MAC-filter experiments (`build_v3qc_{lowmiss,nomiss,nomac_f50}_vcf.sh`,
`count_pg_missingness*.sh`, `count_no_pgmac_missingness.sh`,
`count_mac_on_correct_v2.py`), a one-off v2-index builder
(`build_pangenie_index_v3qc_v2.sh`), and the dead `build_cn_from_pangenie.sh`.

# pangenie_genotyping — building the 231-founder haploid panel VCF

Genotypes the **153 short-read founders** (no long-read assembly) with PanGenie,
QCs and het-masks them, then merges with the **78 long-read cactus founders** to
produce the production haploid panel VCF that feeds both `kmer_pa` (via
`src/build_kmer_pa.py`) and `var_pa` (via the `panel/arch3/` decomposition).

**Production output:** `data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz`
(78 cactus + 153 PG, fully haploid).

> This dir builds the **panel VCF**. The biallelic per-record matrices
> (`var_pa`/`var_called`) are built downstream in `panel/arch3/` — which
> *consumes* this haploid VCF and does **not** rebuild the panel. The k-mer index
> used for *genotyping* (`data/pang_135_pangenie_index`) is distinct from the
> in-house `kmer_pa` k-mer index in `panel/pangenie_index/`.

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
haploidize_v3qc_v3_vcf.sh                     → founders_231_v3qc_v3.haploid.vcf.gz   ← PRODUCTION PANEL
```

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

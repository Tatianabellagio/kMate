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

## Validation (kept, off the build chain)

The PanGenie leave-one-out (LOO) concordance harness checks PG-genotyping
quality against cactus assembly truth for the founders that have both. It does
**not** feed the production panel; kept in place as reusable validation:
`download_loo_one.sh`, `launch_loo_downloads.sh`, `preprocess_loo_one.sh`,
`launch_loo_preprocess.sh`, `launch_loo_pangenie.sh`, `pangenie_loo_one.sh`,
`loo_concordance.py`, `validate_loo.sh`, and `loo_concord_imputed_one.sh` (the
last one targets the hard-rejected Beagle-imputed path — moot, see
`panel/README.md` and `archive/imputation/`).

## `scripts/archive/`

Superseded earlier panel generations and concluded one-offs (history preserved):
the v1 (`build_v3qc_merged.sh`) and v2 (`build_v3qc_v2_phase_b.sh`,
`haploidize_v3qc_v2_vcf.sh`) finals, older haploidizers
(`haploidize_v3qc_vcf.sh`, `haploidize_merged_vcf.sh`), the rejected
missingness/MAC-filter experiments (`build_v3qc_{lowmiss,nomiss,nomac_f50}_vcf.sh`,
`count_pg_missingness*.sh`, `count_no_pgmac_missingness.sh`,
`count_mac_on_correct_v2.py`), a one-off v2-index builder
(`build_pangenie_index_v3qc_v2.sh`), and the dead `build_cn_from_pangenie.sh`.

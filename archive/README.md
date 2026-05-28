# archive/ — superseded work, parked not forgotten

**Archived 2026-05-26** during the production-main cleanup. Nothing here is deleted:
every item is still on disk under `archive/`, and the pre-archive state is captured in
git history (snapshot commit immediately preceding the archive commit — `git log` for the
"snapshot: capture full working-tree state before archive cleanup" commit). Large data
dumps (`*.tsv/*.txt/*.gz/*.jf/*.npz/*.vcf*`) are intentionally git-ignored under `archive/`;
they live on disk only. Code, docs, notebooks, and plots stayed tracked (moved as renames).

This is the "things we tried and dropped" log so future-us doesn't re-litigate them.
Companion archives: `old_docs/` (superseded docs; see `old_docs/METHODS_TRIED.md` for the
full design-space history) and `poolfreq/src/archive/` (superseded solver code).

## What's here and why it was dropped / superseded

| Item | What it was | Why archived |
|---|---|---|
| `imputation/` | Beagle imputation of the SV+SNP panel | **Hard-rejected**: −25 to −30 pp concordance loss on small/medium SVs (LOO). Production keeps PanGenie SV calls unimputed. Never propose Beagle as a fix. |
| `contamination_test/` | AF-based diagnostic-site contamination / mix-in test on SEEDMIX | Planned 2-phase test, never productionized. Uses a **1141-sample diagnostic panel**, *not* the canonical 231-panel SEEDMIX reference — do not reuse its AFs as truth. |
| `panel_overlap_v3_grenenet/` | Overlap of `cn_var_v3` vs GrENE-Net 231 SNPs (the 55%-coverage analysis) | Superseded by the arch3 decomposition (lifted SNP coverage) and by `panel_overlap_135_vs_82/` (current panel-composition analysis). |
| `cactus_panel_overlap/` | Early our-cactus-panel vs xwu-82-accession position comparison | Superseded by `panel_overlap_135_vs_82/`. |
| `pangenome_comparison/` | Record-classification + SV-subgraph topology comparison across graphs/callers | Exploratory; conclusions folded into the cn_var decomposition investigation. |
| `preprocess_qc/` | LD / concordance / missingness QC of the genotyped panel | Conclusions folded into `MISSINGNESS_231PANEL.md` and the v3qc panel-QC decisions. |
| `results/` | Early top-level result tables (e.g. recomb window-size comparison) | Superseded by per-subproject results (`control_p80/results/`, `arch3/chr1/`). |
| `jf_chr1/` | Early jellyfish + founder Jaccard-diagnostic working dir (cactus-vs-1001g heatmaps, carrier-weighted experiments) | Diagnostics folded into the k-mer-imbalance investigation; not a production path. |
| `pangenie_test/` | Early PanGenie test runs + PCA grouping work | Superseded by production `pangenie_genotyping/`. |
| `HapFIRE/` | Local clone of hapFIRE source + early sim outputs | Stale clone; the comparator binary actually used is at xwu's path (see `HANDOFF.md` → Reproducibility). |
| `notebooks/` | Duplicate notebook dir (`blocks_snp_sv_density.ipynb`, `run_plots.py`) | Consolidated into `notebook/`. |
| `logs/` | Old SLURM stdout/stderr (≤2026-04-26) | Ephemeral run logs. |
| `cov10_*_p231_chr1_1_0.output/` (×4) | Stray hapFIRE-format per-window `haplotype_frequencies` outputs from early p231 recomb-sim runs | Superseded by the current `sims/` + `control_p80/` sim framework. |
| `FINAL_RESULTS_cov10_v3_arms_only.ipynb` | Arms-only (centromere-excluded) cov10 v3 eval | Superseded by the full-chromosome `notebook/FINAL_RESULTS_cov10_v3.ipynb`. |

## Added 2026-05-27 — cn_full / cn_var / fasta cleanup (~101 GB)

Front-runner settled (`filt2 + ω=1/m_b`, arch3 cn_var projection) ⇒ experimental zoo archived. Full verdicts in `../METHODS_TRIED_AND_RESULTS.md` §8.

| Item | What it was | Why archived |
|---|---|---|
| `cn_full_superseded/` (56 GB, 32 dirs) | every `cn_full_231_*` except `v3qc_v3` + `v3qc_v3_filt2` | superseded builds + concluded k-mer filter/subsample experiments; front-runner kept in place |
| `cn_var_superseded/` (3.3 GB) | `cn_var_231_{v2,v3,v3qc,v3qc_v2}*` | pre-arch3; production = `arch3/chr1/cn_var_231_arch3_chr1*`, `v3qc_v3` kept |
| `fastas_superseded/` (42 GB) | `unimputed_fastas_v3`, `founder_fastas_231_v3`, `imputed_fastas_v2_DEPRECATED` | current = `unimputed_fastas_v3qc` / `control_p231/fastas_231` |

Some `poolfreq/scripts/{build_cn_full_*,run_seedmix_*}` now reference archived dirs (concluded-experiment scripts; recoverable, kept as records).

## Recovering an archived file into git tracking

```bash
# it's still on disk:        archive/<path>
# its last tracked version:  git log --all --full-history -- <old top-level path>
git show <snapshot-commit>:<old/path> > restored_file
```

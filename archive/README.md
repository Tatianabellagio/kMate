# archive/ — superseded work, parked not forgotten

**Archived 2026-05-26** during the production-main cleanup. Nothing here is deleted:
every item is still on disk under `archive/`, and the pre-archive state is captured in
git history (snapshot commit immediately preceding the archive commit — `git log` for the
"snapshot: capture full working-tree state before archive cleanup" commit). Large data
dumps (`*.tsv/*.txt/*.gz/*.jf/*.npz/*.vcf*`) are intentionally git-ignored under `archive/`;
they live on disk only. Code, docs, notebooks, and plots stayed tracked (moved as renames).

This is the "things we tried and dropped" log so future-us doesn't re-litigate them.
Superseded docs live in `archive/old_docs/` (see `archive/old_docs/METHODS_TRIED.md` for the
full design-space history) — moved here from the repo root 2026-08-24; it was always
described as a companion archive, it just wasn't stored as one.

> **Pointer fixes 2026-08-24.** This file used to reference `poolfreq/src/archive/`
> (superseded solver code) and a singular `notebook/` directory. Neither exists: `poolfreq/`
> was merged into the repo root (`src/ scripts/ tests/ data/`) on 2026-05-28, and the
> notebook directory has always been plural (`notebooks/`). Rows below are corrected.

## What's here and why it was dropped / superseded

| Item | What it was | Why archived |
|---|---|---|
| `imputation/` | Beagle imputation of the SV+SNP panel | **Hard-rejected**: −25 to −30 pp concordance loss on small/medium SVs (LOO). Production keeps PanGenie SV calls unimputed. Never propose Beagle as a fix. |
| `contamination_test/` | AF-based diagnostic-site contamination / mix-in test on SEEDMIX | Planned 2-phase test, never productionized. Uses a **1141-sample diagnostic panel**, *not* the canonical 231-panel SEEDMIX reference — do not reuse its AFs as truth. |
| `panel_overlap_v3_grenenet/` | Overlap of `cn_var_v3` vs GrENE-Net 231 SNPs (the 55%-coverage analysis) | Superseded by the arch3 decomposition (lifted SNP coverage) and by `panel_overlap_135_vs_82/` (current panel-composition analysis). |
| `cactus_panel_overlap/` | Early our-cactus-panel vs xwu-82-accession position comparison | Superseded by `panel_overlap_135_vs_82/`. |
| `pangenome_comparison/` | Record-classification + SV-subgraph topology comparison across graphs/callers | Exploratory; conclusions folded into the cn_var decomposition investigation. |
| `preprocess_qc/` | LD / concordance / missingness QC of the genotyped panel | Conclusions folded into `MISSINGNESS_231PANEL.md` and the v3qc panel-QC decisions. Its SNP/SV LD-decay + SV-SNP tagging analysis (`notebooks/ld_*.ipynb`, `scripts/{compute_ld,aggregate_ld,compute_per_sv_max_r2}.py`, `output/ld/`, `plots/`) and its panel-composition descriptive stats (`notebooks/production_vcf_stats.ipynb`, `scripts/merged_vcf_stats*.{py,sh}`, `output/merged_stats/`) all ran on the **pre-arch3** merged panel (`founders_231_chr.vcf.gz`) — **stale 2026-07-10** now that arch3 is production. Superseded by, respectively: the arch3-vs-GrENE-Net SV tagging comparison in `analysis/grenenet_gea/build_sv_snp_ld.py` → `analysis/grenenet_gea/sv_snp_ld/` (write-up: `analysis/grenenet_gea/notebooks/sv_snp_ld_tagging.ipynb`); and the arch3 panel-composition notebook `analysis/grenenet_gea/notebooks/panel_stats_arch3.ipynb` (backed by `analysis/panel_qc/panel_stats/PANEL_STATS.md` / `scripts/panel_stats_for_paper.py`). Non-LD, non-composition parts of `preprocess_qc/` (dup rates, genotyped concordance) remain live at top level, unarchived. |
| `results/` | Early top-level result tables (e.g. recomb window-size comparison) | Superseded by per-subproject results (`control_p80/results/`, `arch3/chr1/`). |
| `jf_chr1/` | Early jellyfish + founder Jaccard-diagnostic working dir (cactus-vs-1001g heatmaps, carrier-weighted experiments) | Diagnostics folded into the k-mer-imbalance investigation; not a production path. |
| `pangenie_test/` | Early PanGenie test runs + PCA grouping work | Superseded by production `pangenie_genotyping/`. |
| `HapFIRE/` | Local clone of hapFIRE source + early sim outputs | Stale clone; the comparator binary actually used is at xwu's path (see `HANDOFF.md` → Reproducibility). |
| `notebooks/` | Duplicate notebook dir (`blocks_snp_sv_density.ipynb`, `run_plots.py`) | Consolidated into the root `notebooks/`. Note its `blocks_snp_sv_density.ipynb` is the 2026-05-09 version; the live one (2026-07-20) stayed at `notebooks/`. |
| `logs/` | Old SLURM stdout/stderr (≤2026-04-26) | Ephemeral run logs. |
| `cov10_*_p231_chr1_1_0.output/` (×4) | Stray hapFIRE-format per-window `haplotype_frequencies` outputs from early p231 recomb-sim runs | Superseded by the current `sims/` + `control_p80/` sim framework. |
| `FINAL_RESULTS_cov10_v3_arms_only.ipynb` | Arms-only (centromere-excluded) cov10 v3 eval | Superseded by the full-chromosome `FINAL_RESULTS_cov10_v3.ipynb` (now in `notebooks_panel_dev_2026-05/`). |

## Added 2026-08-24 — root-level cleanup (panel-development era)

Companion sweep to the `analysis/` restructure. Code, docs, notebooks and plots stay
tracked and move as git renames. **This round also deleted bulk data for the first time**
— see "Deletions" below; the moves themselves deleted nothing.

| Item | What it was | Why archived |
|---|---|---|
| `notebooks_panel_dev_2026-05/` | The 35 root-level development notebooks (+ their `plots/` subdir) from the panel-construction era: `CN_FULL_COMPOSITION_v3qc_*`, `AF_TRUTH_VS_ESTIMATE_v3qc_*`, `FILT2_RESULT`, `FINAL_RESULTS*`, `RECOMB_SWEEP_RESULTS*`, `V3QC_SANITY_CHECK`, the carrier/bubble diagnostics, and the 2026-06-03 h-imbalance set | All last touched 2026-05-28 / 2026-06-03 and built on the **pre-arch3** `cn_full` / `cn_var` / v3qc panels, which `archive/cn_full_superseded/` and `archive/deprecated_v3qc/` already record as superseded. None was referenced from any doc. **Kept live:** `notebooks/blocks_snp_sv_density.ipynb` (2026-07-20, feeds the LD-block/GEA work). |
| `plots_panel_dev_2026-05/` | The 100 root-level PNGs (`AF_TRUTH_VS_ESTIMATE_*`, `CN_FULL_*`, `H_*`, carrier/bubble diagnostics) | Same era and same subject matter as the notebooks above — these are their rendered figures. Last commit 2026-06-03; not referenced from any doc. |
| `old_docs/` | 20 superseded design/session docs (`ALGORITHM.md`, `METHODS_TRIED.md`, `SESSION_2026-05-*.md`, per-subproject READMEs) | Moved here from the repo root. It was *already* described as a companion archive by this file — it just sat at top level pretending to be live. Its `ALGORITHM.md` is the superseded twin of the root `ALGORITHM.md`; current docs are in `docs/`. |

## Deletions 2026-08-24 (~1.62 TB) — the first bulk removal from `archive/`

Until now `archive/` was append-only: "parked not forgotten". That is still the policy for
**code, docs, notebooks and plots** — none of which were touched. What was removed is
regeneratable output whose *reasoning* is preserved in tracked files.

| Removed | Freed | Why it was safe |
|---|---|---|
| `results_archive_stale/grenenet_kmate_window/` | **1.6 TB** | 23,848 per-sample-per-chrom outputs (13,008 `.tsv` + 10,840 `.npz`) from the **window-mode chain**, which `docs/RERUN_AFTER_FIX.md` records as *"GROUP D — window-mode chain 🅧 RETIRED (decision 2026-07-07, made final 2026-07-10)"*. `analysis/grenenet_gea/GLOBAL_MODE_DECISION.md` is the standing decision that evolved AF is estimated in GLOBAL mode, so this chain will not be rerun. **0 git-tracked files**; every entry was `.tsv`/`.npz` data. The code that produced it survives in `archive/window_hapfreq_retired/`. |
| `archive/contamination_test/` (stale copy) | ~21 GB net | Resolved the root/archive duplication described below. |

### `contamination_test/` — duplication resolved

It existed **twice** (repo root and here) and the two copies had diverged after the
2026-05-26 archiving. Checked before acting:

- every one of the root copy's 236 files was present in the archived copy (0 missing);
- the root copy was **newer on 12 files and older on none**;
- the divergence was a Carnegie→Savio path port (`co_moilab`/`savio4_htc`/`/global/scratch`
  vs `/carnegie/nobackup/scratch/xwu`/`/home/tbellagio`);
- all **24 git-tracked** paths under the archived copy existed in the root copy.

So the **root copy was authoritative and the archived copy was the stale one** — the
opposite of what the directory sizes suggest at a glance. The stale copy was deleted and
the root copy moved into its place, so `archive/contamination_test/` is now the
Savio-ported version (48 files still carrying Carnegie paths, down from 61) and the root
directory is gone. The 24 tracked paths now hold their newer content — expect them to show
as *modified*, not deleted, in git.

Only unique content lost: one `.ipynb_checkpoints` sidecar.

### Remaining bulk (audited 2026-08-24, not removed)

Same character as the above — regeneratable output of superseded panels — but left in place
pending a decision:

| Item | Size | Note |
|---|---|---|
| `cn_full_superseded/` | 220 G | 12+ `cn_full_231_*` panel variants (v2/v3/v3qc/filt2/subsampMedian seeds), 3–27 G each. Pre-arch3. |
| `grenenet_kmate_arch3_oldpanel_archive/` | 154 G | ~340 per-sample AF `.tsv`, uniformly 449 M. Superseded by the current arch3 run. |
| `deprecated_v3qc/` | 112 G | `pangenie_genotyping_data` 62 G, `data` 35 G, `panel_naive_merges` 17 G. |
| `exploration/`, `fastas_superseded/`, `pangenie_test/`, `jf_chr1/` | 199 G | Already documented as exploratory/superseded above. |

## Added 2026-05-27 — cn_full / cn_var / fasta cleanup (~101 GB)

Front-runner settled (`filt2 + ω=1/m_b`, arch3 cn_var projection) ⇒ experimental zoo archived. Full verdicts in `../METHODS_TRIED_AND_RESULTS.md` §8.

| Item | What it was | Why archived |
|---|---|---|
| `cn_full_superseded/` (56 GB, 32 dirs) | every `cn_full_231_*` except `v3qc_v3` + `v3qc_v3_filt2` | superseded builds + concluded k-mer filter/subsample experiments; front-runner kept in place |
| `cn_var_superseded/` (3.3 GB) | `cn_var_231_{v2,v3,v3qc,v3qc_v2}*` | pre-arch3; production = `arch3/chr1/cn_var_231_arch3_chr1*`, `v3qc_v3` kept |
| `fastas_superseded/` (42 GB) | `unimputed_fastas_v3`, `founder_fastas_231_v3`, `imputed_fastas_v2_DEPRECATED` | current = `unimputed_fastas_v3qc` / `control_p231/fastas_231` |

Some `poolfreq/scripts/{build_cn_full_*,run_seedmix_*}` now reference archived dirs (concluded-experiment scripts; recoverable, kept as records).

## Added 2026-07-10 — results/ cleanup: consolidate scattered stale/archived result trees (~1.8 TB)

Part of the `analysis/` vs `results/` centralization pass (kMate results moving under
`analysis/grenenet_gea/`). These three were already flagged stale/superseded and living
under ad-hoc `results/*archive*` names; moved here so there is exactly one archive
convention, not two.

| Item | What it was | Why archived |
|---|---|---|
| `results_archive_stale/` (1.6 TB) | Formerly `results/archive/` — mostly the STALE fixed-bp window-mode cohort (`grenenet_kmate_window` / `_seedmix` / `_smoke`, old multinomial EM, superseded 2026-07-07 by the per_founder fix) plus assorted early seedmix/site04/hapfire-comparison result dumps (`seedmix_231*`, `site04_*`, `session_summary`, `overhang_timing`, etc.). See `docs/PIPELINE_STATE.md` §7/§8 and `docs/RERUN_AFTER_FIX.md` Group D for the staleness call. |
| `grenenet_kmate_arch3_oldpanel_archive/` (154 GB) | Formerly `results/grenenet_kmate_arch3_oldpanel_archive/` — GrENE-Net cohort AF on the old 10.33M (pre-segregating-filter) panel, superseded by the 8.49M-panel rerun. | Superseded panel version; see `analysis/grenenet_gea/WINDOW_UNIT_VALIDATION.md` §3. |
| `seedmix_kmate_arch3_oldpanel_archive/` (7.1 GB) | Formerly `results/seedmix_kmate_arch3_oldpanel_archive/` — SEEDMIX p0 reps on the same old panel. | Same as above. |

## Recovering an archived file into git tracking

```bash
# it's still on disk:        archive/<path>
# its last tracked version:  git log --all --full-history -- <old top-level path>
git show <snapshot-commit>:<old/path> > restored_file
```

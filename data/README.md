# data/ — small tracked inputs (and some large git-ignored ones)

Mixed directory: a handful of **small, tracked, hand-curated inputs** that the pipeline
depends on, plus several **large git-ignored working artifacts** that happen to live here.
The tracked files are the ones that matter — they encode decisions that cannot be
regenerated from code.

Written 2026-08-24 during the root-level cleanup.

## Tracked — curated inputs (these encode decisions; don't regenerate blindly)

**Founder / panel membership**

| File | What it is |
|---|---|
| `vcf_samples_231.txt` | The canonical 231-founder sample list for the production panel. |
| `flag_list.tsv` | Per-assembly QC verdicts with **written reasons** (e.g. mislabeled cactus FASTA identified by k-mer Jaccard, PanGenie substitutes). The audit trail for why founders were excluded or swapped — read this before questioning panel membership. |
| `exclude_list.txt` | Assembly IDs excluded from the panel. |
| `qc_lowcov_exclude.txt` | Samples with Chr1 usable-panel-k-mer fraction < 0.10, including 5 dead/contaminated libraries at <0.01. Excluded from all downstream analysis. |
| `founder_split_cactus_pg.json` | Which founders come from the cactus assemblies vs PanGenie genotyping. |
| `sv_panel_to_accession_id.tsv` | SV-panel ID → accession ID mapping. |

**Cohort manifests** (sample → path lists for SLURM arrays)

`sample_manifest.tsv`, `sample_manifest_savio.tsv`, `sample_manifest_usesample.tsv`,
`sample_manifest_pilot_s4_s54.tsv`. Some sibling manifests
(`sample_manifest_remaining.tsv`, `sample_manifest_timing4.tsv`, `site04_manifest.tsv`) are
git-ignored as regeneratable — see the `data/sample_manifest_*.tsv` rules in `.gitignore`.

**SEEDMIX**

`seedmix_manifest_arch3.tsv` (arch3-era manifest), `seedmix_recipe_normalized.tsv` (the
normalized mixing recipe — the truth the SEEDMIX validation is scored against).

**Block partitions** (used by the GEA/LD-block work)

`greneNet_final_v1.1_density0.5_genomewide_partition.txt`,
`s1_density0.5_fine_genomewide_partition.txt`.

**External reference**

`arapheno.csv` (AraPheno phenotypes), `ASSEMBLIES_Best_version_of_dataset.xlsx` and
`request_assemblies_for_Moi.csv` (assembly provenance correspondence).

## Not tracked — large working artifacts living here

Git-ignored via the `data/kmer_pa_*/`, `data/sim_chr1*/`, `*.npz` rules:

| Path | Size | What it is |
|---|---|---|
| `kmer_pa_231_arch3_filt2inv/` | 7.5 G | Production arch3 `K_pa` matrices (filt2 + invariant-filtered). Regeneratable via `scripts/build_kmer_pa_arch3.sh`. |
| `sim_chr1/`, `sim_chr1_skewed/` | ~750 M | Chr1 simulation fixtures. |
| `hapfire_block_index.npz` | 22 M | hapFIRE block index. |
| `test_chr1_first200.*.npz` | ~2 M | Small test fixtures. |
| `calibration_231_seedmix_s1.json`, `seedmix_S1_h_wls_chr1_first200.tsv` | small | Calibration/working outputs, explicitly re-ignored in `.gitignore`. |

If you add a file here, decide deliberately whether it is a *curated input* (track it, and
add a row above) or a *derived artifact* (add a `.gitignore` rule). The repo's standing
policy is code/docs/small-config tracked, bulk data never.

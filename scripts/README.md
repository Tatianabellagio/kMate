# scripts/ — panel construction, panel auditing, and run accounting

Operational scripts that sit *outside* the `kmate` package (`src/kmate/`) and outside any
one analysis subproject. Rule of thumb: if it builds or audits the **production panel**, or
reports on a **SLURM run**, it lives here. Analysis code belongs in `analysis/<topic>/`;
library code belongs in `src/kmate/`.

Written 2026-08-24 during the root-level cleanup — this directory had been committed to for
months with nothing explaining it.

## Panel construction

| Script | What it does |
|---|---|
| `build_kmer_pa_arch3.sh` | SLURM job that builds the production `K_pa` (k-mer presence/absence) matrices for the arch3 231-founder panel. |
| `apply_monomorphic_filter.py` | Regenerates the production panel as **segregating-only**, dropping monomorphic records (`AC=0` — ALT in no founder — or `AC=AN` — ALT in every called founder). Mirrors the patched jobA4 Step-3 filter `bcftools view -e 'AC=0 \|\| AC=AN'`. |

## Panel auditing

| Script | What it does |
|---|---|
| `invariant_check.py` | Invariant-site audit across VCF + `V_pa` + `K_pa`. Matrix level is definitive, since that's the actual genotype matrix. |
| `compare_index_pg_vs_ours.py` | Level-B equivalence check: our `K_pa` vs PanGenie's index, both built from the same arch3 `merged_231` VCF + REF and filtered identically. Answers whether our indexer is downstream-equivalent to PanGenie's. |
| `compare_kmer_pa_dirs.py` | The generic form of the above — same metrics (k-mer set overlap, per-founder count correlation, carrier agreement on shared k-mers) between any two `kmer_pa` dirs. |
| `nocap_chr1_canary.sh` | Chr1 canary run for the no-cap k-mer indexing regime. |
| `panel_stats_for_paper.py` | Per-chrom and genome-wide panel + kMate-input statistics for the manuscript. Backs `analysis/panel_qc/panel_stats/PANEL_STATS.md`. |

## Validation and accounting

| Script | What it does |
|---|---|
| `validate_seedmix_vs_hapfire.py` | Validates kMate SEEDMIX AF estimates against xwu's hapFIRE SNP-frequency truth, joining the kMate per-record TSV to the per-sample truth. |
| `pilot_compute_stats.py` | Pulls `sacct` accounting for a SLURM array (e.g. the kMate pilot), aggregates the main/`.batch` job span, and reports achieved concurrency/throughput plus a projection to the full cohort. |

## Subdirectories

- `archive/` — superseded scripts, kept per the repo-wide policy in `archive/README.md`.
- `logs/` — SLURM stdout/stderr; git-ignored, prunable.

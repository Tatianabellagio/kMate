# qc — quality control of this analysis

QC of the GrENE-Net selection analysis itself: is the AF input trustworthy, is
the panel what we think it is, and where does kMate lose founders?

Not to be confused with `analysis/panel_qc/`, which is kMate **method/panel**
validation (k-mer index comparison, panel stats, seed-mix duplicate rate).

| script | what it checks | notebook |
|---|---|---|
| `_build_qc_notebook.py` | cohort coverage vs usable-panel-k-mer fraction. Found depth alone is a poor QC signal (corr 0.57), so the cohort is filtered on `nzfrac < 0.10` (17 dropped) rather than depth; 5 libraries look contaminated/mislabelled rather than merely shallow | `qc_coverage_audit.ipynb` |
| `plot_coverage_distribution.py` | sequencing-depth distribution across the 2,415 evolved pool-seq samples | — |
| `_build_support_nb.py` | how many SVs / how much pool mass is lost at each kMate `info` / `n_called` support threshold | `06_sv_support_filter.ipynb` |
| `_compute_panel_overlap_grenenet.py` + `_build_panel_overlap_grenenet_nb.py` | arch3 panel vs the old GrENE-Net 231 SNP catalog — 76.5% exact match (up from 55.05% before the decomposition fix) | `panel_overlap_grenenet.ipynb` |
| `_build_panel_stats_nb.py` | panel composition: class, SV size, AF spectrum, missingness, carriers, density | `panel_stats_arch3.ipynb` |
| `_recreate_density_snp_only.py` | redraws the genomic-density figure SNP-only (and a 2-panel SNP vs SNP+indel+SV version) for the "what SNP-only studies miss" slide. Mirrors `_build_panel_stats_nb.py` section 7 | — |
| `seedmix_identifiability.py` | why kMate under-calls ~19–49 seed-mix founders where hapFIRE loses ~1/231 | — |
| `seedmix_kmer_identifiability.py` | the mechanism: are the under-called founders exactly those with the fewest private k-mers? | — |
| `_build_9977_nb.py` | founder 9977, the worst absorbed founder: no single k-mer twin — diffusely non-identifiable | `founder_9977_kmer_space.ipynb` |
| `compare_window_vs_global_af.py` | per-sample kMate AF, window-mode vs global-mode | — |

## seedmix_validation/

The founder-collapse investigation and its fix arms. `HANDOFF.md` is the
decision trail. **`fix_seedmix/` and `fix_norm/` read as abandoned arms but are
not removable** — later arms (`fix_kfw_fullpanel/`) still read their cached
counts and import `em_variants` from them. The genuinely abandoned arms
(`fix_emreg/`, `fix_local/`) were retired 2026-08-24.

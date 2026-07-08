# RETIRED: SV-haplotype-frequency / Pipeline-B (hapfreq) chain

**Retired 2026-07-08.** Do NOT re-run anything in this directory.

## What this was

The haplotype-frequency ("hapfreq") generator chain, a.k.a. **Pipeline B** — a
pool-seq TEMPORAL climate GEA run on r2>=0.9 LD haploblocks. It projected the
per-sample founder frequency vector `h` onto a founder->hap-cluster membership
to build a per-sample haplotype-frequency matrix, then ran per-block/per-hap
temporal GEA + block-WZA, and fed two downstream SV-haplotype analyses.

Generator chain (orchestrated by `gea_clq90_pipelineB.sbatch`):

    build_hapfreq_matrix.py      # project founder h -> per-sample hap-freq matrix
    build_hapfreq_p0_seedmix.py  # gen-0 seedmix anchor (p0/v0)
    build_hap_trajectories.py    # varlen per-(gen,plot) trajectories
    build_hap_gea.py             # per-hap / per-block temporal GEA (Sidak)
    build_hap_wza.py             # block-WZA
    hapfreq.sbatch               # older 2168-sample window-mode array driver

Phase-1-replication variant:

    phase1_replication/build_hap_lastgen_matrix.py  (+ build_hap_lastgen.sbatch)

Dedicated downstream consumers of the `hapfreq_clq90/pipelineB_varlen` product
(the null SV-haplotype analyses and their plotting/notebook builders):

    sv_adaptive/sv_enrichment_gea.py
    sv_adaptive/sv_temporal_markermatched.py
    _build_sv_selection_nb.py        -> notebooks/sv_selection_currency.ipynb
    _build_svenrich_nb.py            -> notebooks/sv_enrichment.ipynb
    top_block_genes.py               (uses resolve_symbols.py)
    resolve_symbols.py
    site_replicate_concordance.py
    plot_site_scoef.py
    plot_site_scoef_mixed.py
    plot_varlen_manhattan.py
    notebooks/sv_enrichment.ipynb
    notebooks/sv_selection_currency.ipynb

Product data:

    hapfreq_clq90/   # 3.5G pre-fix product of this chain (moved here from
                     # results/grenenet_gea/hapfreq_clq90/)

## Why retired

1. **BROKEN INPUT.** `build_hapfreq_matrix.py` reads
   `results/grenenet_kmate_window/` (WIN store), a **deleted** window-mode store.
   Any re-run reads nothing / stale data and produces garbage.
2. **Both SV threads it fed are RESOLVED NULL.**
   - *SV-selection-currency*: the "SVs are a fitness currency" (x2.07) signal was
     a block-vs-haploblock **unit-mismatch artifact**; at the correct hap-cluster
     unit selected haplotypes are NOT SV-enriched (n.s.).
   - *SV-temporal-purging*: died on the corrected `h`.
3. **Superseded framing.** The production cohort is now `--unit chrom` (global-h),
   which collapses to ~231 haplotypes genome-wide and produces no per-block
   `h_blocks` for this chain to consume. Pipeline B's premise no longer holds.

## LEFT IN PLACE (do NOT confuse with this retirement)

- **build_hap_membership.py** + `hap_membership*.sbatch` +
  `_build_sv_selection_audit_nb.py` + `notebooks/sv_selection_haplotype_audit.ipynb`
  — a SEPARATE, KEPT SV-selection **audit** uses `build_hap_membership`. The
  audit also owns `_sv_hap_rotationnull.py`, `_sv_hap_context.py`,
  `_sv_hap_freqrobust.py`, `_sv_haplotype_axes_sweep.py` (+ their `.sbatch`
  drivers), so those were LEFT even though they carry `sv_hap` names.
- **baypass_build_inputs.py** + `BAYPASS_TEMPORAL_PLAN.md` — a not-yet-run plan.
  It still reads `results/grenenet_gea/hapfreq/hapfreq_matrix.npy` and
  `hapfreq_registry.csv`, which are now **stale/retired** products of this chain.
  Because a kept file references them, those two files were NOT moved; the plan
  should be revisited before use.

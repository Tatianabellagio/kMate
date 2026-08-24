# extras — side investigations, not one of the three results

Everything here asks **one** question in many different ways:

> Are SVs on the fitness-selected haplotypes — or are they passengers riding
> along with the haplotype that actually carries the signal?

It is kept out of r1/r2/r3 because the answer came out **null**, and because
the question is a control on those results rather than a result itself.

## The verdict

- Apparent enrichment of SVs on selected haplotypes (~1.4–2.0×) **collapses
  under proper nulls** — founder-count-matched, regional, and genome-rotation.
- In fine-mapping, **SVs never lead a credible set**: 0/22 blocks at both honest
  (N=31 sites) and inflated (N=355 pools) power; max SV-PIP anywhere 0.165.
- Apparent top-block SV enrichment was a **block-size artifact** (SV-containing
  blocks median 73 records vs 7).

Read it in `../notebooks/sv_selection_haplotype_audit.ipynb`.

## The arms

| script | the version of the question it asks |
|---|---|
| `ecotype_hap_sv_enrichment.py` | the core, at the **correct unit**: are fitness-winning hap-clusters SV-tagged? (founder-count-matched null) |
| `_sv_haplotype_axes_sweep.py` | same, across all three selection axes (JOINT / CLIMATE / TEMPORAL) + a MAF sweep + power |
| `ecotype_hap_sv_regional.py` | is the enrichment a regional confound? |
| `ecotype_hap_sv_rotation2.py` | a null controlling frequency **and** space at once (v1 lost the founder-count control and was retired) |
| `_sv_hap_rotationnull.py` | the spatial block-rotation null |
| `_sv_hap_freqrobust.py` | is it a frequency-alignment artifact? |
| `_sv_hap_context.py` | is it a regional confound at the haplotype unit? |
| `_sv_passenger_test.py` | decisive: in top-JOINT SV blocks, does the SV or the **haplotype** carry the founder-fitness signal? |
| `_sv_founder_direction.py`, `_sv_founder_mechanism.py` | the founder-GWAS is direction-agnostic; mechanism check on the ×2.07 top-0.5% signal |
| `_sv_winning_genetics.py` | are SVs part of the shared winning genetics across sites? |
| `_audit_sv_fitness.py` | audits the original "common SVs enriched on fitness-selected haplotypes" claim |
| `_build_sv_selection_audit_nb.py` | builds the audited notebook |
| `_winners_ratio.py`, `_winners_sv_depletion.py`, `_plot_winner_ratio_allsites.py`, `_plot_winner_ratio_vs_purging.py`, `_plot_winners_sv_depletion.py` | the **panel-artifact control on result 1**: if short-read (PanGenie) founders carry fewer called SVs and those founders rise at hot sites, "SV purging" would be an artifact of panel completeness rather than selection |
| `driver_passenger/` | SuSiE-RSS fine-mapping within climate-selected blocks — the driver-vs-passenger test at the fine-mapping unit |

`ecotype_block_sv_enrichment.py` (the block-unit first attempt) was retired
2026-08-24: its own successor states the block test used the wrong unit.

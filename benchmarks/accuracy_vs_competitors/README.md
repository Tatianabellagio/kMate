# benchmarks/accuracy_vs_competitors/

kMate against other tools. Design rationale and the full argument live in
[`BENCHMARK_DESIGN.md`](BENCHMARK_DESIGN.md) — **read its top banner first**, it
records which half of this directory is retired.

## Current state — one live comparison, one retired

| comparator | status | where |
|---|---|---|
| **vg giraffe** | ✅ **LIVE** | `results/plots/vg_vs_kmate_p231_snp_accuracy.png`, driven by `scripts/score_vg_vs_kmate_p231.sh` |
| **hapFIRE** | ⛔ **RETIRED here** | moved to `../archive/shared_panel_hapfire_comparison_invalid/`; the live hapFIRE comparison is `../speed_vs_hapfire/` |

### Why the hapFIRE half is retired

It scored kMate and hapFIRE on **one shared SNP panel**. hapFIRE errors on missing
GTs, so `scripts/build_shared_snp_panel.sh` builds that panel with
`MISSING '.' -> 0|0` — handing hapFIRE false REF homozygotes at ~7.1% of records,
while kMate's `var_called` mask excludes exactly those sites. A shared panel is not a
shared **input**. The live replacement gives each tool its native panel.

### Why vg is not retired with it

vg's AF comes from **read coverage over graph nodes**, not from inference over panel
genotypes, so the imputed GTs do not enter its estimate the way they enter hapFIRE's
model. Caveat, unquantified and stated rather than assumed: vg's graph *is* built by
`vg autoindex` from the same shared SNP VCF (`build_giraffe_graph_p231.sbatch:25`), so
missing→REF does reach its haplotype paths and could bias mapping. Measure this before
headlining the vg arm.

## Still valid in BENCHMARK_DESIGN

- **§3.1** — hapFIRE/HARP structurally cannot estimate SVs. An argument about HARP's
  per-base substitution likelihood; panel-independent.
- **§4** — speed. Matched-allocation wall/CPU time (maintained in `../speed_vs_hapfire/`).

## Scripts

    build_shared_snp_panel.sh / build_shared_sv_panel.sh   shared site list
    build_truth.py / build_truth_sv.py                     tool-independent truth
    score_snp_fair.py                                      allele-safe SNP scorer (used by the vg path)
    score_sv.py, _score_one_vgsv.py                        SV scoring
    build_giraffe_graph_p231*.sbatch, run_vg_*.sbatch      vg arm (live)
    giraffe_*.sbatch                                       giraffe tuning probes (kept: the
                                                           rescue-off / concurrency findings are
                                                           load-bearing for the vg runtime story)
    run_hapfire*.{sh,sbatch}, align_pool.sh                hapFIRE arm (retired comparison)
    score_competitors.py, ../score_competitor_snp.py       KNOWN-WRONG, hard-guarded

The two KNOWN-WRONG scorers `sys.exit` unless `KMATE_ALLOW_BAD_SCORER=1`; they join on
`(chrom,pos)` only and mispair alleles at multiallelic sites (`../SCORING_RULES.md`
RULE 1). They are left in place, not archived, because their callers still drive vg —
the guard is the quarantine.

## Layout

    scripts/   runners + scorers
    work/      graphs, indexes, BAMs, shared panels (large, gitignored)
    results/   scored tables; results/plots/ for figures
    logs/      SLURM logs (gitignored)

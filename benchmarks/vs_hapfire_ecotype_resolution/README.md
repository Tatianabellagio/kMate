# benchmarks/vs_hapfire_ecotype_resolution/

**Question:** how many of the founders actually in a pool can be *resolved* — named,
not just counted — and how cleanly does the estimated mixture concentrate on them?

This is the **set-recovery / resolution** benchmark. It is deliberately distinct from
`../vs_hapfire/`, which measures founder-frequency *error* (h R², h RMSE) and
runtime on the same native-panel basis. The metrics here — `eff_n`, `recall`,
`precision`, `mass_on_true`, `top_n_jaccard` — are not computed there.

| metric | meaning |
|---|---|
| `true_n` | founders actually in the pool (`weight > 0`) |
| `eff_n` | effective #founders, `1 / Σh²` — kMate's headline readout |
| `n_detected` | founders with `h ≥ 0.5/true_n` (half the equal share) |
| `recall` | true founders detected / `true_n` |
| `precision` | detected that are true / `n_detected` |
| `mass_on_true` | fraction of estimated mass landing on true founders (1 = perfect) |
| `top_n_jaccard` | Jaccard(top-`true_n` by h, true set) — threshold-free |
| `h_rmse` | RMSE(estimate, true weights) over all panel founders |

## The kMate-vs-hapFIRE comparison — use the NATIVE-PANEL one

| | script | output |
|---|---|---|
| ✅ **current** | `scripts/score_ecotype_nativepanel.py` → `scripts/plot_ecotype_nativepanel.py` | `results/ecotype_nativepanel_table.tsv`, `results/plots/ecotype_nativepanel_resolution.png` |
| ⛔ **retired** | shared-panel arm | `../archive/shared_panel_hapfire_comparison_invalid/vs_hapfire_ecotype_resolution_sharedpanel/` |

The retired arm fed hapFIRE the shared panel
`vs_vg_giraffe/work/shared_snps_<panel>_Chr1.vcf.gz`. hapFIRE errors on
missing GTs, so that panel is built with `MISSING '.' -> 0|0` (impute REF) — handing
hapFIRE false REF homozygotes at ~7.1% of records, while kMate's `var_called` mask
excludes exactly those sites. A shared panel is not a shared *input*. Same class of
error as the 2026-07-15 panel mismatch (`../vs_hapfire/PANEL_MISMATCH_BUG.md`).

The current arm gives each tool its own panel and reuses runs that already exist —
**no new simulation or estimator runs**:

- kMate: arch3 panel/reads, `p231/results/kmate_chrom_poolsize_depth/`
- hapFIRE: greneNet panel/reads, `vs_hapfire/results/greneNet_fair/`
- truth: `vs_hapfire/sims_greneNet/*/pool_weights.tsv`

Truth is legitimately **shared** even though the panels are not: a founder's pool
weight is a property of the draw (matched founder set + seed), not of a
variant-calling lineage. This is exactly why the *resolution* question survives the
native-panel split while a shared per-site *AF* truth does not.

Grid: N ∈ {2,5,20,50,150,231} × depth {1,10}, plus N=50 × depth {30,50}, × seeds
42–46 = **70 matched pools**, all present.

### What it shows

Mixed, and worth stating in both directions:

- **kMate concentrates mass far better.** `mass_on_true` 0.91–0.93 vs hapFIRE
  0.74–0.84 across N=20–150 — hapFIRE spreads ~20% of its mass onto off-pool
  founders where kMate spreads ~8%.
- **hapFIRE has better recall at 1× for large N.** At N=150/231, depth 1×, kMate
  recall falls to ~0.85–0.88 while hapFIRE holds ~0.99–1.00. At 10× the gap closes.
- `top_n_jaccard` and `h_rmse` are close between the tools throughout.

## kMate-only work in this dir (unaffected by the above)

`run_h_only.sbatch`, `score_ecotype.py` (still the per-pool kMate scorer),
`bench_block_kmer_floor.py` / `plot_block_kmer_floor.py` /
`run_block_kmer_floor{,_sweep}.sbatch`, `run_localonly_smoke.sbatch`. These score
kMate against simulation truth with no hapFIRE involved.

## Layout

    scripts/    runners + scorers (native-panel comparison, kMate-only h/floor work)
    results/    ecotype_nativepanel_table.tsv, per-pool h npz, ecotype_rows/
    results/plots/   figures (project convention: figures never loose in results/)

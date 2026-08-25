# shared-panel kMate-vs-hapFIRE comparison — ARCHIVED AS INVALID (2026-08-25)

**Do not cite any number in here.** This was the "one shared SNP panel, one shared
truth, scored identically for every tool" design (BENCHMARK_DESIGN §2). It is not a
fair comparison against hapFIRE, because **kMate and hapFIRE do not take the same
kind of input**, and forcing them onto one panel silently handicaps hapFIRE.

## The mechanism

`build_shared_snp_panel.sh` line 10-12 states it plainly:

> haploid founder GT -> phased homozygous diploid (inbred ecotypes): 0->0|0,
> 1->1|1; **MISSING '.' -> 0|0 (impute REF)** so hapFIRE sees a complete phased
> panel (it errors on any unphased/missing GT).

So every missing founder genotype in the shared panel is handed to hapFIRE as a
**false REF homozygote**. kMate is not exposed to this: its `var_pa` ships a
companion `var_called` mask, and its MAR projection `(h·var_pa)/(h·var_called)`
divides those sites out. One tool gets to say "unknown", the other is forced to
guess REF. That asymmetry is baked into the panel, not into the algorithms.

On the p80 shared panel that is 7.1% of records carrying any missing GT
(mean F_MISSING 0.0036, and individual records as high as 0.325).

Note the handicap is **missingness/imputation, not phasing** — a plausible-sounding
but wrong diagnosis. The shared panel is fully phased (`0|0`, 160,000/160,000 phased
separators over the first 2,000 records); haploid→phased-homozygous is exactly what
that builder does. If you revisit this, do not "fix" it by re-phasing.

The design doc half-knew this: §2.1 documents that the winner **flips entirely** with
the truth convention (kMate R² 0.998→0.874, hapFIRE 0.880→0.997) and set a rule to
headline only the physical + fully-called basis. That rule treated the symptom. The
cause is that the shared panel is not a shared input in the first place.

## What replaces it

**`benchmarks/speed_vs_hapfire/`** — matched pools, tool-native panels. Both tools get
the same founder draw and seed, and each is scored against its own natural truth
(kMate on arch3 SNPs, hapFIRE on greneNet SNPs projected through the same
`pool_weights.tsv`). Its scorer says so explicitly: *"This is NOT a shared per-site
truth ... It IS the fair, matched-pool comparison."* It also carries 5 seeds ×
{2,5,20,50,150,231} founders × {1,5,10}× coverage, where this archived design had one
pool at one seed.

Direction differs, which is the whole point — at n231/cov10 the own-panel comparison
puts hapFIRE ahead (RMSE 0.0018 vs kMate 0.0048) while this archived shared-panel one
put kMate ahead (0.0034 vs 0.0058). kMate's own value is stable across both
(R² ≈ 0.9993); what moves is hapFIRE, precisely because the shared panel injects the
false-REF sites its native pipeline never sees.

## What is NOT invalidated by this

- **§3.1 — hapFIRE/HARP cannot estimate SVs.** A structural argument about HARP's
  per-base substitution likelihood. Panel-independent; still stands.
- **§4 — speed.** Wall-clock and CPU-time on matched allocations; does not depend on
  the shared panel. Still stands (`speed_vs_hapfire/`).
- **The kMate estimator refresh itself.** ROADMAP_GLOBAL_REFRESH Phases 1-4 score
  kMate against simulation truth, with no hapFIRE involved.
- **The kMate-vs-vg comparison — LIVE, and deliberately not archived.** The maintained
  version is `accuracy_vs_competitors/results/plots/vg_vs_kmate_p231_snp_accuracy.png`
  (p231, 5 seeds × {2,5,20,50,150,231} founders × {10,50}×), driven by
  `scripts/score_vg_vs_kmate_p231.sh`. Only the *hapFIRE* half of the old 4-tool
  apparatus is retired here.

  **`score_snp_fair.py` is therefore NOT archived** — it was briefly moved here on
  2026-08-25 and restored the same day, because the live p231 vg pipeline calls it. The
  scorer was never the defect: its shared-panel restriction exists to make `(chrom,pos)`
  a unique, allele-safe join key (SCORING_RULES RULE 1), which is correct and still
  needed. The defect was feeding that panel to hapFIRE *as its model input*.

  ⚠️ Open question, stated rather than assumed: vg is not fully immune either. Its
  graph is built by `vg autoindex --workflow giraffe` from this same shared SNP VCF
  (`build_giraffe_graph_p231.sbatch:25`), so the `MISSING -> 0|0` imputation does reach
  vg's haplotype paths. It is much weaker than hapFIRE's exposure — vg's AF comes from
  read coverage over graph nodes, not from inference over panel genotypes, so imputed
  GTs bias mapping rather than the estimate itself. Not quantified. If the vg arm is
  ever headlined, measure this rather than citing this paragraph.

## Contents

    results/benchmark_table_4tool.tsv              pre-refresh table (2026-06-22)
    results/benchmark_table_4tool_chrom.tsv        refreshed kMate arm (2026-08-25)
    results/benchmark_4tool_RMSE_*.png             pre-refresh figures
    results/benchmark_4tool_chrom_RMSE_*.png       refreshed figures
    results/snp_fair_chrom_refresh.tsv             refreshed SNP scores
    scripts/score_all_competitors{,_chrom}.sh      orchestration
    scripts/build_4tool_{table,figure}.py          table + figure builders

(`score_snp_fair.py` is NOT here — see above; it lives at
`accuracy_vs_competitors/scripts/` because the live vg pipeline calls it.)

Historical note: the `_chrom` artifacts are the 2026-08-25 Phase-5 refresh, which
re-ran the kMate arm under `--unit chrom` + per-founder normalization. That refresh
was executed correctly (competitor rows came out byte-identical, confirming only the
kMate arm moved) — it is archived because the *design* is invalid, not because the
run was wrong.

---

## `ecotype_count_sharedpanel/` (added 2026-08-25)

The second site of the same defect, found during the benchmarks cleanup.
`ecotype_count/scripts/run_hapfire_ecotype.sbatch:23` fed hapFIRE the same
`shared_snps_<panel>_Chr1.vcf.gz`, so its kMate-vs-hapFIRE ecotype figures carried the
identical missing→REF handicap.

Archived here: `run_hapfire_ecotype.sbatch`, `ecotype_kmate_vs_hapfire.png`,
`ecotype_radar_vs_hapfire.png`, and the 12 `*_hapfire.tsv` per-pool rows.

**Replaced, not merely retired** — unlike the SNP-AF case, the ecotype *question*
survives the native-panel split, because its truth (`pool_weights.tsv`) is a property
of the founder draw rather than of a variant-calling lineage. The replacement is
`ecotype_count/scripts/score_ecotype_nativepanel.py` (+ its plot script): kMate on
arch3, hapFIRE on greneNet, matched founder draw + seed, 70 matched pools across
N ∈ {2,5,20,50,150,231} × depth {1,10} (+ N=50 × {30,50}) × seeds 42–46. It reuses
existing runs — no new simulation or estimator jobs.

kMate's own ecotype work in that dir (h-only, block-k-mer-floor, localonly smoke) was
never affected and stayed live.

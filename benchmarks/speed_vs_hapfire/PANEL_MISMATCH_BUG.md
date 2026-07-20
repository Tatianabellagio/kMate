# ⛔ FATAL BENCHMARK BUG: kMate and hapFIRE were run on DIFFERENT panels

**Date found:** 2026-07-15. **Status: FIXED, 2026-07-15.** This invalidated all
`hapfire_vs_kmate_*` accuracy results produced before the fix (now archived
under `archive/panel_mismatch_prefit_2026-07-15/`). Speed/compute comparisons
were also affected (hapFIRE was solving against the wrong reference). See "The
fix" below for the corrected design and current results.

## What went wrong

The poolsize×depth accuracy/speed comparison fed the two tools **different
founder panels**, while scoring both against a truth derived from a **third**
source:

| | Panel actually used | Source |
|---|---|---|
| **Simulated reads + truth** | **arch3** (231 founders) | full genome assemblies, `benchmarks/p231/fastas_231/` — contains real SNPs **and indels/SVs** |
| **kMate** | **arch3** | `data/kmer_pa_231_arch3_filt2inv/` + `panel/arch3/chr1/var_pa_231_arch3_chr1.*` |
| **hapFIRE** | **greneNet_final_v1.1** (231 founders) | `benchmarks/speed_vs_hapfire/work/greneNet_chr1.vcf.gz` — old **short-read, SNP-only** VCF, a *different variant-calling lineage* |

kMate was tested on the same panel the reads came from (arch3). hapFIRE was not
— it was handed the old greneNet SNP VCF, whose genotypes do **not** correspond
to the arch3 assemblies the reads were simulated from.

## Why this made hapFIRE look bad (and why it's not hapFIRE's fault)

- The reads carry arch3 indel/SV content that the SNP-only greneNet reference
  cannot represent. **2.46 % of aligned reads (49,304 / 2,000,736) have edit
  distance NM>10** on 150 bp reads — a large population of "unexplained" reads
  HARP's SNP likelihood model cannot place.
- hapFIRE's accession-frequency solve averages over a handful of arm-scale LD
  blocks; the anomalous reads corrupt whichever blocks survive its top-20 %
  filter, so the ecotype (h) estimate collapses (R² −0.02 to −0.20 on Chr1,
  −31 to −210 genome-wide at N=150).
- The prior production/paper benchmark did **not** hit this: there the pool-seq
  reads and the reference VCF came from the *same* short-read protocol, so they
  were self-consistent. Ours were not.

**This is a benchmark-construction error on our side, not a hapFIRE
deficiency.** Any statement of the form "hapFIRE underperforms kMate on
accession recovery" from the pre-fix runs is invalid.

## What is NOT affected

- kMate's own accuracy numbers (kMate was self-consistent: arch3 reads vs arch3
  panel). The kMate-only `poolsize_depth_*` and `kmate_snp_vs_nonsnp_*` figures
  stand.
- The *shape* of the speed/CPU comparison is probably directionally right
  (kMate faster), but the absolute hapFIRE timings were spent solving against a
  mismatched reference and should be re-measured on the corrected panel.

## Correct fix

Run hapFIRE on a **SNP VCF derived from arch3** (the same panel kMate and the
reads use), not the greneNet VCF. Requirements/open questions:

1. **Need an arch3-derived SNP VCF.** (Searching whether one exists; if not,
   export one from `var_pa_231_arch3_chr1` restricted to SNP records
   `ref_len==1 & alt_len==1`.)
2. **Ploidy/phasing.** hapFIRE/HARP expects a diploid *phased* VCF. The arch3
   founders are inbred/haploid (one sequence per founder). An inbred line maps
   trivially to a homozygous diploid genotype (`0`→`0|0`, `1`→`1|1`); phasing is
   trivial because homozygous sites are unambiguously "phased". Need to confirm
   HARP accepts this encoding and that hemizygous/haploid founders are handled
   the same way the production greneNet VCF encoded them.
3. hapFIRE remains SNP-only by construction — the AF comparison can only be
   made on SNP records; kMate's indel/SV accuracy has no hapFIRE comparator
   (already the case, unchanged).

## Superseded outputs (archived under `archive/panel_mismatch_prefit_2026-07-15/`)

- `scripts/score_hapfire_vs_kmate.py` (old version — hapFIRE-on-greneNet-VCF /
  reads-from-arch3 mismatch; scored a shared AF truth that assumed shared reads)
- `results/poolsize_depth/` (per-condition hapFIRE outputs against the
  mismatched greneNet VCF, arch3 reads)
- `results/hapfire_vs_kmate_table.tsv`
- `results/hapfire_vs_kmate_{h,af}_accuracy.png`
- `results/hapfire_vs_kmate_speed.png`, `hapfire_vs_kmate_speed_vs_coverage.png`

`results/genomewide_test/`, `results/oldpanel_test/` (diagnostic runs from the
investigation) were left in place for provenance, not for reporting.

## The fix

Each tool now runs on its own native panel, with reads simulated to match, so
neither tool sees a reference it wasn't built for:

| | Panel + reads | Source |
|---|---|---|
| **kMate** | arch3 panel, arch3-simulated reads (unchanged) | `benchmarks/p231/` (no rerun) |
| **hapFIRE** | greneNet_final_v1.1 panel, reads simulated from **greneNet-derived founder FASTAs** | `scripts/build_greneNet_fastas.sh`, `scripts/run_sim_greneNet.sh`, `scripts/run_hapfire_greneNet.sh`, `sims_greneNet/`, `results/greneNet_fair/` |

Fairness of the comparison rests on one thing: both tools draw the **same
pool** per condition. `run_sim_greneNet.sh` reuses the exact arch3
founders-meta (`data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1.meta.npz`) and the
same `--seed` as the matching arch3 (kMate) sim, so `pool_weights.tsv` is
identical between the two sims for a given (N, cov, seed) — confirmed directly
on the pilot condition (N=50, cov10, s42: identical founder set).

Consequence: kMate's per-variant AF truth (from arch3's own variant calls) and
hapFIRE's greneNet SNP calls are two independent calling lineages at physically
the same positions, so there is no single unambiguous "true AF" shared across
tools anymore (this *did* exist under the old, buggy design only because both
tools shared the exact same reads/truth-generation — which is what made it
unfair). The current comparison therefore covers h-accuracy (founder-mixture
recovery, the metric both truths agree on) and speed/compute, not per-variant
AF. See `scripts/score_hapfire_vs_kmate.py` docstring.

Full N × depth × seed grid (60 conditions, matching kMate's existing arch3
poolsize×depth sweep) submitted and scored 2026-07-15; see `README.md`
"Poolsize × depth (fair...)" section for how to reproduce and current plots.

# realdata_outcross/painting/

**Question:** *where* along the chromosome does recombination happen in each real
individual?

Chromosome-painting pilot. kMate's window mode solves a local `h` per window; for a
**single individual** that local `h` is local ancestry — which founder(s) dominate at
each position, and at what dosage. Switches along the chromosome are recombination
breakpoints.

This is the one place in `benchmarks/` where window mode is used for what it is
uniquely good at. Elsewhere window mode is an accuracy liability on selfing pools
(`../../unit_ldblock_window/README.md`); here the per-window resolution *is* the
signal, not a cost.

Expected read-out, and the built-in sanity check:

| individual | expected painting |
|---|---|
| inbred | one founder, flat across the chromosome |
| F1 | two founders, few or no switches |
| recombinant | switches at breakpoints |
| complex | multiple founders, many switches |

The sample tags in `results/painted/` and `plots/` encode the expectation
(`genome_inbred_*`, `genome_F1_*`, `genome_recomb_*`, `genome_complex_*`), so a flat
painting on a `recomb` sample — or a switchy one on an `inbred` — is a red flag about
the method, not a finding.

## Two window definitions, run side by side

Both on Chr1, per sample, each writing `<SAMPLE>_Chr1_<tag>.h_blocks_per_chrom.npz`:

| tag | definition |
|---|---|
| `_win10k` | fixed 10 kb windows + HMM smoothing (the legacy "star2" recipe) |
| `_dynld` | `dynld_K500` LD blocks (`analysis/grenenet_selection/blocks/results/blocks_mcf90/`) |

Running both is the point — fixed-bp windows and LD-coherent blocks disagree about
where a breakpoint is, and on real data there is no truth to arbitrate, so the two are
shown together rather than one being picked.

kMate-only; panel is production **arch3**
(`data/kmer_pa_231_arch3_filt2inv` + `panel/arch3/…/var_pa_231_arch3`).
No hapFIRE arm — HARP has no per-window local-ancestry output to compare against.

## Layout

    scripts/   run_painting.sbatch (3-sample pilot), run_painting_full.sbatch,
               paint_genome.py (whole-genome view), paint_plot.py (per-chrom view)
    data/      samples.tsv, extremes_samples.tsv (the sample manifests)
    results/painted/   per-sample per-window h (72 npz)
    plots/     genome_<class>_<sample>.png, painting_{win10k,dynld}_Chr1.png
    logs/

Run: `sbatch --array=1-3 benchmarks/realdata_outcross/painting/scripts/run_painting.sbatch`

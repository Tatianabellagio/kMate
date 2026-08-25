# archive_gea — retired GEA tracks

**Retired 2026-08-25.** Everything in here is superseded. Do not build on it, do not
cite its numbers, and do not "fix" it. It is kept because the results are real work and
the reasoning behind retiring them is worth being able to reconstruct — not because
anything here is still in use.

## What replaced it

The live GEA track is **raw per-record LFMM, no block aggregation**:

| | |
|---|---|
| model | LFMM ridge K=16, `run_lfmm_nogif.R` (raw p, GIF logged not applied) |
| unit | clq0.9 **tiling** blocks; a block is a hit if ≥1 of its records clears the bar |
| threshold | **Bonferroni** (BH on raw p calls 60–85% of the genome — unusable) |
| inputs | `r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile/` (88 files = 4 classes × 22 axes; **lfmm only** now) |
| hit calling | `r2_gea_nonsnp/phase1_replication/multiaxis/raw_block_significance.py` |
| figures | `notebooks/raw_manhattan_*`, `gif_manhattan_*`, `newpeak_dotgrid_*` |

⚠ **`wza_in_clq09_tile/` is NOT WZA output.** Despite the name it is the *raw
per-record* input table, and it is the live track's primary data. It stays put. Do not
move it here on a name match — that mistake would break the entire live analysis.

## Why these were retired

**WZA block aggregation** (Booker 2024) is ~170× anti-conservative on this data. The
SNP-number correction standardises the block-Z's first two moments correctly (null mean
−0.001, sd 1.008) but then reads p off a **Normal**, while within-block LD makes the
true null heavy-tailed (excess kurtosis 0.79–0.86; SV 0.32). Standardising cannot change
tail shape, and Bonferroni at ~57k blocks probes ~4.8 SD — exactly where the fat tail
does maximum damage. Under a signal-free permutation null the pipeline returns ~8.5
Bonferroni blocks per scan where 0.05 is expected. Booker et al. document this deviation
themselves and only ever claim the parametric p is *quasi*-uniform.

Correcting it honestly (`empirical_null_hits.py`, in here) collapses ~7,200 nominal
block-hits to **6**, all SNP-class, all in two adjacent Chr4 ~11.9–12.2 Mb loci, with
zero nonsnp/smallindel/sv surviving. That is a defensible result but not a workable
candidate list, which is why the project moved to raw + acknowledged inflation instead.

**quasi-binomial and kendall** were the other two per-record models. They agreed with
LFMM on the axis-level picture and added a 3× multiple-testing and compute burden
without changing any conclusion. One live model is a clearer signal to future work.

## Layout

Mirrors the original tree exactly, so a path in an old note or commit message still
resolves by prefixing `archive_gea/`. Per the tree's LAYOUT.md rule 8, **paths inside
archived scripts were deliberately NOT rewritten** — they point at where things lived
when they were written. That is the record, not a bug.

```
wza/                     the WZA method section: wza_script.py + investigation/
                         (null_tail_shape.py, perm_null_regime.py, honest_hits_check.py,
                          mean_fit_*, the calibration diagnostics)
extras/driver_passenger/ WZA-dependent side investigation
r2_gea_nonsnp/phase1_replication/
    run_wza.py run_kendall.py run_binomial.py run_quasibinom.py + their sbatch/sh
    build_significant_genes*.py compare_*.py sv_*_clq90.py STATUS_clq90.md
    multiaxis/  ma_{kendall,quasibinom,wza_*}.sbatch, empirical_null_hits.py,
                significant_blocks_genes_table.py, class_peaks_overlap_table.py,
                aggregate_multiaxis.py, compare_block_arms.py, _build_*_nb.py
    results/multiaxis/  wza/ wza_arms/ kendall/ quasibinom/ empirical_null/
                        wza_in/ wza_in_clq05/ wza_in_clq05_tile/ wza_in_clq09/
                        wza_in_clq09_tile/  <- ONLY the kendall_*/binomial_* files;
                                               the lfmm_* files stayed live
notebooks/               21 rendered notebooks (WZA regime selection, fit inspection,
                         block-level Manhattans, kendall/binomial comparisons)
MOVED.tsv                every path moved, from -> to
```

Notebooks were classified by **what they actually read in their source cells**, not by
filename — scanning the whole `.ipynb` JSON gives false positives, because base64 PNG
outputs contain arbitrary substrings (`temporal_s_nofilter.ipynb` matched "deg7" inside
an image blob and is an r1 notebook, so it stayed live).

## Result data is gitignored

~66 GB. `.gitignore` rules for these paths were added in the same change as the move —
the pre-existing rules are literal paths and silently stopped matching once the parent
moved (LAYOUT.md rule 3). If you move anything in here again, re-check with:

```bash
git check-ignore -v analysis/grenenet_selection/archive_gea/<new path>
```

## If you need a number from here

Prefer re-deriving it on the live track. If you must cite the WZA-era work, cite the
**empirical-null** figures (`results/multiaxis/empirical_null/summary_empirical.csv`),
never the nominal `Z_pVal` ones — the nominal counts are the ~170×-inflated ones, and
several committed narratives (the "159 non-SNP block-hits → 90 blocks → 102 genes",
HSP26.5, CRK13-16) rest on them and do not survive recalibration.

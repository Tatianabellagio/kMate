# local-only on p231 — benchmark

Benchmarks the global-free window mode (`--local-only`, commit dc90850) on the
**p231 production panel at 10×**, across **all 11 simulation setups, seed s42**
(g0/g1/g3 × n50/n231 × outcross/self97 × balanced/dom500/dom500nr).

`--local-only` abstains (NaN AF) wherever a window is below the k-mer floor /
empty / unassigned, instead of falling back to the chrom-wide h. So **coverage
(call-rate) and accuracy-on-called are separate axes** and are scored separately
by `scripts/score_localonly.py` (NaN-aware; reusing `build_benchmark_table.py`
would have produced all-NaN metrics).

## Run

```
sbatch --array=0-10 scripts/run_localonly_p231.sbatch
```

Per pool, count k-mers ONCE then run 5 modes and score each:
`global` ; `anchored`/`localonly` × {`dynldK500` blocks, `w10kb` fixed windows}.
`min-kmers-per-block 50` for all window runs (so anchored vs localonly differ
ONLY in fallback-vs-abstain). Truth = `recomb_truth_raw.tsv.gz` (`truth_af`).

## Result (ALL-class, call_rate / R² on called)

* **Anchored ≈ global everywhere** (anchor 0.3 + fallback recovers chrom-wide
  accuracy at 100% call): mean R² 0.977.
* **local-only costs BOTH coverage and accuracy**: mean call-rate 66% (dynldK500)
  / 88% (w10kb); R² on called ~0.93 → a **−0.04 to −0.05 R² penalty even where it
  speaks**, on top of the abstained records.
* **10kb windows beat dynld_K500 for local-only on both axes** (88%/0.936 vs
  66%/0.928) — opposite of the coarser=more-k-mers intuition.
* **No selfing advantage**: selfing-balanced setups (where global is near-perfect,
  R²≈0.994) show the LARGEST relative local-only penalty. local-only closes the
  gap to anchored only under g3 recombination+selection (dom), where global
  itself degrades.

Confirms the standing decision ([[gea-blockwindow-benchmark-result]],
[[block-mode-kmer-floor]]): global/anchored is the right AF estimator; pure-local
is a special-purpose mode that trades coverage+accuracy for global-free purity.

Master table: `results/localonly_p231_table.tsv` (220 rows = 11×5×4).

# benchmarks/founder_h_imbalance/

**Question:** does the private-k-mer filter (`filt2inv`) still matter under the
*current* production estimator (`--unit chrom` + per-founder M-step normalization)?

Background: the ac=1 private-k-mer drop in `filt2`/`filt2inv` was introduced as a
**guard against the uneven panel** — the 78 cactus long-read founders carry richer
k-mer fingerprints than the 153 PanGenie-genotyped ones, and the EM over-credited
cactus founders (the "+41% cactus h-bias", `BACKGROUND.md` limitation 2). The
per-founder normalization (2026-07-06) removes that imbalance *at its source*, so the
filter may now be redundant. This dir tests exactly that.

## What it produces

`scripts/plot_h_imbalance.py` renders two views of the same p231 data, both
long-read-vs-short-read coloured, 2 rows (`raw` / `filt2inv`) × 3 regime columns:

| output | view |
|---|---|
| `results/h_imbalance_grid_p231_chrom.png` | per-founder `h` estimate vs truth, scatter |
| `results/h_imbalance_error_boxplot_p231_chrom.png` | distribution of per-founder signed error (est − truth), split by founder class |

The signed-error boxplot is the decisive panel: a cactus-vs-PG offset that is present
in `raw` and absent in `filt2inv` means the filter is still doing work; the two arms
overlapping means the normalization has already absorbed it.

kMate-only; no competitor arm.

Related: `../p80/` is the complementary control (a *homogeneous* 80-founder cactus
panel, where the asymmetry should be absent by construction), and
`../../docs/FOUNDER_NORMALIZATION_FIX.md` documents the normalization itself.

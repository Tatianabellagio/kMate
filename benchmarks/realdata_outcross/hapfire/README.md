# realdata_outcross/hapfire/

**Question:** on *real* outcrossed individuals, do kMate and hapFIRE agree about
which founders are present and in what proportion?

This is a **concordance** measurement, not an accuracy one. Ruth's individuals have no
ground truth — there is no `pool_weights.tsv` to score against — so the two tools are
compared to *each other*, and no row here can be read as "tool X is more accurate."

## Panel basis — each tool on its own panel (correct)

| tool | panel |
|---|---|
| kMate | **arch3** — `data/kmer_pa_231_arch3_filt2inv` + `panel/arch3/…/var_pa_231_arch3` |
| hapFIRE | **greneNet** — `benchmarks/vs_hapfire/work/greneNet_chr1.vcf.gz` |

That is the native-panel setup this repo requires (see `../../vs_hapfire/README.md`).
Both vectors are distributions over the **same 231 founder IDs**, which is what makes
them comparable even though the underlying variant sets differ.

> ⚠️ `scripts/run_hapfire_ruth.sbatch` calls this a "shared panel" in its header
> comment. That wording is wrong and predates the 2026-08-25 shared-panel decision —
> the panels are *not* shared, and it would be a defect if they were. What is shared is
> the **founder ID space**. Corrected inline; flagged here so the old phrasing doesn't
> get taken at face value.

Chr1 only — the hapFIRE-ready panel VCF is Chr1, so kMate's Chr1 `h` is used.

## Metrics — distribution distances, deliberately not r/R²

`scripts/compare_ecotype_vs_h.py` compares with:

| metric | meaning |
|---|---|
| effective #founders | `1/Σf²` from each method |
| **TVD** | `0.5·Σ\|h_kMate − h_hapFIRE\|` — 0 identical, 1 disjoint; literally the fraction of founder mass the two disagree about |
| **Hellinger** | `sqrt(0.5·Σ(√p−√q)²)`, in [0,1] |
| top-founder agreement | do both put the same founder first? |

Pearson r / R² are **avoided on purpose**: these vectors are mostly zeros with a few
large values, so correlation is dominated by the shared tail and is meaninglessly
inflated. Don't add them back.

## Layout

    scripts/   run_hapfire_ruth.sbatch (per-sample hapFIRE), compare_ecotype_vs_h.py
    data/      samples.tsv
    results/   compare_ecotype_vs_h.csv + per-sample hapFIRE _ecotype_frequency.txt
    plots/     compare_ecotype_vs_h_{scatter,effn}.png
    logs/

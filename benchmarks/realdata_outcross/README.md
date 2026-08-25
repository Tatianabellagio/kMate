# `benchmarks/realdata_outcross/` — outcrossing / MOI detection (kMate vs hapFIRE)

Can kMate tell an **inbred** individual from an **outcrossed** one, the way Ruth's
hapFIRE workflow does?

## The question

Ruth Epstein grows plants in the greenhouse from GrENE-net seeds, **individually**
(not pooled) deep-sequences them, and runs each through **hapFIRE** to ask whether the
individual matches **one** of the 231 founders (→ inbred / selfed lineage) or **several**
(→ outcrossed hybrid). This is a multiplicity-of-infection (MOI) / outcrossing call.

kMate's EM already solves exactly this object: the **founder-mixture vector `h`** on the
231-founder simplex. For a single individual,

- **inbred** → `h` collapses onto ~1 founder (`max_h ≈ 1`, `eff_n_founders = 1/Σh² ≈ 1`),
- **outcrossed** (e.g. F1) → `h` spreads across ≥2 founders (`max_h ≈ 0.5`, `eff_n ≈ 2`).

So no new estimator is needed — we run kMate global mode and read out `h`.

## Validation set (hapFIRE-labelled, from Ruth)

`manifest.tsv` — 5 confirmed **outcrossed** + 5 confirmed **inbred** 231-parent
individuals. Input = the **deduped BAMs** under
`…/adaptest_N800_deepseq/dedup/` (these are PCR-based deep-seq libraries, so the
already-deduped reads are the right input). kMate counts k-mers straight from the BAM
via `samtools fastq`.

| label | samples |
|---|---|
| outcrossed | 754, 320, 711, 308, 62 |
| inbred | 101, 113, 93, 117, 202 |

## Panel

Production 231-founder arch3 panel — same as the GrENE-net cohort:
`data/kmer_pa_231_arch3_filt2inv` + `panel/arch3/chr{N}/var_pa_231_arch3_*`.

## Pipeline

```
run_moi.sbatch     array 1-10: per sample, build a node-local k-mer DB once from the
                   BAM, then global EM for all 5 chroms -> out/<sample>_<chrom>.tsv
                   + .h_per_chrom.npz (the h vector we read out)
analyze_moi.py     average per-chrom h -> genome-wide founder composition; compute
                   eff_n_founders / max_h / top founders; compare to hapFIRE labels;
                   write moi_summary.csv + moi_summary.png
```

Run:

```bash
sbatch --array=1-10 benchmarks/realdata_outcross/scripts/run_moi.sbatch
# when done:
python benchmarks/realdata_outcross/scripts/analyze_moi.py
```

## Output

- `moi_summary.csv` — per-sample `eff_n_founders`, `max_h`, top-3 founders, hapFIRE label.
- `moi_summary.png` — eff_n_founders by label (left) + genome-wide founder composition (right).
- Console prints whether the two classes separate cleanly and the threshold.

## Sub-scopes

Two further analyses on the same individuals, each self-contained with its own README:

| dir | question |
|---|---|
| [`hapfire/`](hapfire/README.md) | do kMate and hapFIRE **agree** about founder composition? A concordance measurement (TVD / Hellinger) — real data has no ground truth, so no accuracy claim is possible. Each tool on its **native** panel: kMate arch3, hapFIRE greneNet. |
| [`painting/`](painting/README.md) | **where** does recombination happen along the chromosome? kMate window mode as local ancestry; two window definitions (fixed 10 kb + HMM, and dynld_K500) run side by side. |

Layout note: this dir's own outputs moved 2026-08-25 — `out/` → `results/h_per_chrom/`,
figures → `results/plots/`, scripts → `scripts/`, manifests → `data/`. Both sub-scopes
follow the same shape.

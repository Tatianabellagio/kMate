# Phase-1 GEA replication on kMate AF — SNPs, small indels & SVs

**Goal.** Re-run the GrENE-Net **phase-1** genotype–environment-association (GEA)
pipeline — **Kendall-τ + LFMM (K=16) + binomial regression → WZA block
aggregation** — but on **kMate** allele frequencies, which (unlike the SNP-only
phase-1 paper) cover **SNPs, small indels, *and* SVs**. Question: do we recover the
same significant blocks (headline: **CAM5 / AT2G27030, Chr2 ~11.53 Mb**), and do
indels/SVs reveal blocks the SNP-only scan misses?

Owner: T. Bellagio. Started 2026-06-15.

---

## Decisions (locked 2026-06-15)

| Decision | Choice | Rationale |
|---|---|---|
| Variant classes | **Separate WZA per class** (SNP / smallindel / SV run independently, block lists compared) | Cleanest answer to "do indels/SVs add blocks SNPs miss"; isolates each class's contribution to CAM5. |
| Climate axis | **bio1 only** (annual mean temperature) | The CAM5-relevant axis; keeps the first pass focused. |
| Generations | gen1 · gen3 · **gen9 = last-gen** (flower-weighted timepoint merge) | gen9 is the phase-1 analysis unit; lfmm/binomial run on gen9, kendall on all three. |
| Analysis unit | **site_gen_plot pool** (flower-weighted), 745 pools = phase-1 columns exactly | `build_pool_matrix.py`; verified vs phase-1 merged_hapFIRE to ~1e-15. |
| Low-freq filter | **MAF ≥ 0.05** + min-pool-coverage (finite AF in ≥ N pools) | Phase-1 used `maf05 mincount05`; low-freq alleles flagged by user. |
| LD blocks | phase-1 hapFIRE SNP blocks (16,674; `lib.assign_ld_blocks`) | indels/SVs inherit the nearest genotyped-SNP block. |

Variant class defs (`lib.py`): SNP = `ref_len==1 & alt_len==1`; SV =
`|alt_len-ref_len| > 50`; smallindel = non-SNP with `|alt_len-ref_len| ≤ 50`.

---

## WZA code review (vs Booker et al. 2024, Mol. Ecol. Resour.)

Reviewed `../wza_script.py` (our in-repo copy) against the **canonical
`general_WZA_script.py` from github.com/TBooker/WZA (master)** and the phase-1
GrENE-Net variant `general_WZA_script_mod_polynomial_order7.py`.

**Verdict: our copy is the canonical Booker WZA and the core math is correct.**

- Core weighted-Z is identical to canonical:
  `z = Φ⁻¹(1−p)`, weight `pq = MAF(1−MAF)`,
  `Z_w = Σ pq·z / sqrt(Σ pq²)`. ✓
- SNP-number correction identical to canonical: rolling (roller=50,
  minEntries=40) mean/SD of Z vs SNP-count, **deg-2** polynomial fit,
  empirical `Z_pVal = 1 − Φ(Z; loc=mean̂, scale=sd̂)`. ✓

Two intentional deviations from canonical, both fine:
1. **`--maf_filter` exposed as a CLI arg** (canonical hardcodes 0.05 in `main`).
   Our default is `0.0` → **must pass `--maf_filter 0.05`** to match canonical,
   *or* (as we do) apply the MAF≥0.05 cut upstream in the class-matrix builder.
2. **SD safety-floor** clip on the polynomial's predicted SD (canonical has none →
   `norm.cdf` returns NaN in sparse large-window tails). Strict improvement; only
   touches windows that would otherwise be NaN.

**Which correction we use (UPDATED 2026-06-17, supersedes the original deg-2 plan).**
Phase-1 used a **deg-7** polynomial (`chosen_degree = 7`, `minEntries = 10`), uncapped.
The `../wza_investigation/` study (RESULTS.md, `wza_investigation.ipynb`) showed **deg-7
genuinely fits the SNP-count→SD support curve better** (RMSE 0.94 vs deg-2's 2.06) — a
quadratic is too rigid for our 1→9,158-SNP LD blocks. deg-7's only danger is *tail
extrapolation*, which a SNP cap removes. So the canonical/primary correction is now
**deg7-cap2000** (`--poly-deg 7 --min-entries 10 --sample-snps 2000`), with **deg7-nocap**
(= phase-1 exactly) and **deg2nocap** (canonical Booker) as sensitivities. The cap is
load-bearing only under deg-2; under deg-7 cap vs nocap barely moves CAM5. The earlier
"deg-2 primary, deg-7 sensitivity" plan was too strong and is retired.

---

## Status

See **`STATUS.md`** for the live state (branch matrix, coverage, CAM5 result, NEXT).
As of 2026-06-17: all 3 models × 3 classes built for gen9; kendall also gen1/gen3;
WZA at deg7-cap2000 (primary) + deg7-nocap + deg2nocap; CAM5 reproduces in all 3 models.

## Layout (this folder)
- `README.md` — this file (plan, decisions, WZA review, status).
- *(scripts added as built)* — `build_class_matrices.py`, `run_kendall.py`,
  `run_wza.py`, `run_lfmm.*`, `run_binomial.py`, `compare_blocks.py`.
- Outputs → `analysis/grenenet_gea/r2_gea_nonsnp/phase1_replication/results/`.

## Provenance
- Original phase-1 pipeline: `/global/scratch/users/tbellg/gea_grene-net/`
  (`kendall_tau/`, `lfmm/`, `binomial_regression/`, `wza/`, `cam5/`,
  `signficant_intersection_GEA_models/`).
- Canonical WZA: github.com/TBooker/WZA · Booker et al. 2024, *Mol. Ecol.
  Resour.* (`papers/Booker et al. 2024 - Mol. Ecol. Resour_.pdf`).
</content>
</invoke>

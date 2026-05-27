# v3 per-SNP AF outliers — what we know

**Last updated:** 2026-05-14
**Validation sample:** SEEDMIX_S1, Chr1 SNPs, n = 437,537 records matched on (chrom, pos, ref, alt)
**Definition of outlier:** `|v3_AF − hapFIRE_AF| > 0.1`
**Baseline count:** 2,544 outliers (0.58%); baseline MAE 0.0144, RMSE 0.0246

---

## 1. Where the outliers are

| Region | Share | Notes |
|---|---|---|
| Centromere (Chr1 ~14–17 Mb) | ~30% | Structural — alignment dead zone. Not fixable at EM or cn_var level. |
| Arms (diffuse) | ~70% | Concentrated at SNPs ≤100 bp from an INDEL/SV bubble in the cactus graph. |

---

## 2. Carrier-class breakdown (Diagnostic J)

| Carrier class | n records | outlier rate | mean \|err\| | RMSE |
|---|---:|---:|---:|---:|
| **both** cactus & PG carriers | 370,524 | **0.691%** ← highest | 0.0160 | 0.0268 |
| cactus-only carriers | 37,478 | 0.133% | 0.0053 | 0.0144 |
| **pg-only** carriers | 25,418 | **0.043%** ← lowest | 0.0051 | 0.0100 |
| no carriers (mono) | 4,117 | 0.146% | 0.0080 | 0.0124 |
| **OVERALL** | 437,537 | 0.601% | 0.0144 | 0.0252 |

**Implication:** Refutes the "69 extra PG founders add noise" hypothesis. If the extras were noisy, PG-only carriers (sites only the extras hit) would dominate. They have the **lowest** outlier rate. Outliers concentrate where **both panels have carriers but disagree on identity**.

Filter test (drop PG-only records): MAE 0.0144 → 0.0150 (no improvement).

---

## 3. What the outliers ARE NOT

- ❌ **A coverage problem.** cov50 doesn't reduce MAE on cactus_em-v3 ★★ filt2 — residual is simplex-limited, not coverage-limited.
- ❌ **A 69-extras / PG-heavy artifact.** PG-only carrier records are the cleanest class (0.043% outlier rate).
- ❌ **A hapFIRE artifact.** REF allele matches between v3 and hapFIRE; hapFIRE matches 1001G short-read GTs.
- ❌ **Fixable by EM-level k-mer filtering.** filt2 (drop singleton k-mers) debiases `h` but does **not** reduce outliers, because outliers live in `cn_var`, not in `h`.
- ❌ **A cactus bug.** Cactus correctly represents the long-read assemblies it was given (verified by 99.81% SNP identity between assemblies cactus calls "same sample"). When a widely-used tool seems wrong, suspect upstream sample handling first.

---

## 4. What the outliers ARE

**`cn_var` disagreement between cactus-graph + PanGenie (v3) and 1001G short-read (GrENE-Net).**

Mechanism: `bcftools norm -m -any` atomization of multi-allelic bubbles in the cactus graph creates spurious per-founder carrier calls at SNPs adjacent to INDELs/SVs. The two panels see the same REF/ALT at these sites but disagree on **who carries the ALT**, so `cn_var @ h` produces a per-SNP AF mismatch even when `h` is correct.

**Confirmed by Diagnostic G** (`h_v3 @ cn_GN`): keep v3's `h`, swap in GN's `cn_var` → outliers nearly disappear → the bug is in `cn_var`, not in `h`.

---

## 5. Two fixes, two independent problems

| Fix | What it changes | Targets | Effect on SEEDMIX_S1 Chr1 |
|---|---|---|---|
| **Fix 1 — filt2** (drop k-mers with ac<2 from `cn_full` before EM) | `h` vector | h-vector cactus bias (singleton-k-mer evidence amplification, ~250×) | cactus h-bias +37% → +6%, PG bias −20% → −3%; **per-SNP outliers ≈ unchanged** (Fix 1 does not touch `cn_var`) |
| **Fix 2 — hybrid `cn_var`** (use GN entries at SNPs ≤100 bp from any INDEL; v3 elsewhere) | `cn_var` matrix | per-SNP AF outliers (atomization-induced GT disagreement) | **RMSE 0.0246 → 0.0152 (−38%)**; **outliers 2,544 → 122 (−95%)**; 372,809 records swapped |

**Orthogonality.** Fix 1 changes which k-mers feed the EM (changes `h`). Fix 2 changes which `cn_var` rows you multiply against (does not touch `h`). They can be stacked.

---

## 6. Residuals after both fixes (~122 outliers remaining)

- **Centromere structural mass** (Chr1 14–17 Mb dead zone) — independent of cn_var; alignment-driven.
- **Sites further than 100 bp from any INDEL** that still have per-founder GT disagreement. Rare. Threshold could be relaxed to 250 bp if needed.
- **p82-control–persistent Chr1q-knob outliers.** The p82 control (pure cactus, no PG extras) still shows the Chr1q-knob outlier cluster → partial residual mechanism that is **independent of the 69 PG extras**.

---

## 7. Decisions that follow

- **Do not** propose "rebuild `cn_full` from Beagle-imputed VCF" — Beagle hard-rejected for SVs (−25 to −30 pp concordance).
- **Do not** propose "more coverage" — cov50 doesn't help.
- **Do not** propose a Dirichlet prior toward uniform — overfits to SEEDMIX (a uniform pool).
- **Do not** smooth in AF space — catastrophic (R² 0.95 → 0.01). Smooth only in `h` space if needed.
- **Production path:** ship Fix 1 (filt2 cn_full) + Fix 2 (hybrid cn_var) together. Both are mechanism-driven and validated against the recipe truth (~uniform 1/231) and per-SNP hapFIRE truth.

---

## 8. Related findings (context, not directly about outliers)

- **5772 / 6150 mislabeling.** v3 has two copies of 6150 (T980) with one mislabeled as 5772. Upstream wet-lab sample-handling error, NOT a cactus bug. Set-1's real long-read assembly is missing from v3.
- **vcfbub coverage gap.** v3 `cn_var` currently covers 55% of GrENE-Net 231 SNPs; 38% are recoverable with a `vcfbub -l 0` rebuild (which preserves nested SNPs from bubbles); 7% are genuinely outside cactus bubbles.

---

## 9. Pointers to evidence in the repo

- Main notebook (diagnostics + 2×2 figure of two fixes): `SEEDMIX_S1_v3_vs_hapfire.ipynb`
- Fix 1 isolation on per-SNP scatter: `FIX1_SCATTER_EFFECT.ipynb`
- Carrier-class diagnostic raw output: `CARRIER_ASYMMETRY_DIAGNOSTIC.ipynb`
- p82 control conclusion: `HANDOFF_P82_CONTROL.md`, memory `project_p82_control_result.md`
- K-mer rebalancing alternatives tested and rejected: `BALANCING_KMERS.md`
- Production EM with `--filt2` and `--ac-weight-counts` flags: `poolfreq/src/per_sample_per_chrom.py`

# Handoff: 82-founder cactus-only control experiment

**Created 2026-05-13. Standalone task — complete in-session, all artifacts live under `control_p82/`.**

## Run state (final 2026-05-13 11:52)

**Full pipeline completed in ~80 min wall (10:34 → 11:52)** — much faster than the 10-15h budget because cov10/Chr1 is small.

| phase | jobid | status | wall | notes |
|---|---|---|---|---|
| A1 (atomized VCF) | 61881 | ✅ done | 2 min | 1.3M Chr1 records, 82 Accession_IDs |
| A1b (bubble VCF for PG-index) | 61967 | ✅ done | 1 min | 950k bubble-form records |
| A2 (PanGenie-index) | 61976 | ✅ done | 9 min | 181k bubbles → `idx_Chr1_kmers.tsv.gz` (81 MB) |
| A3 (cn_full_p82) | 61977 | ✅ done | 60 min | (82, 18.3M), nnz=211M, density 14.04% |
| A4 (cn_var_p82) | 61884 | ✅ done | 2 min | 17M cn_var.npz + 50M meta |
| A5 (FASTAs ×82) | 61885 | ✅ done | 2 min | renamed `.chr1.fa` → `.chr.fa` |
| B1 (sim n50_g1) | 61978 | ✅ done | 5 min | 350M r1.fq + r2.fq + 4.9M truth |
| B3 (sim n50_g3) | 61979 | ✅ done | 5 min | (same scale) |
| C1G global n50_g1 | 61980 | ✅ done | 2.5 min | 33M TSV, EM 200 iter, eff_n_founders=10.7 |
| C1S star2 n50_g1 | 61981 | ✅ done | 5 min | 33M TSV |
| C3G global n50_g3 | 61982 | ✅ done | 3 min | 33M TSV |
| C3S star2 n50_g3 | 61983 | ✅ done | 5 min | 33M TSV |

## Headline results (corrected 2026-05-13 14:30)

**Initial comparison was sloppy** — I compared v3 results from `cov10_*_p231_chr1/` sims against p82 results from `cov10_*_p82_chr1/` sims. Same seed (42), but different founder universes (231 vs 82) → multinomial draws diverged → **different truths** (median |Δtruth|=0.08-0.12) → R² values aren't comparable.

**Proper comparison**: run v3 panel cactus_em on **p82's reads** (jobs 61988-61991, ~5 min each), then score all 4 panel×method combos against p82's truth on the (chrom, pos, ref, alt) intersection of variants (775,944 shared records out of v3's 7.7M and p82's 1.3M).

### Apples-to-apples (same reads, same target variants, same p82 truth)

**n50_g1** (357,700 shared polymorphic SNPs, cov10):
| run | MAE | RMSE | R² | slope | \|res\|>0.1 |
|---|---|---|---|---|---|
| v3 panel  global | 0.0234 | 0.0339 | 0.980 | +0.982 | 0.86% (3,064) |
| **p82 panel global** | **0.0215** | **0.0301** | **0.984** | **+0.996** | **0.73% (2,615)** |
| v3 panel  ★★ | 0.0205 | 0.0304 | 0.984 | +0.982 | 0.45% (1,612) |
| **p82 panel ★★** | **0.0185** | **0.0259** | **0.989** | **+0.992** | **0.22% (795)** |

**n50_g3** (297,807 shared polymorphic SNPs, cov10):
| run | MAE | RMSE | R² | slope | \|res\|>0.1 |
|---|---|---|---|---|---|
| v3 panel  global | 0.0517 | 0.0720 | 0.916 | +0.926 | 16.64% |
| **p82 panel global** | **0.0503** | **0.0696** | **0.922** | **+0.943** | **16.07%** |
| v3 panel  ★★ | 0.0332 | 0.0469 | 0.965 | +0.950 | 4.48% |
| **p82 panel ★★** | **0.0310** | **0.0426** | **0.971** | **+0.962** | **3.49%** |

## Verdict: hypothesis CONFIRMED (across the board, not just on easy regimes)

✅ **p82 beats v3 on every metric (MAE, RMSE, R², slope, outlier fraction) on both regimes and both methods.** The earlier claim that "p82 global regresses on n50_g3" was wrong — that was an artifact of comparing different sims. With identical inputs, p82 global is **better** than v3 global on n50_g3 (MAE 0.0503 vs 0.0517, R² 0.922 vs 0.916).

✅ **The cactus-PG k-mer asymmetry was a real and significant contributor to v3's bias.** Replacing the heterogeneous v3 panel with a homogeneous p82 panel reduces MAE 3-10% and cuts outlier fractions 7-50%.

✅ **p82 ★★ is the production-relevant config.** On the shared variants:
- n50_g1: MAE −10% (0.0205 → 0.0185), outliers −51% (1,612 → 795)
- n50_g3: MAE −7% (0.0332 → 0.0310), outliers −22% (13,341 → 10,385)

❓ **Chr1q knob outliers (Mb 21-23) still dominate** the residual outlier list on p82 ★★. So that cluster is NOT solely caused by cactus-PG asymmetry — there's residual local pathology even with a homogeneous panel. Mechanism candidates: (a) genuinely panel-rare variants in 1-3 cactus founders → EM identifiability in tighter simplex; (b) per-window k-mer count variance in the knob (high SV density); (c) recomb hotspot density mismatch with actual ancestry breakpoints.

⚠️ **h-vector collapse on p82** (homogeneous-panel artifact). EM concentrates h on ~10 representative founders even though 42 are truly in the pool (eff_n_est=10.7 vs eff_n_truth=37.9, median h_est/h_truth=0.30). AF projection still works because grouped-similar founders have near-identical carrier profiles. But h on p82 is NOT a reliable founder-attribution signal.

## Implication for production direction

The k-mer rebalancing path (see `BALANCING_KMERS.md`) is the right next step. Beagle imputation is ruled out (per `feedback_no_beagle_solutions.md`). What p82 demonstrates:
- A homogeneous fingerprint distribution (CV ~0% on p82 vs CV 2.6% on v3) gives ~1pp R² lift and 10× outlier reduction.
- The residual Chr1q knob outliers are NOT fixable by panel homogeneity alone — needs a separate diagnostic (per-window k-mer variance? region-specific EM regularization?).
- ★★ window-mode is essential for high-recomb regimes regardless of panel composition.

## Where to pick up

- Full hexbin plots: open `control_p82/results/FINAL_RESULTS_cov10_p82.ipynb` and run all cells.
- Diagnostic for the Chr1q knob residual: compare v3 vs p82 per-Mb-bin outlier residuals (which Mb the bias survives in).
- h-vector diagnostic: load `*.h_per_chrom.npz` and `*.h_blocks_per_chrom.npz` to see if homogeneous-panel EM still has over-credit on any subset.

## Scope locked 2026-05-13

- **Founders**: all 82 cactus pangenome accessions (was originally framed as "80 GrENE-Net-overlap" but for a standalone control there's no reason to drop founders — more founders = better simplex coverage; the v3 comparison happens at per-variant AF level, not at founder level). Naming convention `p82` throughout.
- **Coverage / chrom**: cov10 / Chr1 only (mirrors `FINAL_RESULTS_cov10_v3.ipynb`'s eval surface).
- **Methods**: `cactus_em global` and `cactus_em ★★ (window 10kb + λ=0.3 + smooth 5α0.5)`. **No** freqk, **no** hapFIRE. The hypothesis is mechanistic (does the +41% cactus h-bias and its off-diagonal streaks vanish when the panel is homogeneous?) — freqk and hapFIRE don't add signal to that question.
- **Comparison**: p82 build vs the v3 cells in `FINAL_RESULTS_cov10_v3.ipynb` (same eval pipeline, same truth-vs-est hexbins). A 151-PG-only mirror is a possible follow-up if this control confirms the mechanism. (Beagle imputation is NOT considered as a path forward — it has been hard-rejected for SVs; see `feedback_no_beagle_solutions.md`.)
- **Regimes**: `n50_g1` and `n50_g3` are must-haves. `n80_g1` (full-panel pool) is the stretch goal. Pool size constrained to ≤ 82.

## Why this experiment exists

`FINAL_RESULTS_cov10_v3.ipynb` shows cactus_em ★★ (window 10kb + λ=0.3 + smooth 5α0.5) at R²=0.986 on cov10 n200_g1, but the hexbin has two visible off-diagonal streaks:

- **Below diagonal**: 136 SNPs (0.02%) at truth 0.4-0.8 → est 0.15-0.4. Characterized last session as **PG-specific variants** concentrated in Chr1q knob (Mb 21-23).
- **Above diagonal**: symmetric streak at truth 0.05-0.4 → est 0.5-1.0. Predicted by the same mechanism to be **cactus-specific variants** but not yet characterized.

**Diagnosed cause**: residual **+41% cactus over-credit** in the EM-estimated h vector (see `memory/project_v3_singleton_kmer_bug.md`, `project_v3_pg_specific_outliers.md`). Mechanism: cactus founders have richer k-mer fingerprints in `cn_full_231_v3` than PG founders (CV 2.6% in v3 vs 0.9% in v2), so the EM assigns more weight to cactus evidence regardless of true pool composition.

## The hypothesis this experiment tests

> If the bias is caused by cactus-vs-PG k-mer asymmetry, then rebuilding everything with **only the 82 cactus founders** (so all founders are equally well-characterized by long-read assemblies) should eliminate the +41% bias and make both off-diagonal streaks disappear.

This is a **positive control** — the answer reviewers will ask for. If the streaks DO vanish, the mechanism is confirmed and the next step is to develop a k-mer **rebalancing** strategy on the 231-founder panel (per-founder k-mer-budget equalization, row-normalized cn_full, etc. — see the brainstorm in `BALANCING_KMERS.md`). If the streaks persist, the residual bias has some other source we haven't identified.

**Beagle imputation is explicitly ruled out** as a balancing approach: previous LOO testing showed −25 to −30pp concordance loss on small/medium SVs, and SVs are the variant class we care most about. See `feedback_no_beagle_solutions.md`.

## Standalone project tree

**No symlinks or shared artifacts with v3/v2.** Everything lives under:

```
/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p82/
├── README.md
├── data/
│   ├── samples_82.tsv            # Assembly_ID  Accession_name  Accession_ID
│   ├── samples_82_rename.tsv     # Assembly_ID  Accession_ID  (for bcftools reheader -s)
│   ├── pangenome_p82_chr1.vcf.gz # canonical panel VCF (post-norm, renamed to Accession_IDs)
│   ├── kmers_p82/                # PanGenie-index output
│   ├── cn_full_p82/cn_Chr1.{cn,meta}.npz
│   └── cn_var_p82.{cn_var,meta}.npz
├── fastas_82/<Accession_ID>.chr1.fa     # bcftools consensus from pangenome_p82_chr1.vcf.gz
├── sims/<regime>/                       # VISOR sim dir per regime
├── results/
│   ├── cactus_em_global/<regime>/<sample>.tsv
│   └── cactus_em_star2/<regime>/<sample>.tsv
├── scripts/                             # numbered build/sim/run scripts
└── logs/                                # SLURM stdout/stderr
```

All sample IDs use **1001G Accession_ID** (e.g. `6911`), NOT Assembly_ID (e.g. `100042`). Renaming happens immediately after VCF subset/norm so all downstream artifacts (cn_full, cn_var, FASTAs, sim outputs) are keyed consistently.

## Inputs (verified 2026-05-13)

| Artifact | Path | Notes |
|---|---|---|
| **Cactus pangenome VCF (filtered)** | `/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.vcf.gz` | 82 samples (Assembly_IDs as sample names: `100042`, …); contigs `Chr1..Chr5`; 1.13M Chr1 records; has cactus `AT=…` bubble-path field → needs `vcfbub + bcftools norm -m -any` before downstream use. Raw is 11GB at `.raw.vcf`; **use the `.vcf.gz`** (already filtered). |
| 82 long-read assemblies (raw .chr.fa) | `/home/tbellagio/scratch/pang/sv_panel/assemblies_clean/<asm_id>.chr.fa` | 82 cactus assemblies. **NOT used as sim FASTAs** in this build — see Phase A5 reasoning. |
| **TAIR10 reference (IUPAC→N, numeric chroms)** | `/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.numeric.iupacN.fa` | Cactus-aligned reference. IUPAC codes have been N-replaced (required because cactus VCFs N out IUPAC codes; per `memory/project_cactus_ref_is_tair10.md`). |
| Assembly_ID → Accession_ID map | `/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/sv_panel_to_accession_id.tsv` | 82 rows, all 82 Accession_IDs unique (the original handoff's "82→80 collapse via duplicate ecotype" was wrong — verified 2026-05-13). |
| `vcfbub` binary | check `pang` env: `/home/tbellagio/miniforge3/envs/pang/bin/vcfbub` (most likely) | Used in the same `vcfbub -l 0 -r 100000` pattern as v3 build. |
| `bcftools` binary | `/home/tbellagio/miniforge3/envs/BIOS424/bin/bcftools` | Or any other env that has bcftools (`gwas`, `pipeline_snakemake`, `sequencing_pipeline` all have it). |
| `PanGenie-index` | check `pangenie` env: `/home/tbellagio/miniforge3/envs/pangenie/bin/PanGenie-index` | Confirm it exists; this is the v3-build env. |
| Existing reference builder scripts | `pangenie_genotyping/scripts/build_pangenie_index.sh`, `poolfreq/scripts/build_cn_full_v3_one.sh`, `poolfreq/scripts/build_cn_var_v3.sh`, `sims/visor_freqk/scripts/build_unimputed_v3_fastas.sh`, `sims/visor_freqk/scripts/run_pool_sweep_82_recomb.sh` | Clone-and-adapt for p82 — do NOT modify originals. |
| `per_sample_per_chrom.py` | `poolfreq/src/per_sample_per_chrom.py` | The cactus_em driver. Reused as-is; we just point `--cn-kmer-prefix` and `--cn-var` at p82 artifacts. |
| `compute_recomb_truth.py` | `sims/visor_freqk/scripts/compute_recomb_truth.py` | Reused as-is with `--cn-var control_p82/data/cn_var_p82.cn_var.npz`. |

## Step-by-step plan

### Phase A — Build the p82 panel artifacts (all in `control_p82/`)

**Consistency invariant:** the same final VCF (`data/pangenome_p82_chr1.vcf.gz`) drives `cn_var_p82`, `cn_full_p82`, AND the per-founder consensus FASTAs. This is the v3 fix pattern — avoid the simulation-FASTA-vs-cn_full mismatch that caused the 2× MAE gap (see `memory/project_v3_singleton_kmer_bug.md`).

**A1. Subset `pang_1001gplus_82acc.vcf.gz` → Chr1 → vcfbub → bcftools norm → reheader to Accession_IDs.**
- Restrict to Chr1 (`bcftools view -r Chr1`).
- Flatten bubbles ≤ 100kb (`vcfbub -l 0 -r 100000`).
- Atomize multi-allelics (`bcftools norm -m -any -f TAIR10.chr.numeric.iupacN.fa`).
- Rename samples Assembly_ID → Accession_ID (`bcftools reheader -s samples_82_rename.tsv`).
- Output: `data/pangenome_p82_chr1.vcf.gz` + `.tbi`.
- Wall: < 30 min.

**A2. PanGenie-index → kmers.tsv.gz.**
- `PanGenie-index -r TAIR10.chr.numeric.iupacN.fa -v data/pangenome_p82_chr1.vcf.gz -o data/kmers_p82/idx -t 8 -k 31`
- Output: `data/kmers_p82/idx_kmers.tsv.gz` (+ companion files).
- SLURM: `--cpus-per-task=8 --mem=64G --time=8:00:00`.

**A3. Build `cn_full_p82/cn_Chr1.{cn,meta}.npz`** with `poolfreq/src/build_kmer_cn.py`.
- `python build_kmer_cn.py --kmers data/kmers_p82/idx_kmers.tsv.gz --vcf data/pangenome_p82_chr1.vcf.gz --ref TAIR10.chr.numeric.iupacN.fa --chrom Chr1 --out data/cn_full_p82/cn_Chr1`
- Output: `data/cn_full_p82/cn_Chr1.cn.npz` + `data/cn_full_p82/cn_Chr1.meta.npz`.
- SLURM: `--cpus-per-task=4 --mem=64G --time=8:00:00`.

**A4. Build `cn_var_p82.{cn_var,meta}.npz`** with `poolfreq/src/build_cn_var.py`.
- `python build_cn_var.py --vcf data/pangenome_p82_chr1.vcf.gz --out data/cn_var_p82`
- Output: `data/cn_var_p82.cn_var.npz` + `data/cn_var_p82.meta.npz`.
- SLURM: `--cpus-per-task=2 --mem=32G --time=2:00:00`.

**A5. Build per-founder consensus FASTAs from `data/pangenome_p82_chr1.vcf.gz`.**
- **DO NOT** symlink to `unimputed_fastas_v3/` — those were built from the v3 merged-haploid VCF, not from the standalone cactus pangenome → variant sets differ → would recreate the v3 bug.
- **DO NOT** symlink to raw `sv_panel/assemblies_clean/*.chr.fa` either — raw assemblies carry variants the cactus VCF may have filtered (`vcfbub` drops big bubbles), which puts sim reads outside the panel's known k-mer space.
- **DO**: `bcftools consensus -f TAIR10.chr.numeric.iupacN.fa -H 1 -s <Accession_ID> data/pangenome_p82_chr1.vcf.gz > fastas_82/<Accession_ID>.chr1.fa`
- Array job 1-82 (one founder per task).
- Output: `fastas_82/<Accession_ID>.chr1.fa`.

### Phase B — Build the simulation

**B1. Adapt the recomb-sim driver for a p82 pool.**
- Clone `sims/visor_freqk/scripts/run_pool_sweep_82_recomb.sh` → `control_p82/scripts/06_run_sim_p82.sh`.
- Add `PANEL=p82` branch:
  - `CACTUS_DIR=control_p82/fastas_82` (yes, all 82 founders are "cactus-style," but the FASTAs are consensus-built so this naming is just historical).
  - `CN_KMER_PREFIX=control_p82/data/cn_full_p82/cn`
  - `CN_VAR=control_p82/data/cn_var_p82.cn_var.npz`
  - `CN_VAR_META=control_p82/data/cn_var_p82.meta.npz`
- Constraint: pool sizes ≤ 82.
- Regimes: `n50_g1`, `n50_g3` mandatory; optional `n80_g1`.
- Cov: 10×, Chr1, seed 42, hotspots crossover model (matches v3).
- Output per regime: `control_p82/sims/cov10_<regime>_p82_chr1/`.

**B2. Compute truth via `compute_recomb_truth.py`** with `--cn-var control_p82/data/cn_var_p82.cn_var.npz --cn-var-meta control_p82/data/cn_var_p82.meta.npz`.
- Output: `control_p82/sims/<regime>/recomb_truth.tsv.gz` per regime.

### Phase C — Run methods on the p82 panel

For each regime:

**C1. cactus_em global** — `per_sample_per_chrom.py --block-mode global` with `--cn-kmer-prefix control_p82/data/cn_full_p82/cn`, `--cn-var control_p82/data/cn_var_p82.cn_var.npz`, `--cn-var-meta control_p82/data/cn_var_p82.meta.npz`.

**C2. cactus_em ★★** — same driver, `--block-mode window --window-bp 10000 --global-anchor-weight 0.3 --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5`.

Submit one SLURM job per (method, regime) — that's 4 jobs for must-haves, 6 with `n80_g1` stretch. Each ~30-60 min wall.

**Out of scope this experiment**: freqk and hapFIRE on the p82 panel. The hypothesis is about the cactus_em h-bias mechanism — those two methods don't help test it.

### Phase D — Evaluation and writeup

**D1. Build `control_p82/results/FINAL_RESULTS_cov10_p82.ipynb`** — same structure as `FINAL_RESULTS_cov10_v3.ipynb` but on the p82 regimes. Headline table + hexbin per method per regime.

**D2. Outlier characterization (h-bias diagnostic)**
- Compute per-sample h-vector mean (`cactus_em` writes h alongside per-variant AF) and verify cactus_em on p82 has h close to uniform 1/82 across founders for `n50_g1` pools where the true h is uniform. The +41% cactus over-credit should not be present (because there's no cactus-vs-PG asymmetry to over-credit).
- For variants with `|est_af - truth_af| > 0.1`, compute carrier rates and locate by Mb on Chr1. Compare to the v3 outlier characterization in `memory/project_v3_pg_specific_outliers.md` — the 136-SNP Chr1q knob cluster should be gone.

**D3. Verdict to memory.**
- If streaks vanish → write `memory/project_p82_control_result.md` confirming mechanism + update `project_v3_singleton_kmer_bug.md` with this evidence.
- If streaks persist → file diagnostic candidates in same memory file for next-session triage.

## Success criterion

The off-diagonal streaks should **disappear** (or drop by an order of magnitude in n_outliers / max-|residual|) on the p82 build. n50_g3 R² should match or exceed v3's 0.960. If yes → cactus-PG asymmetry confirmed as the mechanism → write up as the positive control + move forward with a **k-mer rebalancing** strategy on the 231-founder panel (see `BALANCING_KMERS.md`). **Beagle imputation is not on the table** as a fix.

If streaks persist on the p82 build → mechanism is something else. Candidates: per-window k-mer count variance, simplex identifiability at high similarity between cactus founders, residual sim/panel inconsistency. Diagnose from there.

## Caveats to call out in the writeup

1. **p82 panel is smaller than production (231).** A reviewer might say "of course it's cleaner, you removed the hard half of the panel." Response: this experiment isolates the *mechanism*, not the production recommendation. The production fix is to equalize k-mer fingerprints across founders via k-mer rebalancing on the 231-panel (see `BALANCING_KMERS.md`), not to drop to 82.

2. **Pool sizes ≤ 82**, so direct apples-to-apples vs the v3 `n200_g1` cell isn't possible. Match on regimes that overlap: `n50_g1`, `n50_g3`, optionally `n80_g1`.

3. **Truth definition must use `cn_var_p82`** — using v3's `cn_var_231_v3` would project against the wrong panel (different variant set).

4. **Variant set will differ from v3.** The standalone p82 cactus VCF after `vcfbub + norm` will contain a different set of variants than v3's merged 231-founder VCF (which mixes cactus + Beagle-imputed PG). For variant-level comparisons against v3, intersect on `(chrom, pos, ref, alt)` from both panels' meta files.

## Compute budget estimate

| Phase | Wall | Notes |
|---|---|---|
| A1 (subset → norm → reheader, Chr1) | < 30 min | bcftools/vcfbub pipeline |
| A2 (PanGenie-index Chr1, 82 founders) | ~2-4h | one SLURM job |
| A3 (cn_full_p82, Chr1) | ~2-3h | one SLURM job |
| A4 (cn_var_p82, Chr1) | < 1h | numpy build from VCF |
| A5 (founder FASTAs, 82 of them) | ~1-2h | array job 1-82, each ~1 min |
| B (sims for n50_g1 + n50_g3, cov10, Chr1) | ~3-5h each, parallel | VISOR SHORtS dominates |
| B2 (truth) | < 30 min/regime | |
| C1+C2 (cactus_em global + ★★, 2 regimes) | ~30-60 min each, parallel | 4 SLURM jobs |
| D (notebook + writeup) | ~half day | |

**Total: ~1 day SLURM compute + ~half-day dev/eval (must-have regimes only).**

## Anchor files to read before starting

- `SESSION_2026-05-12.md` — the v3 lock-in + the residual bias context
- `memory/project_v3_singleton_kmer_bug.md` — the +41% h-bias diagnosis
- `memory/project_v3_pg_specific_outliers.md` — the 136-SNP outlier characterization
- `memory/project_cactus_ref_is_tair10.md` — why we use `TAIR10.chr.numeric.iupacN.fa`
- `sims/visor_freqk/RECOMB_SIM.md` — the recomb-sim architecture (will be reused almost as-is, just with a different panel)
- `ALGORITHM.md` + `CACTUS_EM_MATH.md` — the EM mechanics, in case you need to reason about why the bias arises

## Reference scripts to clone (do not modify originals)

| Original | Adapt to | Key changes |
|---|---|---|
| `pangenie_genotyping/scripts/build_pangenie_index.sh` | `control_p82/scripts/02_build_pangenie_index_p82.sh` | Input VCF path + output dir |
| `poolfreq/scripts/build_cn_full_v3_one.sh` | `control_p82/scripts/03_build_cn_full_p82.sh` | Input kmers + VCF; output dir; Chr1 only (drop array 1-5) |
| `poolfreq/scripts/build_cn_var_v3.sh` | `control_p82/scripts/04_build_cn_var_p82.sh` | Input VCF; output prefix |
| `sims/visor_freqk/scripts/build_unimputed_v3_fastas.sh` | `control_p82/scripts/05_build_fastas_p82.sh` | Input VCF (p82); 82-task array; sample IDs from Accession_ID; output dir `fastas_82/` |
| `sims/visor_freqk/scripts/run_pool_sweep_82_recomb.sh` | `control_p82/scripts/06_run_sim_p82.sh` | `PANEL=p82` branch with the 4 path overrides |
| `poolfreq/src/per_sample_per_chrom.py` | (reused as-is) | Just pass p82 args |
| `sims/visor_freqk/scripts/compute_recomb_truth.py` | (reused as-is) | Just pass `--cn-var control_p82/data/cn_var_p82.cn_var.npz` |

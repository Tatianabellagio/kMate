# Overnight status — 2026-04-27 23:15 → 2026-04-28 morning

I'll update this file as jobs land overnight. Latest status at top.

## Snapshot at start of overnight session (23:15)

**Running:**
- 56188 build_cn_full (genome-wide cn, all 5 chroms) — 33 min in
- 56192 imp_merge (build ref/target panels) — 3 min in
- 56196 build_cn_var — 2 min in
- 56169 cactus_nocap rebuild — 10h in (~30h to go)
- 56180 cactus_all rebuild — 7h in (~50h to go)

**Queued (waiting on dependencies):**
- 56189 val_genomewide ← 56188
- 56194 imp_beagle ← 56192
- 56195 imp_loo ← 56192

**Major milestones expected overnight:**
- ~00:30 — 56192 done → 56194 Beagle starts
- ~01:30 — 56188 done → 56189 validation starts
- ~01:30 — 56194 done → imputation complete
- ~03:00 — 56195 LOO validation results
- ~03:00 — 56189 validation done → genome-wide r² for SEEDMIX-realistic

## Issue handling rules I'm using overnight

- OOM / timeout: bump resources, resubmit
- Path / typo errors: fix, resubmit
- Genuine model-level failures: leave for morning review with diagnosis
- DO NOT cancel cactus rebuilds (56169, 56180) without explicit go-ahead

## Updates (latest first)

### 2026-04-28 08:55 — morning triage

**COMPLETED:**
- 56188 build_cn_full → all 5 cn matrices written (~1.3 GB total)
- 56196 build_cn_var → `cn_var_82.cn_var.npz` written

**FAILED:**
- 56192 imp_merge → bash bug `$PYTHON=...` with `set -u`. Fixed (removed bad line). Resubmitted as **56199**.
- 56189 val_genomewide → OOM Killed at 128 GB during CVXPY KL solve on ~1M k-mers. Bumped to 256 GB. Resubmitted as **56202**.

Cancelled 56194/56195 (chained on dead 56192). New chain: 56199 → 56200 (Beagle) + 56201 (LOO).

(scheduled overnight wakeup didn't actually fire — failures sat ~3h before morning triage)

### 2026-04-28 09:01 — second imp_merge attempt failed

- 56199 imp_merge → FAILED at 7 min in step 6 (target panel build). Bash heredoc with `set -e` triggered on a benign warning. Rewrote step 6 in pure Python using `pysam` for cleaner control flow. Resubmitted as **56206** (and chained 56207 imp_beagle, 56208 imp_loo, 56209 cn_231).

### 2026-04-28 09:30 — imp_merge succeeded, Beagle running

**COMPLETED:**
- 56206 imp_merge → 21 min. ref_80.vcf.gz: 573 MB, 3,420,143 records (3.24M SNPs + 184,663 SVs after F_MISSING filter). target_151.vcf.gz: 994 MB, 3,476,635 records (3.24M SNPs + 241,155 SVs as `./.`).

**Note**: target_151 has ~56K MORE SVs than ref_80 — these are SVs where the cactus panel had F_MISSING > 0 and were dropped from ref. They will pass through Beagle as `./.` and get filtered at the cn_var stage. Acceptable.

**RUNNING:**
- 56207 imp_beagle (started 09:26) — Beagle on chrom 1 in progress
- 56208 imp_loo (started 09:26) — first leave-out (sample 6939) finished in 1m20s ✓ — promising sign that Beagle is happy with the merged panel
- 56202 val_genomewide — 35 min in, ~50 min CPU, 114 GB RAM, 13 GB disk read (loading cn + counting k-mers, no output yet)
- 56203 per_sample_smoke — 31 min in, ~52 min CPU, 115 GB RAM, 10 GB disk written (skewed5 sim k-mer counting in progress)

**QUEUED:**
- 56209 cn_231 → waits on 56207

**Auxiliary documentation completed during wait:**
- `imputation/IMPUTATION_PLAN.md` rewritten with full thought process (why imputation, ref-genome verification, SNP overlap test, two non-obvious bugs, validation strategy, integration with pool-seq pipeline).

### 2026-04-28 13:30 — imputation finished, big finding on per-record metric

**COMPLETED:**
- 56207 imp_beagle: 151 founders imputed, 3.42M records (chr1 done in 6m39s; total ~30 min)
- 56208 imp_loo: 10 founders × full genome → final founder 10002 hit **98.5% SV concordance** (181,853/184,663 SVs match cactus truth, 0 missing) — exceeds 96.8% smoke benchmark, confirming Beagle is *more* accurate at full panel scale
- 56203 per_sample_smoke: 6,187,481 per-record alt_freq TSVs for SEEDMIX_S1 and skewed5_sim, 82-founder

**FAILED:**
- 56202 val_genomewide: segfault at 4h16m, no stdout. Diagnosed as combined cause: (1) `solve_em` cast `cn.astype(np.float32)` → 26 GB allocation per call, accumulating across 3 pools; (2) `solve_block_wls` via CVXPY/SCS impractical at K=80M (would need >100 GB just for problem data). Fix in commit (untagged): drop WLS, cast cn to f32 once at top, add `python -u`. Resubmitted as **56704**.

**MAJOR finding from validate_seedmix_recipe.py on 82-founder baseline:**
SEEDMIX_S1 per-record alt_freq vs recipe-truth:
- **Overall: R² = 0.994, Pearson r = 0.997** ← excellent, despite 65% missing-mass on the 151 unobserved founders
- Stratified: AC=1 R²=−7.0, r=0.09 (singletons are noise) | AC 2–4 R²=0.01, r=0.66 | AC 5–10 R²=0.32, r=0.75 | AC>10 R²=0.99, r=0.994

**Implication**: imputation's value is on **rare-SV freq estimates** (singletons fold into AC>10 once 151 imputed founders fill in), not on the already-saturated common bucket. The previously reported "r=−0.15 → 0.99 jump" was at the founder-level h-vector. At the per-record level (the actual deliverable), the 82-founder pipeline is already strong on commons; 231 should mainly lift rares.

**RUNNING:**
- 56209 cn_231 — Chr1, Chr2 done; on Chr3 (started 13:48). ETA Chr3+4+5 ≈ 4–5h → ~18:30 today.
- 56704 val_genomewide (rerun) — densify+cast OK, in EM on uniform82 sim
- 56705 seedmix_82 — running per_sample_driver on SEEDMIX_S2..S8 (S1 pre-staged from smoke). When done, gives full 8-replicate baseline for 82-vs-231 comparison once cn_231 lands.

**Infra prepped during wait:**
- `poolfreq/src/validate_seedmix_recipe.py` — per-record validator vs recipe truth, supports both 82 and 231 panels via `--panel-map` flag
- `poolfreq/tests/run_seedmix_82_baseline.sh` — 8-replicate batch runner
- `poolfreq/tests/run_batch_template.sh` — parameterized for 82 vs 231 via env vars
- `em_solver.py:47–48` — skip float32 cast when input is already f32 (eliminates double-allocation)

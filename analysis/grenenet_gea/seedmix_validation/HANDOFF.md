# kMate founder-h non-identifiability — investigation handoff (2026-07-06)

> **RESOLVED (2026-07-07, refined same week).** Despite this doc's "NEXT AGENT" framing
> below picking **EM-REG (Dirichlet/prior_h) as the front-runner**, that is NOT what won.
> The actual fix came from the **NORM arm**: `fix_norm/em_variants.py`'s "poisson" M-step
> (normalize by each founder's own k-mer content `Kf_w`, not a global count total —
> a missing RNA-seq-style effective-length correction). Validated in
> `../kmate_founder_fix_results.ipynb`: poisson + filt2inv + no-ω beats production 2.7×
> on AF-MAE, 0 vs 33 absorbed founders (drop ω=1/mb, keep filt2inv). It was then integrated
> into the real `solve_em()` and refined one day later — `Kf_w` must sum over the FULL
> k-mer unit, not observed-only k-mers (the `fix_norm` prototype never hit this second bug);
> see `fix_kfw_fullpanel/kmate_kfw_fullpanel_results.ipynb`, commit `9669be7`. `fix_emreg/`,
> `fix_local/`, and `fix_seedmix/` are the abandoned/superseding-candidate arms, kept below
> as decision-trail, not live guidance.

Memory: `kmate-founder-h-nonidentifiability.md` (auto-loads). This is the actionable version.

## THE PROBLEM
kMate GLOBAL-mode EM under-estimates the seed-mix frequency `h` of ~19-49 of the 231 founders
(worst 9977 at 3e-15), while all 231 are present at ~1/231 by design. **Allele frequency is fine
(r~0.98); only the founder DECOMPOSITION h is wrong.** This biases every p0-anchored quantity —
the selection coefficient `s = slope from p0` reads spurious selection on absorbed founders — so the
whole founder-level GEA/GWAS (Family A & B) is affected.
kMate MUST stay standalone: **hapFIRE was only a validation ruler; no fix may depend on it.**

## CAUSE — SETTLED
1. **Pure EM non-identifiability (headline, DIAG arm 1).** Ideal *noiseless* counts `c = K^T·h_uniform`
   fed to the production-panel EM STILL collapse **36 founders to ~0** (9507→5e-31, 9977→4e-11) and
   over-credit cactus 1.58×; eff_n ~155/231, EM did not converge in 400 it. → the collapse is panel
   geometry + EM sliding off the uniform saddle to a sparse vertex, NOT read noise/coverage.
   Fix must be **regularization or panel geometry**; better reads won't help.
2. **Two mechanisms:**
   - (A) absorbed founders carry few *discriminative* (low-`ac`) k-mers; their k-mers are shared by
     ~175-202/231 (9977 = rank1 fewest ac<=5 & highest mean_ac). Deep absorbers (p0<1e-4) are pure A.
   - (B) panel imbalance: 80 long-read(cactus)+151 short-read(PG). LR have ~1.76× more discriminative
     k-mers → raw EM over-credits LR (cactus mass 1.58×). Production `ac==1` drop + `omega=1/m_b`
     bubble weighting exist to fix (B). Some PG founders w/ many discriminative k-mers still absorbed.
   - `src/kmate/filter_kmer_pa_production.py` keeps `2<=ac<=230`; its comment calls ac==1 "zero
     discrimination" — BACKWARDS (ac==1 = max discrimination), but ac==1 is ALSO the LR-over-credit
     source, so can't naively restore.
3. `src/kmate/em_solver.py` `solve_em()` ALREADY has `dirichlet_alpha` (pull→uniform), `prior_h`+
   `prior_weight` (anchor), `omega` — production uses none of the priors.

## YARDSTICK — BUILT
`benchmarks/h_accuracy/score_h_vs_truth.py` scores per-founder h vs sim truth
(`benchmarks/p231/sims/cov10_n231_g0_s42_hotspots_p231_chr1/pool_weights.tsv`), tracks LR/SR mass.
**BASELINE to beat (n231_g0 uniform sim):** production filt2mb **RMSE 0.00332, 27 absorbed,
cactus mass Δ −0.017**. Raw arm: RMSE 0.00388, 38 absorbed, cactus Δ +0.215 (the imbalance).
A valid fix cuts RMSE/absorbed while keeping cactus Δ near −0.017 (don't reintroduce LR over-credit).
Also verify on a NON-uniform truth (n50 sim) that the fix doesn't flatten real signal.

## RUNNING JOBS (land on disk regardless of session)
- DIAG ladder: `diag/reproduce_founder_collapse.py` → `diag/diag_full.out`. Arm 1 done; arms 2-4
  (filt2inv+omega, raw panel ±omega) still computing (~900s each on node).
- **EM-REG (front-runner): SLURM job 35541214 `emreg_sweep` RUNNING** → `fix_emreg/` (Dirichlet-alpha
  + prior_h sweep). Check `squeue -u $USER`; results land in fix_emreg/.
- LOCAL: `fix_local/poc_local_em.py` → `poc.log` (window-estimate + identifiability-weighted aggregate;
  standalone version of why hapFIRE wins; backed by per-chrom within-absorbed Spearman(h,mean_ac)=-0.79).
- NORM: `fix_norm/{em_variants.py,run_experiment.py}` written but run OOM-FAILED — **rerun via sbatch**.

## NEXT AGENT — DO THIS
1. Collect results from `fix_emreg/ fix_local/ fix_norm/` + finish `diag_full.out`. Score each candidate
   h with the yardstick.
2. Pick the winner. **Regularization (Dirichlet / prior_h) is the front-runner** given the pure-
   identifiability result. Open question: regularize WITHOUT erasing selection signal on EVOLVED
   samples → use `prior_h` = previous generation / p0, NOT uniform (uniform trivially fixes the seed
   mix but would flatten real evolved differences).
3. Integrate the winner into `src/kmate/em_solver.py` (currently UNMODIFIED — agents prototyped in
   fix_*/ only), re-run seedmix + a couple evolved samples, re-score.
4. Then revisit the downstream GEA: `analysis/grenenet_gea/build_selection_trait.py` still drops to
   212 founders on `p0>1e-3` (should keep ~230 / anchor p0 at 1/231). Re-derive the selection trait
   with corrected h.

## OPS / GOTCHAS
- k-mer panels are huge (Chr1 ~10.9M k-mers, 461M nnz). Loading multiple concurrently OOMs the 251G
  node (happened this session). **SBATCH heavy k-mer jobs.** Chr1-only for prototyping.
- Env: `conda activate kmate` (compute) / `basic` (matplotlib+nbconvert; matplotlib HANGS in `plotting`).
- Compute-node only (hostname n*.savio*; login-node hook blocks compute). No pip (mamba). No co-author trailers.

## SEPARATE uncommitted edits this session (not the h problem)
- Genomic control REMOVED (raw kinship-corrected z) from `analysis/grenenet_gea/class_split_gwas.py`
  + `founder_gwas_multisite.py` (were self-inflating where per-site λ<1). `class_gwas_*.npz` regenerated
  no-GC (raw λ_JOINT ~1.34 snp/nonsnp, 1.20 sv; SNP↔nonSNP concordance unchanged). Downstream overlap
  tables / Manhattan PNGs now stale. Decision pending: calibrate JOINT/CLIMATE via site-permutation null
  (the code's planned Stage 2), NOT GC.

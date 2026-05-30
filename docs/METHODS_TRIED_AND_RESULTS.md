# Methods tried & final results — authoritative log

**Last updated 2026-05-27.** The current, consolidated record of every method/artifact
we tried and its verdict. Supersedes `old_docs/METHODS_TRIED.md` (taxonomy as of 05-07)
and complements `RESULTS_LOG.md` (chronological numbers) and `archive/README.md`
(what's parked). Status legend: ★ production · ✓ kept · ◯ neutral/optional · ✗ lost · ⋯ in-progress.

---

## 0. Production pipeline (as of 2026-05-27)

kMate (legacy name `cactus_em`): per-sample Poisson k-mer EM on the 231-founder simplex
→ h → project through var_pa → per-record AF.

- **EM kmer_pa (k-mer index):** `kmer_pa_231_v3qc_v3_filt2` (drop ac=1 singletons). ★
- **EM weighting:** `--kmer-weight inv_mb` (ω_k = 1/m_b, per-bubble de-replication), GLOBAL mode. ★ for the heterogeneous 231 panel — but **panel-conditional** (see §3).
- **Projection var_pa:** arch3 — `var_pa_231_arch3_chr1_atomized` (SNP-level GEA) + raw `var_pa_231_arch3_chr1` (SNP/indel/SV). ★
- **Canonical VCF:** `panel/arch3/chr1/merged_231_chr1_final.vcf.gz` (v3qc QC: GQ≥20 + 5772/9947 dropped; arch biallelic; unimputed; het→missing).
- **Driver:** `src/per_sample_per_chrom.py` (`--block-mode global --kmer-weight inv_mb`); solver `src/em_solver.py` (`solve_em(omega=)`).

---

## 1. kmer_pa lineage (k-mer presence/absence matrix $K_{\mathrm{pa}}$) — what each version was & verdict

| Version | What changed | Verdict |
|---|---|---|
| `kmer_pa_231` | v1 first build | ✗ superseded |
| `kmer_pa_231_v2` | PanGenie-index-derived k-mers | ✗ superseded (in-house builder replaced PG-index) |
| `kmer_pa_231_oursHAP/DIP/CAP2X` | in-house build: haploid / diploid / per-allele-cap-2× experiments | ◯ haploid+cap chosen as the build convention; dirs superseded |
| `kmer_pa_231_v3` | in-house `build_kmers_tsv.py`; fixed missing-GT→N handling | ✗ superseded by v3qc |
| `kmer_pa_231_v3_filt2/3/5/10` | singleton-threshold sweep (drop ac<2/3/5/10) | ✓ **filt2 (ac≥2) won** the threshold sweep; filt3/5/10 over-prune |
| `kmer_pa_231_v3qc` | + QC (GQ≥20, exclude 7 assemblies) | ✗ superseded by v3qc_v3 |
| `kmer_pa_231_v3qc_filt2` | v3qc + filt2 | ✗ superseded by v3qc_v3_filt2 |
| `kmer_pa_231_v3qc_lowmiss/nomiss` | missingness-threshold builds | ✗ failed/empty builds (512 B) |
| `kmer_pa_231_v3qc_v2*` (+ rownorm, mixed*) | v3qc rev2 + row-normalization + mixed-loose/strict/conserv k-mer-set experiments | ✗ all superseded by v3qc_v3 |
| `kmer_pa_231_v3qc_v3` | + het→missing masking (May 18); CURRENT base | ✓ base for filt2 |
| **`kmer_pa_231_v3qc_v3_filt2`** | v3qc_v3 + drop ac=1 | ★ **front-runner base** |
| `…_v3qc_v3_mixedloose/strict/conserv` | mixed cactus∪PG k-mer-set rules | ✗ lost to filt2 |
| `…_v3qc_v3_classmatchPG_refilt2` | class-match cactus k-mers to PG count | ✗ lost (see [[project_kmer_imbalance_investigation]]) |
| `…_v3qc_v3_subsampMedian_refilt2` (+seed1/7/100) | per-AC subsample to panel median | ✗ lost to filt2+1/m_b on AF MAE (over-corrects) |
| `…_v3qc_v3_subsampProtect1/2_refilt2` (+seeds) | subsample but protect low-AC discriminating tags | ◯ protect1 best RMSE at dense pools, but lost overall to filt2+1/m_b |
| `…_v3qc_v3_filt2_subsampMedian` | filt2 then subsample | ✗ lost |

---

## 2. k-mer filtering / subsampling strategies (the h-imbalance fight)

Goal: remove the cactus-vs-PG h over-credit (cactus founders over-attract due to more
rare/discriminating k-mers). See [[project_h_imbalance_rootcause]].

| Strategy | Verdict |
|---|---|
| filt2 (drop ac=1 singletons) | ✓ **always on** (removes founder-private sequencing-error tags) |
| subsampMedian (per-AC delete to median) | ✗ over-corrects; lost on AF MAE |
| subsampProtect1 (protect ac=2 tags, subsample high-AC) | ◯ −22% RMSE at dense pools, but not the overall winner |
| effective-n / correlation down-weighting | ✗ no better than filt2 |
| inverse-AC, balanced-SNP, down-wt-rare, per-bubble-cap, uniform-anchor | ✗ all lost (NONE beat filt2; some hurt AF 1.5–2×) |
| **ω_k = 1/m_b global (per-bubble de-replication) on filt2** | ★ **WINNER (heterogeneous panel)** — turns "h∝k-mer count" into "h∝locus count" |

---

## 3. Final method decision & the panel-conditional caveat

**filt2 + ω_k=1/m_b global** wins AF MAE in every regime on the heterogeneous 231 panel
(g0 sims, replicate-validated). DECISIVE result.

**BUT — benchmarks/p80 (homogeneous 80-cactus panel, 2026-05-27):** with no cactus-vs-PG
imbalance, `filt2 + 1/m_b` LOSES to plain `filt2` in every regime (+2% to +41% MAE).
⇒ 1/m_b is NOT a clean de-replication correction; it's a cactus-suppression that only
helps when the over-credit exists to cancel. It down-weights large-bubble k-mers that
carry real per-founder discrimination. So: **1/m_b is opt-in for the imbalanced panel;
plain filt2 is the right default on a balanced/homogeneous panel.** See
[[project_p80_control_1mb_hurts]]. Implementation: `solve_em(omega=)` (omega=None
byte-identical to the old MLE solver).

---

## 4. var_pa lineage (variant presence/absence $V_{\mathrm{pa}}$ for projection)

| Version | Verdict |
|---|---|
| `var_pa_231_v2` | ✗ Beagle-IMPUTED (incl. SVs) — rejected ([[var_pa_231_v2_is_beagle_imputed]]) |
| `var_pa_231_v3`, `_v3qc`, `_v3qc_v2` (+OLD_no_called_mask) | ✗ superseded; pre-arch decomposition |
| **`panel/arch3/chr1/var_pa_231_arch3_chr1`** (raw) | ★ SNP/indel/SV classes (2.62M records) |
| **`…_arch3_chr1_atomized`** | ★ per-base SNP-level (7.46M); use for SNP GEA |

Arch3 = graph-annotated + `convert-to-biallelic` decomposition; lifted hapFIRE SNP
coverage 62.6%→79.8%; atomization closed MNP-vs-SNP encoding outliers (RMSE −23%).
Never use Beagle imputation ([[feedback_no_beagle_solutions]]).

---

## 5. Panels

| Panel | What | Verdict |
|---|---|---|
| 231 heterogeneous | 78 cactus long-read + 153 PG short-read | ★ production; the real GrENE-Net panel |
| p82 control | 82-accession standalone cactus | ✓ control; p82 beat v3 apples-to-apples ([[project_p82_control_result]]) |
| **p80 control** | 80 all-long-read (homogeneous; drop 2 flagged) | ✓ isolates the imbalance effect ([[project_p80_control_1mb_hurts]]) |

---

## 6. Recombination model (sims)

- **OLD:** crossovers forced at hapFIRE BigLD block boundaries (`--crossovers-from-ld-blocks`). ✗ **circular** for benchmarking the block-based method.
- **NEW (benchmarks/p231):** RANDOM crossovers at the A. thaliana rate (4 cM/Mb, Poisson, uniform positions). ★ unbiased. Implementation: drop the LD-block flag.

## 6b. Gen-0 founder sampling — **multinomial-w/-replace is BANNED for g0 sims**

Established 2026-05-27 after comparing the new p231/p80 g0 results to H_FIX_METHODS_COMPARISON's g0 plots, which showed cleaner founder recovery. The H_FIX sims came from `build_g0_uniform_sim.py` which used `rng.choice(replace=False)` — truth_h = 1/n uniform per sampled founder. The new `make_recomb_mosaics_*.py` defaulted to **with-replacement** at gen-0, which injects Poisson(n/F) count quantization into truth_h (some founders count=0/1/2/3+, truth_h ∈ {0, 1/n, 2/n, 3/n, ...}). The resulting "staircase" truth muddies h-recovery comparisons.

**Policy (★ production for sims):** always pass `--gen0-no-replace` to `make_recomb_mosaics_p{80,231}.py`:
- n ≤ F → sample without replacement → truth_h = 1/n for the n sampled founders, 0 elsewhere (this is the SEEDMIX-mimicry case).
- n > F → balanced allocation (`floor(n/F)` or `ceil(n/F)` per founder, evenly distributed) → at most 2 distinct truth_h values.

Applies to BOTH g0 AND g≥1 when the regime models a uniform pool (e.g. n=F at any gen — recomb permutes ancestry within individuals but the chrom-average per-founder count is preserved by the gen-0 draw).

`06_run_sim_p{80,231}.sh` were updated to **always** pass `--gen0-no-replace`. The legacy with-replacement gen-0 path is preserved in the builders (default), but the sim drivers no longer invoke it.

See [[feedback_gen0_no_replace_required]] for the full rationale and do-not-revisit note.

---

## 7. Benchmarks

- `benchmarks/p80/` — homogeneous-panel control. filt2 vs filt2+1/m_b A/B, 6 regimes. DONE; results in `benchmarks/p80/results/filt2_mb_vs_uniform_summary.tsv` + notebook.
- `benchmarks/p231/` — ⋯ **headline 231 benchmark, in progress** (random crossovers; both var_pa arms; kmer_pa rebuilt from arch3 for single-source). See `benchmarks/p231/README.md`.
- Prior eval notebooks: `notebooks/FINAL_RESULTS_cov10*.ipynb`, `panel/arch3/chr1/AF_TRUTH_VS_ESTIMATE_arch3_chr1.ipynb`.

---

## 8. Archive manifest — cleanup 2026-05-27 (~101 GB moved, nothing deleted)

Moved (reversible `mv`) to `archive/`; large data is git-ignored there (on-disk only).

**`archive/kmer_pa_superseded/` (56 GB, 32 dirs):** every `kmer_pa_231_*` EXCEPT the
two kept in place — `kmer_pa_231_v3qc_v3` (front-runner base) and
`kmer_pa_231_v3qc_v3_filt2` (front-runner). Archived: v1, oursCAP2X/DIP/HAP, v2, v3,
v3_filt2/3/5/10, v3qc, v3qc_filt2, v3qc_lowmiss/nomiss, v3qc_v2 (+rownorm/mixed*),
v3qc_v3_{classmatchPG, filt2_subsampMedian, mixed*, subsampMedian(+seeds),
subsampProtect1/2(+seeds)}. Verdicts in §1–2.

**`archive/var_pa_superseded/` (3.3 GB):** `var_pa_231_{v2,v3,v3qc,v3qc_v2}*`
(+ `.OLD_no_called_mask`). Kept: `var_pa_231_v3qc_v3` (data) and the
production `panel/arch3/chr1/var_pa_231_arch3_chr1*`.

**`archive/fastas_superseded/` (42 GB):** `unimputed_fastas_v3` (pre-QC),
`founder_fastas_231_v3`, `imputed_fastas_v2_DEPRECATED`. Kept: `unimputed_fastas_v3qc`
(current) + `benchmarks/p231/fastas_231` (rebuilt from arch3).

**Held (not archived):** old `pool_sweep_82_recomb*` p231 sims — archive only after
`benchmarks/p231` sims validate (they supersede these). The earlier 65 GB of
`*.v2_imputed_BACKUP / OLD_BUGGY / DELETE_ME / macbad` dirs were DELETED 2026-05-27
(explicitly stale-marked, current counterparts verified).

**Dangling refs (recoverable):** `scripts/{build_kmer_pa_*,run_seedmix_*}`
for the archived variants now point at archived dirs — they are the build/run scripts of
the concluded experiments; not archived (kept as records), will fail only if re-run.

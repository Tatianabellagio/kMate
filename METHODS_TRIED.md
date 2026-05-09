# Methods tried — what we explored, what won, what didn't, and the literature that grounds each

**Last updated 2026-05-07.** Companion doc to `RESULTS_LOG.md` (which is the chronological numerical record). This file is the *taxonomy* — read it to get a map of the design space we've covered, what our chosen recipe is *vs.* its alternatives, and how each piece relates to existing published work.

For headline numbers see `FINAL_RESULTS.ipynb` (cov50) and `FINAL_RESULTS_cov10.ipynb` (cov10). For the per-sweep reasoning see `RESULTS_LOG.md` 2026-05-06 entries.

---

## Status legend

| Symbol | Meaning |
|---|---|
| ★ | **Production** — current shipped recipe |
| ✓ | Tried, kept (component of the production recipe) |
| ◯ | Tried, neutral or marginal — kept as option |
| ✗ | Tried, lost — discarded with postmortem |
| ⋯ | In-progress / parked |

---

## 1. Per-sample inference architectures

### 1.1 Per-individual genotype calling (NOT our problem)

Listed for context only — this is the "wrong tool" baseline.

| Method | What it does | Why it's not us |
|---|---|---|
| PanGenie (Ebler 2022, *Nat Genet*) | Pangenome k-mer index + per-individual diploid inference | Per-individual — no pool-seq notion of fractional ancestry |
| KAGE / KAGE 2 (Grytten 2022/2023) | Alignment-free k-mer pangenome genotyper, Bayesian variant-correlation model | Per-individual; their variant-correlation idea is borrowed in spirit by our HMM smoother |
| PanTax (Du 2025, *Genome Research*) | Pangenome-graph + abundance estimation for metagenomics | Strain-level metagenomics, not founder pool-seq |

We use PanGenie's k-mer index format (its `kmers.tsv.gz` output is the source of `cn_full_231_v2`), but we run our own EM on top.

### 1.2 Pool-seq frequency estimation — read-level EM (the "HARP family")

| Method | Status | Notes |
|---|---|---|
| **HARP** (Kessner 2013) | wrapped via hapFIRE | base-level multinomial EM with QV; alignment-required; ~30 min/block on Chr1 |
| hapFIRE (xwu, unpublished) | wrapped | HARP per coarse LD block + analytical projection; **prior production** |
| **Microhaplotype pool method** (Tang & Anderson 2023, *BMC Bioinformatics*) | not tried | Read-level EM on 2-5-SNP "microhaplotypes"; closest published method to our exact problem. Worth re-reading if we ever want to derive a principled phase prior |
| PoolHap (Kessner & Novembre 2010) | not tried | Older read-level EM; superseded by HARP |
| HAF-pipe (Tilk 2019) | not tried | Fixed window per-SNP EM; no read-level phase. Conceptually a strict subset of `cactus_em window` |

### 1.3 Pool-seq frequency estimation — k-mer-level EM (the "cactus_em family")

This is our addition. No published precedent for read-level/k-mer-level EM on a *founder simplex* (pool-seq pop-gen); the architecture is canonical in metagenomics but new for our problem.

| Method | Status | Where it lives |
|---|---|---|
| **cactus_em global** | ✓ baseline | single chrom-wide h; `--block-mode global` |
| **cactus_em window 200 kb** | ✓ prior production | `--block-mode window --window-bp 200000` |
| **cactus_em window 100 kb / 10 kb (no anchor)** | ✗ underperforms | per-window EM at fine scale is too noisy without regularization (n50_g3 SNP R² 0.96 → 0.90 going 100 kb → 10 kb) |
| **cactus_em window 10 kb + λ-anchor (Route 1)** | ✓ keystone | `--global-anchor-weight 0.3` — Dirichlet pseudocount toward chrom-wide h_global |
| **cactus_em window 10 kb + post-EM HMM smoothing** | ✓ keystone | Li-Stephens-style across-window smoothing (`--hmm-smooth-passes 5 --hmm-smooth-alpha 0.5`) |
| **★ window 10 kb + anchor 0.3 + smooth 5α0.5** | **★ production winner** | The stack of the two above. See `FINAL_RESULTS.ipynb` |
| cactus_em bigld_panel | ◯ tried, mediocre | EM on BigLD blocks using PanGenie bubble-allele k-mers (16/32 cap). R² 0.84-0.87 SNP; loses to window_200kb |
| cactus_em clean_cc5 (block-haplotype EM) | ◯ tried | Reconstruct full-block FASTA per founder, dedup at class level, k-mers carried by ≥5 founders |
| cactus_em clean_smooth | ◯ tried | clean_cc5 + HMM smoothing (passes=10, α=0.2) |
| cactus_em dose_aware | ◯ tried | clean_cc5 with k-mer multiplicity (count, not 0/1) |
| **bigld_haplotype** | ⋯ in-progress (2026-05-03) | Per-block sequence-unique founder dedup; ~1000× more discriminative k-mers than bigld_panel; full-Chr1 cn build pending |

---

## 2. Three "routes" we explored on top of the cactus_em window baseline

Documented so we don't re-try them. All numbers from n50_g3 cov50 SNP R² @ 10 kb; full sweeps in `RESULTS_LOG.md` 2026-05-06.

### Route 1 — global-anchor KL prior  ✓ KEPT

**Idea.** Per-window EM is high-variance in thin-evidence windows (sparse panel coverage, low cov, repetitive). Add a Dirichlet pseudocount centered on chrom-wide `h_global` so high-evidence windows still adapt locally but low-evidence windows fall back to global.

**Math.** `h_new[f] ∝ h[f]·(cn @ cw) + λ·N·h_global[f]`, then renormalize (`em_solver.py:78-99`).

**Result.** λ=0.3 lifts n50_g3 SNP R² from 0.902 → 0.963 at 10 kb (+6.1 pp). U-shape across λ ∈ {0, 0.1, 0.3, 0.5, 1.0}; sweet spot 0.3-0.5. See `RESULTS_LOG.md` 2026-05-06 "Route 1" entry.

**Literature analogs.** Closest in spirit is **StrainFacts** (Smith 2022, biorxiv) — Bayesian regularization for strain inference, but they regularize on genotypes not on a previously-estimated h. KAGE 2 (2023) uses variant-correlation priors per-individual. The exact "anchor a per-window simplex EM toward a chrom-wide simplex EM" appears to be novel for pool-seq founder freqs; no need to credit.

### Route 2 — overlapping windows  ◯ NEUTRAL

**Idea.** Disjoint 10 kb windows mean each k-mer enters exactly one window. Overlap by step < window_bp so each k-mer enters multiple windows; per-record AF = inverse-distance-weighted average over covering windows.

**Result.** Marginal (≤ 0.1 pp on R²); the projection averaging over-smooths and slightly attenuates real signal. See `block_em.project_blocks_to_records_overlap`.

**Literature analogs.** Sliding/overlapping windows are standard in GWAS haplotype tests but not in HAF-pipe/hapFIRE. Easy to defend as "creeping window approach analogous to GWAS sliding-window haplotype tests" if asked.

### Route 3 — kallisto-style read-level pseudoalignment + EC EM  ✗ LOST (at fine scale)

**Idea.** Per-k-mer Poisson EM throws away within-read phase. A read carries ~120 linked 31-mers from one founder; pseudoalign by intersecting the per-k-mer "compatible-founder" sets across the read's k-mers, group reads by identical compatibility set into equivalence classes (ECs), run EM on EC counts per window.

This is **literally kallisto** (Bray 2016) with founders-as-transcripts. The EC counts are sufficient statistics for the multinomial; per-iter cost is O(unique_ec) instead of O(reads) or O(kmers).

**Result.** **6 pp WORSE than per-k-mer Poisson at the same resolution.** R²=0.842 SNP at 10 kb (cov50_n50_g3) vs `window_10kb`'s 0.902. Per-k-mer Poisson uses BOTH high- and low-count k-mers (zero counts constrain `h^T cn[:, k]` downward); EC-multinomial only uses observed reads → loses absence signal + count gradient. Confirmed across coverage scaling 0.4% → 100% reads. Postmortem in `RESULTS_LOG.md` 2026-05-06 "kallisto-EM" entry.

**Implementation.** `poolfreq/src/per_sample_kallisto_em.py`, `poolfreq/scripts/precompute_panel_u64_index.py`. Kept on disk as a reference implementation; not in the production recipe.

**Literature analogs.**
- **kallisto** (Bray 2016) — original equivalence-class EM for transcript abundances.
- **EMIRGE** (Miller 2011) — full-length 16S reconstruction via read-level EM.
- **GRAMMy** (Xia 2011) — k-mer/alignment EM for metagenomic abundance.
- **Mora** (Zheng 2024, *BMC Bioinformatics*) — read re-assignment for closely-related references with EM + set-cover regularization. Closest in spirit to our problem (similar references → under-determined). Not yet revisited.
- **PanTax** (Du 2025, *Genome Research*) — pangenome-graph strain profiling at metagenomic strain level; conceptually closest.

**Why we'd revisit.** EC-EM is the only architecture that scales O(unique_ec) per iteration — kallisto routinely runs >10⁸ reads in minutes on one core. If we ever need to push to whole-genome 5-chrom 8-kb-window scale on cluster-wide concurrency, the kallisto-style index becomes attractive again *despite* the absence-signal loss. Worth re-examining with: (i) higher coverage where absence-signal contributes less; (ii) stacking with Route 1 anchor on EC EM; (iii) Mora-style set-cover regularization to handle the under-determinedness.

---

## 3. Smoothing variants

| Variant | Status | Notes |
|---|---|---|
| **post-EM Li-Stephens HMM (5 passes, α=0.5)** | ★ production | Per-window h_b ← α·h_b + (1−α)·exp(−recomb·gap)-weighted neighbor average. Adds +0.3-2 pp on top of Route 1. `block_haplotype_em.smooth_h_across_blocks` |
| HMM smoothing (10 passes, α=0.2) | ◯ tried | Slightly over-smooths; gives back ~0.1 pp |
| HMM smoothing (20 passes, α=0.1) | ◯ tried | More over-smoothing; gives back ~0.15 pp |
| HAFpipe-style linear interpolation between adjacent windows | ✓ on by default for non-HMM-smoothed, off when HMM is on | `project_blocks_to_records(smooth=...)`. Auto-disabled when `--hmm-smooth-passes > 0` to avoid double-smoothing. |
| **AF-space spatial Gaussian smoothing** | ✗ catastrophic | R²=0.95 → 0.01 when smoothing global est_af in AF-space. Smooth in h-space, never in AF-space. Memory: `feedback_no_smoothing_on_global_af.md` |

---

## 4. Block-partitioning schemes

| Partition | Status | When to use |
|---|---|---|
| **fixed window 10 kb** | ★ production | Used by the winning recipe |
| fixed window 200 kb | ✓ prior production | Coarse but very robust on n50_g3 (0.973 SNP R²); kept as fallback |
| fixed window 100 kb | ◯ tried | Intermediate; not used |
| BigLD panel-derived blocks | ◯ tried | Median 963 bp / 7 kb mean on Chr1; signal-poor under PanGenie's 16/32 cap |
| Gabriel-style adjacent-r² LD blocks | ◯ tried | `compute_ld_blocks_gabriel` |
| Complete LD partition (hapFIRE-style, fully-independent blocks) | ◯ tried | `compute_ld_blocks` |
| **bigld_haplotype** (sequence-unique founder dedup per BigLD block) | ⋯ in progress | Should give ~1000× more discriminative k-mers per block; full Chr1 cn pending |

---

## 5. Pre-/post-processing decisions

> **Hard rejections — do NOT re-test these without new evidence.** The rows below
> marked ✗ are decisions we have audited and locked. The full numerical evidence
> is preserved on disk (paths in the table); the workdirs themselves were
> archived to free space, but the conclusion is documented. Before proposing to
> revisit any of these, read the linked README first.

| Decision | Status | Reason / where the evidence lives |
|---|---|---|
| Use cactus pangenome (vs syri) for the SV panel | ✓ locked 2026-04-27 | freqk drops 88% of syri SVs at the index step; cactus retains 100% of bubbles |
| Use `vcfbub -l 0 -r 100000` on raw cactus VCF | ✓ | Recovers 19 pp more syri SVs vs the cactus pre-filtered output |
| Cap cactus contig length (`--maxLen`) | ◯ irrelevant | maxLen never gates SV size in this panel — verified by 42-h nocap rebuild |
| **Beagle 5.5 imputation of merged 231-founder VCF** | **✗ HARD REJECT for SVs** | **−29.9 pp on small_sv, −24.7 pp on medium_sv vs the pre-imputation catalog (75-ecotype LOO; 0/75 ecotypes improved post-imputation). Independently corroborated: Crysnanto 2024 (cattle pangenome) + Roboubi 2026 (French dairy SV imputation) + Yan 2024 (rare-disease SV imputation review). Locked: keep PanGenie SVs unimputed for production. Full numerical evidence: `RESULTS_LOG.md` 2026-05-03 (per-class concordance table) + 2026-05-06 (literature lock-in). Production decision: `pangenie_genotyping/data/merged/README_GOLDEN_STANDARD.md`** |
| **KMC 3 / rapidgzip k-mer counting** | **✗ HARD REJECT, jellyfish stays** | **KMC: 80% memory regression (23 → 41 GB) for ~10% speedup. rapidgzip: ~7% speedup once `jellyfish -t 8` already saturates CPU. Earlier "1.45–1.88×" rapidgzip claim was a CLI bug (`rapidgzip -d -c R1 R2` only processed one file). Full table + postmortem: `HANDOFF.md` 2026-05-02 entry. Next real lever (if revisiting per-sample wall) is SSHash or LP-MPHF — see § 6 below, NOT a re-test of KMC** |
| Per-chrom driver vs genome-wide driver | ✓ per-chrom | 3× lower peak RAM (137 GB → 39 GB), 41% faster, bit-identical output. Unlocks 64 GB memex nodes. `per_sample_per_chrom.py` |
| Single calibration linear fit on 231-eco panel | ✓ | 1.43× scale bias from imputation noise; calibrated once via SEEDMIX rep S1 |
| posterior_combine of hapFIRE-proj + freqk | ◯ implemented but unused | weight on hapFIRE ≈ 0.997-1.000 across 23 sim cases; combined MAE ≡ hapFIRE-proj MAE. Useful only on real evolved samples where the recomb assumption is violated and freqk's per-SV variance becomes competitive |
| Polymorphic-only filter (`0 < truth_af < 1`) for the headline R² | ✓ default | `build_session_dataset.py:397`. Two effects: (a) drops records where the regime's founder-subsample has no carrier (truth=0); (b) drops records fixed in the pool (truth=1). Without the filter, R² inflates by 0.013-0.022 because methods correctly predict ~0 at truth=0 records — easy points. The verdict is unchanged; the polymorphic-only filter is just the harder, more honest test. Quantified in FINAL_RESULTS.ipynb "Evaluation coverage" cell |
| Evaluating on n50 regimes (vs full 231) | ✓ stress-test by design | n50 pools 50 of 231 founders → 47-55% of Chr1 records become truth=0 (alt unique to un-sampled 181 founders) and drop. n200 keeps 200 founders → only 20% drop. **Real evolved GrENE-Net samples have all 231 founders contributing**, so the truth=0 fraction is small there; we expect to evaluate ~700-800k Chr1 SNPs (vs ~377k on n50_g3). Documented in FINAL_RESULTS "Evaluation coverage" cell |

---

## 5b. Panel / cn_var version history (data, not algorithms)

For full audit-trail when someone wonders "why is there a v2?". Older builds
were superseded as the founder set, the pangenome graph, or the build code
changed; the current production cn_var is `poolfreq/data/cn_var_231_v2.*`.

| panel | era | superseded because |
|---|---|---|
| `cn_var_82.{cn_var,meta}.npz` | hapFIRE-projection era (pang_82 cactus, 82 founders) | replaced by 231-founder build once 151 PanGenie genotypes landed (2026-05-03 merge) |
| `cn_var_82_renamed_to_1001g.*` | same era, ID-renamed for compatibility checks | downstream tools now consume 1001G IDs natively; renaming step retired |
| `cn_var_231.{cn_var,meta}.npz` (v1) | first 231-founder build | rebuilt as **v2** with the corrected Beagle-imputed founder FASTAs; v1 had a per-record carrier-status bug fixed in v2 |
| `cn_var_decomposed_chr1.*` | early Chr1-only smoke test | replaced by full-genome v2 |
| **`cn_var_231_v2.{cn_var,meta}.npz`** | **★ production** | current |

The four superseded matrices total ~350 MB; on 2026-05-08 they were moved out
of `poolfreq/data/` into `archive/cn_var_old/`. Active code paths
(`per_sample_per_chrom.py`, the cactus_em SLURM templates) hard-code v2 only.
A handful of historical references in `notebooks/run_plots.py`,
`poolfreq/tests/run_seedmix_*`, and the archived hapFIRE-projection notebook
still mention older names — those are historical / smoke-test paths, not
production.

**`cn_var_231_v2` IS Beagle-imputed** — see project-memory note
`project_cn_var_231_v2_is_beagle_imputed.md`. This is intentional: the cactus_em
sim pipeline uses the imputed cn_var because per-sample driver was validated
against it. **The 2026-05-03 RESULTS_LOG decision to NOT impute SVs applies to
the production pool-seq deliverable** (`founders_231_chr.vcf.gz`) — a
different artifact. Both the imputed cn_var (sim-validation) and the
unimputed VCF (production) are correct; they serve different consumers.

---

## 6. Things we identified but **have not tried** (parked, with literature pointers)

| Idea | Pointer | Why it's parked |
|---|---|---|
| **Read-level EC-EM with Route 1 anchor** | kallisto + StrainFacts | Could rescue Route 3 if absence-signal is the only thing it lost. ~1-2 days of glue: wire `cn_full` into a kallisto-compatible index format |
| **Mora-style set-cover regularization** | Zheng 2024 BMC Bioinformatics | Designed for exactly our pathology (similar references, under-determined). Worth a focused half-day if we revive Route 3 |
| **KAGE 2's variant-correlation Bayesian prior** | Grytten 2023 | Per-variant correlation prior; conceptually parallel to our HMM smoother but at the variant level, not the founder-h level |
| **StrainFacts gradient-based fuzzy genotypes** | Smith 2022 biorxiv | If we ever want to make our solver gradient-based / GPU-friendly. 100× speedup over predecessors |
| **Microhaplotype pool method** | Tang & Anderson 2023 | Closest published method to ours; their derivation is useful if we want to add a principled within-read phase prior beyond what HMM smoothing approximates |
| **LD-aware k-mer reweighting** (the "fix freqk" idea) | — | For each SV, find 50 nearest tight-LD SNPs in the panel; use hapFIRE per-SNP freq as expected AF for SV-flanking k-mers; drop k-mers that disagree (off-target hits). Inverse of HAF-pipe. Listed in HANDOFF.md as option 4 of the original 5-step plan |
| **LP-MPHF + custom Rust k-mer counter** | various | 3-10× counter speedup. Defer until per-chrom + memex concurrency proves insufficient at full GrENE-Net scale |
| **SSHash specialized k-mer dictionary** | Pibiri 2022 | ~10× lookup speedup. Same trigger as above |
| **Recomb-aware HMM transition tuning** | Li & Stephens 2003 | Current smoother uses flat α; principled λ tuned to AT 4 cM/Mb may close the remaining 0.6 pp on n50_g3. See HANDOFF "What's pending" #4 |
| **Coverage-saturation curve** | — | We assert "panel-EM saturates above ~5×" but haven't measured. cov∈{1,2,3,5,10,20,50}× on n50_g3 sweep. See HANDOFF "What's pending" #5 |

---

## 7. Citation cheat-sheet

For when we write any of this up. (Items 1–7 are the ones we wrap or directly compare against; the rest are conceptual or methodological context.)

1. **HARP** — Maximum likelihood estimation of frequencies of known haplotypes from pooled sequence data. Kessner et al. 2013, *MBE*.
2. **hapFIRE** — Wu et al., GitHub repo (unpublished as of 2026-05-07).
3. **PanGenie** — Pangenome-based genome inference. Ebler et al. 2022, *Nat Genet*.
4. **HAF-pipe** — Accurate haplotype frequency estimation from low-coverage pool-seq. Tilk et al. 2019.
5. **VISOR** — Variant simulation in single haplotypes for read pools (HACk + SHORtS modules used here).
6. **kallisto** — Near-optimal probabilistic RNA-seq quantification. Bray et al. 2016, *Nat Biotechnol*.
7. **Beagle 5.5** — Phasing + imputation. Browning et al.
8. **Microhaplotype pool method** — Tang & Anderson 2023, *BMC Bioinformatics*.
9. **PoolHap** — Kessner & Novembre 2010.
10. **EMIRGE** — Miller et al. 2011.
11. **GRAMMy** — Xia et al. 2011.
12. **Mora** — abundance-aware metagenomic read re-assignment. Zheng 2024, *BMC Bioinformatics*.
13. **PanTax** — strain-level pangenome-graph profiling. Du et al. 2025, *Genome Research*.
14. **KAGE / KAGE 2** — Grytten et al. 2022 *Genome Biology* / 2023.
15. **StrainFacts** — scalable strain inference. Smith 2022, biorxiv.
16. **Crysnanto et al. 2024** — Cattle pangenome workflow. *Genome Research* 34:300, PMC10984387.
17. **Roboubi et al. 2026** — French dairy cattle SV imputation benchmark. bioRxiv.
18. **Yan et al. 2024** — Pangenome graphs in rare-disease SV analysis. *Nat Commun*.
19. **k-mer-based approaches bridging pangenomics and population genetics** — review, *MBE* 2025.
20. **Li & Stephens 2003** — modeling LD/recomb (basis for HMM smoothing).

---

## How to add to this doc

When you try a new method, before discarding or shipping it, add a row to the relevant section above. For each entry:

- **What it does** — one or two lines, math notation OK.
- **Status** — ★ / ✓ / ◯ / ✗ / ⋯.
- **Result** — headline number(s) on a named regime, ideally a single R² or MAE; full numbers go in `RESULTS_LOG.md`.
- **Where the code lives** — file path or commit ref.
- **Literature analog** — single sentence + citation; if novel, say so.

Add to `RESULTS_LOG.md` for the chronological reasoning; this file is the *map* across the design space, not the lab notebook.

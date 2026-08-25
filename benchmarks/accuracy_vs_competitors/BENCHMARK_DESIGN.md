# kMate vs competitors — benchmark design & rationale

> ⛔ **THE SHARED-PANEL hapFIRE COMPARISON IS RETIRED AS INVALID (2026-08-25).**
> Everything in this doc that scores kMate against hapFIRE **on one shared SNP panel**
> — §2, §2.1, §2.2, and the §5 head-to-head — has been **archived to
> [`../archive/shared_panel_hapfire_comparison_invalid/`](../archive/shared_panel_hapfire_comparison_invalid/README.md)
> and must not be cited.** Numbers are left inline below only so the archived reasoning
> stays readable; treat every kMate-vs-hapFIRE row as withdrawn.
>
> **Why:** kMate and hapFIRE do not take the same kind of input, so one shared panel is
> not a shared input. hapFIRE errors on missing GTs, so `build_shared_snp_panel.sh`
> imputes `MISSING '.' -> 0|0` — handing hapFIRE **false REF homozygotes** at 7.1% of
> records. kMate never eats that: its `var_called` mask divides exactly those sites out
> of the MAR projection. One tool may say "unknown", the other is forced to guess REF.
> (The handicap is missingness/imputation, **not** phasing — the shared panel *is* fully
> phased. Don't re-diagnose this as a phasing bug.) §2.1's "the winner flips with the
> truth convention" was the symptom of this, not a separate issue.
>
> ✅ **The live, legitimate comparison is [`../speed_vs_hapfire/`](../speed_vs_hapfire/)**
> — matched pools, **tool-native panels** (kMate on arch3 SNPs, hapFIRE on greneNet
> SNPs), 5 seeds × {2,5,20,50,231} founders × {1,5,10}× coverage. It points the other
> way at n231/cov10 (hapFIRE RMSE 0.0018 vs kMate 0.0048), which is expected: kMate's
> own value is stable across both designs (R² ≈ 0.9993); what moves is hapFIRE, because
> the shared panel injected sites its native pipeline never sees.
>
> **Still valid in this doc, unaffected by the above:** §3.1 (hapFIRE/HARP structurally
> cannot estimate SVs — an argument about HARP's per-base likelihood, panel-independent)
> and §4 (speed). The vg comparison is not subject to this specific handicap (vg
> consumes the graph, not a phased VCF) but was entangled in the same 4-tool apparatus;
> rebuild it standalone if wanted.
>
> Terminology note: "global (chromosome-wide) mode" = `--unit chrom` (selfing
> production, CLI default); "window / LD-block mode" = `--unit bp` / `--unit ld`; the
> old `--block-mode global|window` are deprecated aliases.

> ⛔ **Before scoring any AF comparison, read [`../SCORING_RULES.md`](../SCORING_RULES.md).**
> Two repeat-prone mistakes (a `(chrom,pos)`-only join that mispairs multiallelic alleles,
> and the MAR truth being a closed loop that flatters kMate) once faked a "kMate worse at AF"
> result. Use `scripts/score_snp_fair.py`; the old `score_competitors.py` /
> `score_competitor_snp.py` are quarantined (they refuse to run). **SNP AF is a tie with
> hapFIRE** — headline kMate's real wins (speed, alignment-free, SVs, ecotype/h resolution).

2026-06-19. How we benchmark kMate against the two relevant competitor classes, the
fairness decisions behind each panel, and the exact claims we can defend. Companion
to the scripts in `benchmarks/accuracy_vs_competitors/scripts/` and the speed runs in
`benchmarks/speed_vs_hapfire/` (`benchmarks/archive/speed_benchmark_abandoned/` is an
earlier, incomplete precursor — vg was never run there — kept only for provenance).

---

## 0. The competitive landscape — two classes, two claims

kMate's value proposition has two separable halves; each needs a *different* comparator.

- **Parity claim (SNPs):** kMate is as accurate as the established founder-aware pool-seq
  AF tools → comparator = **hapFIRE** (a HARP wrapper; the GrENE-net phase-1 tool).
- **Novelty claim (indels/SVs + speed):** kMate estimates SNP+indel+SV frequencies in a
  single alignment-free pass → comparator = **vg giraffe** (the only pangenome-native
  alternative; AF from graph read-support).

There is **no pangenome-native *pool-seq* AF estimator** in the literature (vg/Giraffe
genotypes SVs per-individual then aggregates — never from a pool). That gap is kMate's
novelty.

Competitors we deliberately do NOT run: bare HARP / HAF-pipe (redundant — hapFIRE wraps
HARP; HAF-pipe is Drosophila-tuned), PoolHapX (de-novo, no founder panel), grenedalf /
PoPoolation2 (compute stats *from* AF, don't estimate it).

---

## 1. Panel: p80 first (all long-read founders)

We run the whole comparison on the **p80** panel (80 cactus long-read-assembly founders),
not the production p231 (78 cactus + 153 PanGenie), because:

1. The mixed 231 panel's imbalance is *our* problem; kMate was designed for all-long-read
   founders, so p80 is the clean home-turf test.
2. p80 SNP missingness ≈ 2.2% (vs 4.5% on 231). Low missingness keeps hapFIRE's imputed
   input ≈ kMate's input — on 231 the heavy missingness would diverge the inputs and we
   couldn't tell algorithm from imputation apart.
3. The p80 cactus graph (`pang_1001gplus_82acc.gbz`) doubles as vg giraffe's input, so all
   three tools share ONE panel.

---

## 2. Accuracy — fair design

**One shared SNP panel + one shared truth, scored identically for every tool.**

- **Shared SNP panel** (`work/shared_snps_p80_Chr1.vcf.gz`, 1,391,487 biallelic Chr1 SNPs):
  extracted from the p80 founder VCF, one record/pos, haploid→phased-homozygous diploid,
  missing→REF (hapFIRE needs a complete phased panel). Both tools estimate AF for *these*.
- **Truth** (`build_truth.py`): founder-genotype-weighted alt-AF from the *un-imputed*
  GTs + pool weights, computed straight from the VCF (NOT via kMate's `var_pa`) → a
  tool-independent yardstick. Emits two conventions (see §2.1).
- **Tools, same pool, same sites:** kMate (`var_pa` output subset to shared SNPs);
  hapFIRE (shared VCF + sim BAM aligned to TAIR10 + `TAIR10.chr.iupacN.fa`); vg
  (giraffe→surject→allele-depth AF at shared sites).
- **Metrics:** MAE + RMSE + R² (kMate convention) AND Pearson r (freqk-paper convention).

### 2.1 The missing-data truth convention (decides the "winner" — read this)

The founder FASTAs are built `bcftools consensus -H 1`, which puts **REF at missing GTs**
(`benchmarks/p80/scripts/05_build_fastas_p80.sh:53`). So the simulated reads physically
contain the **missing-as-REF** allele frequency. Therefore:

- `truth_af_phys` (missing→REF) = the **physical** AF actually in the reads → the CORRECT
  ground truth for any read-based estimator. **Primary metric.**
- `truth_af` (MAR; AF among *called* founders) = kMate's `var_pa` estimand. Scoring against
  MAR is a subtle closed-loop that flatters kMate and unfairly sinks hapFIRE/vg.

The winner flips entirely with this choice (driven by a ~1% high-missingness SNP tail):

Refreshed 2026-08-25 on `cov10_n231_g0_s42_hotspots_p80_chr1` (n=1,391,485 panel SNPs;
fully-called n_called==F=80 → n=1,186,489):

| truth | kMate R² | hapFIRE R² |
|---|---|---|
| MAR (kMate estimand) | 0.9981 | 0.8802 |
| **physical** (missing→REF) | 0.8739 | **0.9969** |
| **fully-called only** (n_called==F; convention-free) | **0.9985** | 0.9978 |

Against the physical truth hapFIRE matches the phase-1 paper (R²≈0.99 across 2–150 ecotypes).
**RULE:** headline the **physical + fully-called** numbers; never headline MAR.

The refresh did **not** weaken this caveat — the winner still flips completely with the
convention (kMate 0.998→0.874, hapFIRE 0.880→0.997), so scoring against MAR remains a closed
loop that must never be headlined. What changed is only the convention-free tie-break: on the
fully-called basis kMate now edges hapFIRE (0.9985 vs 0.9978) where it previously trailed
(0.995 vs 0.998). Pre-refresh values were 0.994 / 0.867 / 0.995 for kMate.

### 2.2 SNP results (p80, n231_g0, `--unit chrom`, fully-called) — REFRESHED 2026-08-25

Basis: `info ≥ 0.99` (the convention-free fully-called cut, §2.1), n = 1,188,444 SNPs.

| tool | type | MAE | R² | Pearson r |
|---|---|---|---|---|
| **kMate** | founder-aware | **0.0036** | **0.9993** | 0.9996 |
| hapFIRE | founder-aware | 0.0057 | 0.9979 | 0.9991 |
| vg-giraffe | founder-naive coverage | 0.0672 | 0.656 | 0.857 |

**kMate now leads hapFIRE at g0** (MAE 1.6× lower). Pre-refresh these were kMate 0.0093 /
0.995 vs hapFIRE 0.0061 / 0.998 — i.e. the old "hapFIRE is marginally the most accurate SNP
tool" conclusion was an artifact of the pre-2026-07-07 estimator (observed-only `Kf_w`), not a
property of the methods. Both founder-aware tools remain ≫ founder-naive coverage (vg). vg is
**unbiased but ~10× noisier** at 10× depth (per-site coverage AF; verified not a bug — coverage
9.3×, per-bin calibrated; QC filtering doesn't help).

⚠️ Scope: this is **one pool, one seed, at g0** (an unrecombined founder mixture — the regime
`--unit chrom` is exactly right for). It does **not** generalize to recombinant pools, where
hapFIRE still wins — see §5. State the regime whenever quoting this row.

**Coverage panel (p80 n80_g0, fully-called SNP R²):**
| tool | 10× | 50× |
|---|---|---|
| **kMate** | **0.9993** | **0.9993** |
| hapFIRE | 0.9978 | 0.9979 |
| vg-giraffe | 0.659 | 0.837 |

**Founder-aware methods are coverage-robust (flat at both depths); vg-coverage is
depth-dependent** (0.66→0.84, still trailing at 50×). The low-coverage robustness is the core
value of founder-aware AF for pool-seq E&R (typically ~10×). (NB vg must run **rescue-off**
`--rescue-attempts 0` — default/`fast` giraffe is rescue-bound on these sims, >5 h→timeout;
rescue-off is AF-equivalent and ~35× faster: 8:49 at 10×.)

---

## 3. SV benchmark — kMate + vg, **hapFIRE = blank**, and WHY

The SV panel compares **kMate (direct)** and **vg (graph-native `vg pack`/`vg call`)**, with
**hapFIRE shown as a blank/greyed cell**. The "why" is load-bearing, so it's spelled out:

### 3.1 hapFIRE/HARP cannot estimate SVs from read evidence
From the HARP paper (Kessner et al. 2013, `papers/mst016.pdf`), HARP's entire signal is a
**per-base substitution likelihood**:

  P(read | haplotype) = ∏_i P(read_base[i] | haplotype_base[i], baseQ[i]),  alphabet {A,C,G,T}

— it scores reads by **base matches at SNP positions** the read covers. Consequences:

1. An SV (deletion/insertion) is a *length* change, not a base substitution. HARP has no
   term for the split-read / read-depth / breakpoint evidence SVs require, and no "absent"
   or "+Nbp" symbol a haplotype can carry at a position.
2. Reads spanning a deletion align with a gap/clip → no read base to compare at the SV → the
   SV contributes nothing (or noise) to the likelihood.
3. Its `vcf_processing` expects single-base `0|1` SNP alleles; symbolic/long SV alleles
   aren't parsed into a usable haplotype matrix.

So feeding SVs into hapFIRE's VCF does **not** yield an SV frequency from the data. hapFIRE
fundamentally estimates **founder (ecotype) frequencies** from SNP base-matches; its per-SNP
AF is itself just that founder-frequency vector projected onto the founders' SNP genotypes.

The blank cell is therefore the honest statement: **hapFIRE has no native SV estimation**
(its model is SNP-base only). The only route to SVs is an *external, LD-based aftermath
projection* — take hapFIRE's `_ecotype_frequency` vector and multiply by the founders' SV
genotypes — which is NOT hapFIRE measuring the SV and structurally **fails for SVs not in LD
with SNPs**.

### 3.2 Why we do NOT claim a kMate>hapFIRE SV advantage in Chapter 1 (global mode)
Critical subtlety: in **GLOBAL (chromosome-wide) mode**, kMate also estimates ONE h per
chrom (dominated by the millions of SNP-ish k-mers; SV k-mers are negligible to h) and
**projects it onto `var_pa` SV rows** — i.e. kMate-SV = (chrom-wide ĥ)·(SV genotype). That
is *the same operation* as the hapFIRE founder-projection. Since hapFIRE's chrom-wide ĥ is
excellent (its SNP R²=0.998 *is* that ĥ projected onto SNPs), **hapFIRE-projection would tie
kMate on SVs in global mode** — and in these founder-known sims the SV frequency *is* exactly
the founder projection by construction (g0 = whole founders), so there is no residual SV
signal to separate them. Claiming "kMate uniquely does SVs" from global-mode sims would be
refuted by anyone who runs the projection.

So the Chapter-1 SV panel claims only what's true:
- **kMate ≈ hapFIRE-projection ≫ vg-direct** (founder-aware projection beats founder-naive
  graph-coverage for SVs at low depth), and
- **hapFIRE has no *native* SV output** (blank cell, §3.1).

kMate-SV is accurate here (R²=**0.9987** on fully-called SVs; 54,968 SVs, |indel|≥50 bp;
refreshed 2026-08-25, was 0.993).

### 3.2a SV scoring mechanics + results (p80 n80_g0 cov10, `score_sv.py`)
Three SV gotchas, all load-bearing:

1. **Join is POSITIONAL, not by content key.** `(pos,ref_len,alt_len)` is NOT unique for
   SVs — two distinct ALT *sequences* of equal length can share a position (e.g. Chr1:5654574
   1→86 appears twice, AF 0.0125 and 0.95). Joining kMate↔truth on that key cross-joins and
   mismatches truth↔est (R² collapses to 0.58). The panel VCF, kMate output, and `truth_sv`
   are all the **same panel meta in identical row order**, so the join key is the panel row
   index over the SV subset (`svidx`); `score_sv.py` validates pos/len agreement after merge.
2. **vg-SV = coverage (AD), genotype-given-variants — the Ash in-silico-pool method**
   (`danielwood1992/Ash_pangenome` 8_2: `vg pack -Q5` → `vg call -k pack -v <known.vcf>` →
   AF = AD_alt/(AD_ref+AD_alt)). NOT de-novo `vg call`: de-novo discovers only ~3 k SVs at 10×
   and snaps AF to 0/1 (genotyper, 1–2 reads/snarl) → matched 12 k panel SVs at MAE 0.52,
   R²=−12.9 (rejected; kept only as `--vg-denovo` sanity). The `-v` panel SV VCF
   (`work/panel_sv_p80_Chr1.vcf.gz`, 54,968 records) is subset from the graph's own deconstruct
   (`pangenome_p80_chr1.vcf.gz`), so its `(chrom,pos,REF,ALT)` sequences map 1:1 onto `svidx`.
   ⚠️ vg #4443: `vg filter` corrupts AD on ~2.5% of sites — do NOT filter the GAM/pack before AD.
3. **Truth convention as in §2.1** — headline physical + fully-called, never MAR.

REFRESHED 2026-08-25 (kMate arm re-run; vg arm reused). Values as `allrec / fullcalled`:

| tool | basis | n (allrec / fullcalled) | MAE | R² | r |
|---|---|---|---|---|---|
| kMate | founder-aware projection | 54,968 / 45,786 | 0.0025 / **0.0023** | 0.998 / **0.9987** | 0.999 / 0.9994 |
| vg-giraffe | founder-naive AD coverage | 4,156 / 3,799 | 0.0714 / 0.0659 | 0.482 / 0.561 | 0.804 / 0.849 |
| hapFIRE | — (no native SV, §3.1) | blank | — | — | — |

The vg-SV cell (previously "pending job 35136778") is now filled: vg lands
**~29× worse MAE than kMate** on the fully-called basis, and recovers only **3,799 of 45,786**
fully-called panel SVs (8.3%) — so its headline weakness is *coverage of the SV catalog*, not
just noise on what it does call. Direction matches the SNP result (vg R²=0.66 at 10×).

Note the refresh also **removed the old all-SV vs fully-called gap**: pre-refresh kMate scored
0.756 (allrec) vs 0.992 (fullcalled), and that spread was attributed to the missing-call tail
(the physical-vs-MAR mechanism of §2.1). Post-refresh both bases sit at ~0.998, so most of that
gap was the pre-fix estimator's error on high-missingness records, **not** the truth convention.
§2.1's convention caveat still stands on its own evidence — but this particular spread is no
longer an example of it.

### 3.3 Where kMate's real SV advantage lives (Chapter 2)
kMate beats founder-projection only when **SV frequency ≠ founder-haplotype projection** —
i.e. when projection breaks: SVs not in LD with SNPs (~7% per our SV-SNP-tagging result),
recurrent/homoplastic SVs, panel SV-genotype errors/missingness, or recombination decoupling
SV from SNP haplotype. That requires **window/LD-block mode** (local k-mer evidence directly
on the SV) and/or the real heterogeneous panel — the **second chapter** of the benchmark.

---

## 4. Speed & CPU — alignment-free is the point

Speed is paper-critical, so it gets **clean, dedicated, sole-occupancy runs** (`--exclusive`,
matched 8 threads, `/usr/bin/time -v` → wall + user/sys CPU + %CPU + peak RSS), separate from
the (contended) accuracy runs. `benchmarks/speed_vs_hapfire/` (maintained; see its README).

**End-to-end fastq→AF for every tool** — this is what makes kMate's alignment-free design
visible:
- kMate: fastq → k-mer count → EM → AF (**no alignment, ever**).
- hapFIRE: fastq → minimap2 align → HARP → AF (alignment included).
- vg: fastq → giraffe map → AF (mapping included).

Including alignment is the fair choice: kMate pays for k-mer counting, vg pays for giraffe,
so hapFIRE must pay for its alignment too — otherwise it gets a "free" BAM kMate/vg don't.

### 4.1 Two numbers, stated precisely (the framing that survives review)

HARP runs **single-threaded** (`haplotype_generation.py:290` is a serial Python loop over
independent LD blocks, each a `harp like`/`harp freq` subprocess; `import multiprocessing`
is present but unused). So given a matched 8-core allocation, kMate uses ~3 cores, HARP ~1.

**CLEAN exclusive-node numbers (n80_g0, end-to-end fastq→AF, matched 8 threads):**
| | kMate wall | hapFIRE wall | WALL × | kMate CPU-s | hapFIRE CPU-s | CPU × |
|---|---|---|---|---|---|---|
| 10× | 76 s | 659 s | **8.7×** | 226 | 718 | 3.2× |
| 50× | 92 s | 799 s | **8.7×** | 431 | 1184 | 2.7× |

- **Wall-clock, matched 8-core node: ~8.7×** (consistent 10×/50×). Cross-validates the earlier
  independent `speed_vs_hapfire` bench (8.6× wall, 2.9× CPU). **Headline.**
- **CPU-time (user+sys core-seconds): ~3×.** Parallelism-invariant per-core floor (HARP's
  per-read×per-SNP likelihood vs kMate's k-mer-count+EM; both read every read once).
- ⚠️ A naive contended measurement gave hapFIRE 37 min → a spurious ~30×; the **clean
  exclusive-node run is the valid one (~8.7×)**. Do NOT quote 30× — it was node contention.

**Claim it as:** "On a matched 8-core node, kMate is **~8.7× faster wall-clock** than hapFIRE
*as distributed* (HARP is single-threaded); **~3× faster in CPU-time**." State the condition +
footnote the CPU-time so the "just parallelize it" referee is pre-answered ("8.7× as-shipped,
still ~3× per-core"). hapFIRE *could* be parallelized (independent block loop, `mp` imported)
→ that's why CPU-time is the parallelism-invariant floor.

NB the production cohort gap is *larger* than the per-sample 8.7× because the cohort adds
kMate-only scale-out wins: count-once (~1.7×) + node-local/RAM k-mer DB I/O fix (query
1944 s → 39 s). That stack — per-sample 8.7× × scale-out — is the weeks→afternoon experience.

### 4.2 vs vg — the cleaner alignment-free story
vg giraffe is genuinely multi-threaded (like kMate), so kMate-vs-vg wall-clock is apples-to-
apples — no single-threading caveat. giraffe mapping is the unavoidable expensive step, so the
alignment-free advantage holds here without qualification. (Note: default giraffe is
pathologically slow on these sims — a few hundred rescue-alignment reads at 10–40 s each; we
verify the `fast` preset is AF-equivalent and use a sane config so the vg speed reflects
normal giraffe, not the rescue pathology.)

---

## 5. Two chapters

- **Chapter 1 (this doc): chromosome-wide (global) h.** SNP parity (kMate≈hapFIRE) + SV
  (kMate≈hapFIRE-projection≫vg, hapFIRE blank) + speed (alignment-free).
- **Chapter 2: window / LD-block mode.** Per-block h on LD-coherent blocks (not fixed 10 kb).

  *Ghost-block problem (old clq0.9 map):* on median-7-variant LD blocks only ~42% get a local
  fit (≥50 observed k-mers); ~30% coverage-limited, ~28% zero-panel-k-mer. p80≈p231 → block
  thinness, NOT panel imbalance.

  *Fix = dynld_K500 map* (`analysis/grenenet_selection/blocks/results/blocks_mcf90/chr{N}_units_dynld_K500.tsv`):
  grow each CLQ0.9 block along LD until ≥500 panel k-mers. **Corrected/audited map (2026-06-19):
  22,939 units genome-wide / 6,123 Chr1, 72% k-mer-covered, median 32 var / 2.1 kb.** Local-fit
  rate → **~80.5%** (p231 80.7%, p80 80.3%; ~10% fallback, ~9% empty). Panel-independent.

  *Benchmark (cov10 n50_g3, dynld-window vs global, SNP R², via `score_ldblk.py`):*
  ⚠️ **BOTH COLUMNS ARE PRE-REFRESH** (window mode was out of scope for Phase 5, and these
  "global" numbers predate the corrected estimator). Kept for the record; re-run both arms
  before citing. The refreshed `--unit chrom` value for the p80 outcross row is **0.9694**.
  | regime | global | dynld-window |
  |---|---|---|
  | outcross (high recomb), p231 | 0.955 | **0.963** |
  | outcross, p80 | 0.966 | **0.970** |
  | selfing (low recomb), p231 dom500, 6 seeds | ~0.993 (s42) | 0.987 ± 0.004 |

  **Window BEATS global on the high-recomb outcross arm (both panels); global slightly beats
  window on the low-recomb selfing arm** (no recombination to recover → windowing adds noise).
  ⟹ use **window for recombinant/outcross pools, global for selfing/inbred** (why the evolved
  GrENE-net cohort runs global). The dynld map makes window mode viable + better-where-it-should-be.
  Figure: `benchmarks/ldblock_window_test/blockpanel_dynld_outcross_g3.png`.

  *Block-mode head-to-head vs hapFIRE* (hapFIRE is natively LD-block-based — the proper
  comparator). p80 outcross g3, **fair convention-free basis** (info>0.9 well-informed SNPs;
  the raw `var_pa`-MAR truth deflates hapFIRE to 0.899 by penalizing its physical estimates at
  high-missing sites — the SAME convention artifact as §2.1, removed here):
  | tool | MAE | R² | status |
  |---|---|---|---|
  | hapFIRE (LD-block HARP) | 0.0125 | **0.9872** | reused (unchanged) |
  | kMate-window (dynld) | 0.0215 | 0.970 | ⚠️ **STALE** — window mode was out of scope for the Phase-5 refresh |
  | kMate-`--unit chrom` | 0.0226 | 0.9694 | refreshed 2026-08-25 |

  **On the recombinant g3 arm, hapFIRE still clearly leads kMate** (0.9872 vs 0.9694) — and
  unlike g0, the refresh barely moved kMate here (0.966 → 0.9694). That contrast is the
  informative part: the per-founder/`Kf_w` fix removes a *normalization* error, which was
  dominating at g0, but it cannot fix the *structural* limitation at g3 — one chromosome-wide
  `h` cannot represent a recombinant mosaic. That is what window mode is for.

  **Revised honest ordering — it is regime-dependent, not global:**
  - **g0 (unrecombined founder mixture): kMate > hapFIRE** (0.9993 vs 0.9979, §2.2).
  - **g3 (recombinant): hapFIRE > kMate-chrom** (0.9872 vs 0.9694).

  So the blanket claim "SNP accuracy is NOT kMate's win" is **retired** — it was true only under
  the pre-refresh estimator. The defensible statement is: kMate wins SNP AF on unrecombined /
  selfing-like pools (which is the GrENE-Net production regime) and trails hapFIRE on genuinely
  recombinant pools until window mode is refreshed. kMate's other advantages are unchanged and
  regime-independent: **~8.7× wall-clock speed, alignment-free, and SVs (hapFIRE cannot do SVs
  at all)**.

  ⚠️ The `kMate-window (dynld) 0.970` row above predates the refresh, so the old
  "kMate-window > kMate-global confirms window mode helps on recombinant data" conclusion is
  **no longer supported by these two rows** — 0.970 (stale estimator) vs 0.9694 (refreshed) is
  not a like-for-like comparison. Re-run window mode before restating it.

  *(SV value-add of window mode — direct k-mer evidence on SVs not in LD with SNPs — is the
  next sub-test; see §3.3.)* Handoff: `analysis/grenenet_selection/BLOCKS_HANDOFF.md`.

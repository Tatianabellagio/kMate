# GrENE-net evolved cohort: why kMate runs in GLOBAL mode (not window/block)

**Decision (2026-07-01): estimate evolved allele frequencies with kMate in GLOBAL
mode (one founder-frequency vector `h` per chromosome per pool). Treat founder
haplotype linkage as preserved, and state that as an explicit modeling assumption.
Do NOT run window / coarse-block mode to "detect recombination."**

This closes a recurring question ("it's a selfing organism, shouldn't we give
window/block mode a chance to catch recombination?"). The short answer: selfing —
together with the ~3-generation timescale — is exactly *why* we can't and needn't.

---

## The question

kMate can run per-chromosome (GLOBAL) or per-window (window/block, optionally
`--local-only`). Window mode could in principle reveal new recombination that broke
founder haplotypes during the experiment. Is it worth running, or should we assume
everything stayed linked and run GLOBAL?

## Why the obvious diagnostic is circular

kMate projects `AF(variant v, pool s) = Σ_founders h[s,f] · G[f,v]` where `G` is the
**intact** founder genotype matrix. Two consequences:

1. **A chimeric (recombinant) haplotype is not in the model's vocabulary.** It cannot
   be represented as a mixture of intact founders.
2. **Any two variants with identical founder columns** (`G[·,A] ≡ G[·,T]`, i.e. in
   perfect founder LD) are **algebraically forced to equal projected AF for any `h`.**

So PC1-variance-explained / any linkage metric computed on kMate's *projected* evolved
AF reflects the **founder panel's** covariance, not the evolved population. It cannot
detect recombination. This makes the older
`analysis/grenenet_selection/archive/window_hapfreq_retired/gen9_window/breakage_global_vs_window.png` and the r≈0.99
`hapfreq/site4_block_vs_global_h.png` "everything is linked" reads **partly circular** —
consistent with the null, but not proof of it. (They were also built with the
prior-anchored + global-fallback window run, so the "window" side was partly global by
construction — a second, independent reason those panels agree.)

## The escape, and why it still fails here

Circularity is a property of putting both variants under **one** `h`. If A and T sit in
separate, independently-fit windows (`--local-only`), then `h₁ ≠ h₂` is allowed and a
breakpoint *between* the windows can show up as a spatially-coherent shift in local `h`.
The founder panel is not the fundamental barrier — window k-mers carry local haplotype
background. So recombination is, in principle, recoverable as an adjacent-window `|h₁−h₂|`
shift. The binding constraint becomes the **k-mer floor** (each window needs ≳1–2k
nonzero k-mers or `|h₁−h₂|` moves from estimation noise alone — see `FLOOR_DERIVATION`).

An adversarial audit (2026-07-01) found this test **not trustworthy to run and act on**,
for two compounding reasons:

### 1. Underpowered by construction
At ~97% selfing, a selfing event is a recombination no-op. Effective new cross-ancestry
crossovers ≈ 0.13 per lineage per generation. Over **~3 generations** (note: the "gen9"
label is a timepoint code, **not** the number of meioses) that is ≈ **0.4 events per
lineage**. After the detectable-frequency filter (a recombinant must reach ~10–15% of a
window's ancestry to move `h` above noise), the realistic yield is a handful of loci
genome-wide, plausibly **zero at most sites**. The most likely honest output is "no
detectable new recombination" — which GLOBAL already assumes.

### 2. The simulation noise floor does not transfer to real data
Closed-loop sims generate reads *from* the panel, so they miss the dangerous confounds —
the ones that are **both spatially coherent and replicate-reproducible**, and therefore
masquerade as recombination and survive both proposed filters:

- **Panel incompleteness (the killer).** A true evolving ecotype not in the 231-panel is
  approximated by a founder *mixture*; that best-fit mixture differs between adjacent
  windows (the missing haplotype resembles different founders in different regions) →
  a spatially-coherent local-`h` shift that is **pure model misspecification**,
  indistinguishable from a real breakpoint, and **deterministic → reproduces across every
  replicate**. Same signature as the previously-investigated spurious ecotype 9761.
- Reference/mapping bias, panel-assembly error, companion-species contamination — all
  spatially structured and deterministic.

And the **replicate-reproducibility filter runs backwards**: deterministic artifacts
reproduce (kept in); a drift-driven real recombinant arises in one lineage and need not
recur (thrown out). Reproducibility only survives under *parallel selection* — at which
point the test measures selection, not recombination.

Finally, the estimand `|h₁−h₂|` is the per-founder `h` difference that the h-certainty
analysis proved **no variance method (Fisher / bootstrap / CRLB) can calibrate** — so the
sim floor stands in for a quantity known to be non-identifiable at the `h` level.

## Decision

**Run GLOBAL mode.** It is both the practical choice (more k-mers per estimate → better
precision; local ≈ global where recombination is rare) and the honest one (no panel-based
pooled-AF method can separate the rare new recombination from panel-incompleteness
artifacts in this system). Intact linkage is an **assumption to state**, not a bug to
"fix" with a window run that would mostly manufacture false breakpoints.

### Methods sentence to use

> Evolved-pool allele frequencies were estimated with kMate in global mode, projecting
> observed k-mer counts onto the intact founder panel (one founder-frequency vector per
> chromosome per pool). Because the population is ~97% selfing and evolved for only ~3
> generations, new recombination is expected to be negligible and, where it occurs, to
> fall below the detection limit of panel-based pooled-sequence allele-frequency
> estimation; founder haplotype linkage is therefore treated as preserved.

## Update (2026-07-08): machinery archived; GLOBAL = `--unit chrom`

"GLOBAL mode" is now the unified estimator invoked as **`--unit chrom`** (byte-identical
to the old `--block-mode global`). The window/hapfreq *diagnostic* machinery referenced
throughout this doc — cross-chrom agreement, `site4_block_vs_global_h`, the `gen9_window/`
breakage script subtree — has been **archived** to `archive/window_hapfreq_retired/`
(the window stores `results/grenenet_kmate_window[_seedmix]` it read were deleted, and the
`--unit chrom` cohort emits no per-block `h_blocks`). Window-vs-chrom is settled in favor of
**chrom**: haploblocks (r²=0.1, eps=0) collapse to ~231 ≈ chromosome-wise, so the reframe is
"framing only." (The `gen9_window/` DATA directory under `results/` is kept.)

## Related
- `FLOOR_DERIVATION.md` / `benchmarks/localonly_p231/` — the k-mer floor (AF error ∝ nnz;
  ~1–2k k-mers to plateau; 5–10 kb windows sit below it).
- `WINDOW_UNIT_VALIDATION.md`, `gen9_window/` — the window-unit machinery (kept for GEA
  *test units*, not for AF estimation).
- Memory: `cohort-uses-global-mode`, `gea-blockwindow-benchmark-result`,
  `gea-blocks-redundant-ecotype-level`, `block-mode-kmer-floor`, `h-certainty-framework`,
  `grenenet-window-vs-global-decision`.

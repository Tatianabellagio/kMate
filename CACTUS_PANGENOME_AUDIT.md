# Cactus pangenome audit — findings

**Status as of 2026-04-27 14:25.** Three jobs pending; this doc updates as they complete.

## Bottom line

**The cactus pangenome is the right foundation for this project.** Earlier "10% syri-cactus match" was a metric artifact (position-match instead of interval containment). Real coverage is ~71–93% depending on metric, with the residual gap concentrated in cactus alignment dead zones (mostly pericentromeric).

## What the cactus pangenome contains

Built by `cactus-pangenome` (commit 16fae1c8) at `/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/`:

- **Source:** TAIR10 + 82 long-read 1001G+ collab assemblies, all chr-only (Chr1–Chr5, no organelles)
- **Built parameters:** `--reference TAIR10 --haplo --vcf --gfa --gbz --giraffe --maxLen 10000`
- **Outputs:** `raw.vcf` (6.58M records, all LV levels), `vcf.gz` (4.45M filtered), `gfa`, `gbz`, `giraffe` index
- **Run completed cleanly:** ~39h, "Workflow stopped. Success: True"

By size (Chr1 LV=0):
- 341,366 SNPs/MNPs
- 67,419 indels (1–50bp)
- 2,559 SVs (≥50bp)

## The 90% gap was a metric artifact

Initial finding: only 10% of syri SVs matched cactus records within 200bp position-tolerance. **This was wrong.** The right metric is *interval containment* (does the syri SV's position fall inside the reference span of any cactus bubble?).

Position-match misses cactus's giant multi-allelic bubbles that span kilobases. Example: Chr1:12640885 syri DEL is inside cactus bubble [12625993, 12646638] (refl 20,646bp, 79 alts). The bubble's POS is 14.9kb away from syri's call → position-match fails. But containment correctly catches it.

## Actual coverage numbers

**Distance from syri SV to nearest cactus record (Chr1):**
| Distance | Coverage |
|---|---|
| Within 100bp | **69.4%** |
| Within 1kb | 81.9% |
| Within 5kb | **92.9%** |
| >5kb (genuinely missed) | 7.1% |

For high-AC syri SVs (≥20 carriers): 78% within 100bp, only 4% >5kb away.

**Median distance: 4bp.** The two catalogs agree on positions to within a few bp for most SVs.

## What the residual gap is

**Cactus has alignment dead zones.** On Chr1, 277 gaps of >5kb between consecutive records, totaling 4.3Mb (14.1% of Chr1). Largest gaps are pericentromeric:
- Chr1:14.6–14.9Mb (256kb gap = the centromere)
- Chr1:3.3–3.4Mb (95kb)
- Chr1:14.2–14.3Mb (90kb)
- Chr1:16.0–16.1Mb (88kb)
- Chr1:12.3–12.4Mb (82kb)

**18.6% of syri Chr1 SVs fall in these gaps**, 10.5% of high-AC SVs. These are regions where minigraph couldn't chain alignments through (mostly heterochromatin / TE-rich).

These regions are difficult for short-read genotyping with *any* method — k-mer methods like freqk also fail there. So losing them isn't catastrophic.

## vcfbub on raw vs vcfbub on cactus-filtered

This is the actionable finding for the pipeline.

| Approach | LV=0 records | syri SV containment |
|---|---|---|
| `vcfbub -l 0 -r 100000` on **filtered** `vcf.gz` | 3,364,026 | **52.0%** |
| `vcfbub -l 0 -r 100000` on **raw** `raw.vcf` | 4,457,012 | **71.2%** |

Cactus's own pre-filtering step drops ~1.1M LV=0 records — many of which contain real syri SVs. **Use raw + vcfbub, not the pre-filtered output.**

Decomposition of the gap:
- Filtered → vcfbub → 52% containment
- Raw → vcfbub → 71.2% containment
- Improvement: 19 percentage points recovered by bypassing cactus's filter

By size, the vcfbub-on-raw containment is uniform at 68–76% across all bins:
- 0–50bp: 71.1%
- 51–200bp: 71.7%
- 201–500bp: 71.6%
- 501–1kbp: 68.0%
- 1–5kbp: 73.7%
- 5–10kbp: 67.9%
- >10kbp: 76.1%

## Pipeline implications

1. **Use `raw.vcf` + `vcfbub -l 0 -r 100000`** as PanGenie input, not the pre-filtered `vcf.gz`. Recovers 19% more syri SVs. (Job 56173 building this index now.)
2. **Cactus pangenome retains 100% of bubbles in PanGenie's k-mer indexing** (verified on Chr1, vs freqk's 12% on the syri panel). Major win.
3. **Heterochromatic / pericentric regions are inaccessible to both pipelines** — accept the ~20% loss; those would not produce reliable AF estimates from short reads anyway.
4. **Don't trust syri as ground truth.** Syri's pairwise calls without graph-consistency enforcement can include alignment artifacts; cactus's graph-aware approach is more conservative (potentially more reliable).

## Pending jobs (will update when they land)

- **56169** — cactus rebuild without `--maxLen 10000` cap (~40h ETA). Tests whether the cap was affecting bubbles ≤10kb cascadingly. Also recovers the 10/155 syri >10kb SVs we currently miss.
- **56172** — `vcfbub -l 0` on raw without `-r 100000` filter. Tests whether very large bubbles (>100kb) contain syri SVs in the dead zones.
- **56173** — PanGenie-index on the vcfbub-from-raw output (the better input). The new k-mer index for the production pipeline.

## What this means for the larger pipeline

The decision tree is now collapsed:

```
Use vcfbub -l 0 -r 100000 on cactus raw.vcf  ← THE INPUT
  ↓
PanGenie-index → unique k-mer set per multi-allelic bubble
  ↓
PanGenie kmer counts per pool-seq sample (or generalize for pool-seq)
  ↓
Custom pool-seq frequency model on top (Poisson emission + LD via blocks)
  ↓
Decomposition: bubble alt-allele frequency → per-SV/SNP frequency
```

The cactus + PanGenie path is sound. The earlier discussion about hybrid panels or PGGB rebuilds is no longer needed.

## What we learned about benchmarking

Two methodological lessons from this audit:

1. **Position-match against multi-allelic bubble starts is meaningless.** Always use interval containment for VCF-vs-VCF comparison. A bubble's POS is the snarl start; the variants encoded within span the entire reference window.
2. **Cactus has at least three filtering steps**: minigraph alignment (creates dead zones), `--maxLen` cap (drops large bubbles), and a post-deconstruct filter that drops more LV=0 records. Each loses different things. The raw VCF preserves more than the published `vcf.gz`.

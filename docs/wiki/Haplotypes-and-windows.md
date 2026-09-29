# Haplotypes and windows

Two ideas decide how you run kMate. Both come from **your species' biology**, not
from software defaults.

1. What a **founder haplotype** is — one column of your panel.
2. Over what **span** kMate estimates their mixture — a whole chromosome, or a window.

---

## What is a founder haplotype?

A kMate panel column is **one haplotype**: a single continuous sequence, one allele at
every variant. Not a diploid genotype.

That is forced by the method. To decide which k-mers a founder contributes, kMate
reconstructs that founder's *sequence* — and a sequence has one base per position. A
heterozygous genotype describes two sequences, so it cannot be one column.

Where your founder haplotypes come from depends on the organism:

| your founders | haplotypes per sample | how |
|---|---|---|
| **inbred lines** — *Arabidopsis* accessions, MAGIC/RIL founders, NAM parents | 1 | the two copies are near-identical, so the line *is* a haplotype (`--het missing`) |
| **outbred, phased individuals** — HPRC-style assemblies | 2 | split the sample into `sample.h1` and `sample.h2` (`--het split`) |

So 100 inbred accessions give 100 founder haplotypes; 100 phased outbred individuals
give 200. kMate estimates a frequency for each. See
[Building a panel](Building-a-panel#choosing---het).

---

## Over what span is the mixture estimated?

kMate estimates one founder-haplotype mixture **`h` per unit**, where a unit is a stretch
of the genome. This is the `--unit` setting, and it is the most consequential choice
you make.

The question it answers is: **over how long a stretch is ancestry constant?**

```
   one mixture per chromosome            one mixture per window
   ──────────────────────────            ──┬───┬───┬───┬───┬───
   ancestry constant along it              └───┴───┴───┴───┴──
                                          ancestry changes along it
```

### It depends on recombination

- **Little or no recombination** — a selfing species, inbred lines, or a founder (F0)
  mixture that has not yet been crossed. Each individual in the pool carries whole
  founder chromosomes. Ancestry does not change along the chromosome, so estimating one
  mixture from **all** of that chromosome's k-mers gives the best-determined answer.
  → `--unit chrom`

- **Real recombination** — an outcrossing species, or a few generations of crossing.
  Each chromosome is a **mosaic** of founder segments. A single chromosome-wide mixture
  would average that mosaic away, which is exactly the signal you wanted.
  → `--unit bp --window-bp 10000`

This is not a small effect. On a selfing benchmark, fitting per LD-block instead of per
chromosome moved the error from **0.0033 to 0.0080** and outliers from **0.001% to
0.625%** — because low-diversity regions such as centromeres cannot tell founder
haplotypes apart from local k-mers alone, so the local fit drifts. In the opposite
direction, a chromosome-wide fit on a genuinely recombinant pool reports one average
ancestry for a chromosome that does not have one.

### It also depends on sequencing depth

A window is only fittable if enough of its panel k-mers were **actually seen in the
reads**. Shrinking the window divides the same coverage into smaller pieces:

- a **big** unit has many observed k-mers → the mixture is well determined
- a **small** unit has few → the estimate gets noisy, and below a floor it is not
  attempted at all

`--min-kmers-per-block` (default **200**) is that floor, counted in **observed**
k-mers. Blocks under it get no local fit: they are left as missing, or fall back to the
chromosome-wide mixture.

So depth sets how fine you can go. At high coverage small windows are viable; at low
coverage the same windows are mostly floor-outs and you are better off with a larger
unit. **Check the fallback rate** — on a real 8-founder MAGIC run at 10 kb windows,
94.3% of windows were fitted locally and 5.7% fell back, and those were concentrated
15.8× at centromeres, where the panel simply has no distinguishing k-mers. That is
structural, not a depth problem — but a fallback rate that is high *everywhere* means
your windows are too small for your coverage.

---

## The options

| `--unit` | one mixture per | use when |
|---|---|---|
| `chrom` *(default)* | chromosome | selfing, inbred, or F0 founder pools |
| `bp` | fixed window (`--window-bp`) | recombinant pools |
| `ld` | LD block (`--ld-r2`) | recombinant pools, where you want blocks from the panel's own structure |
| `tsv` | your own blocks (`--blocks-tsv`) | you have a partition you trust |

`--unit ld` is **not** a safer default: on low-diversity regions it collapses, which is
what produced the 625× outlier rate above. Use it only on genuinely recombinant pools.

---

## Choosing, in practice

1. **Ask how much recombination separates your pool from the founders.** None or
   almost none → `chrom`. Several generations of outcrossing → `bp`.
2. **If windowed, start at `--window-bp 10000`** and look at the fallback rate.
   Widespread fallback means the windows are too small for your depth — widen them.
3. **Check the founder-haplotype mixture**, not just the frequencies. If it has
   collapsed onto a handful of haplotypes where you expected many, the unit is likely
   too small (or coverage too low) for a local fit to be identifiable.

Longer discussion of why local fits fail in low-diversity regions:
[`docs/EM_UNIT_CHOICE_AND_NONIDENTIFIABILITY.md`](https://github.com/Tatianabellagio/kMate/blob/master/docs/EM_UNIT_CHOICE_AND_NONIDENTIFIABILITY.md).

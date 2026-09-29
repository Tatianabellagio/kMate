# Haplotypes and windows

Two choices determine how kMate is run, both set by your species' biology:

1. What a **founder haplotype** is: one column of the panel.
2. Over what **span** their mixture is estimated: a chromosome, or a window.

---

## What is a founder haplotype?

A panel column is **one haplotype**: a continuous sequence with one allele at every
variant.

The method requires it. To determine which k-mers a founder contributes, kMate
reconstructs its *sequence*, which has one base per position. A heterozygous genotype
describes two sequences, so it cannot occupy one column.

How many haplotypes a sample yields depends on the organism:

| your founders | haplotypes per sample | how |
|---|---|---|
| **inbred lines**: *Arabidopsis* accessions, MAGIC/RIL founders, NAM parents | 1 | the two copies are near-identical, so the line *is* a haplotype (`--het missing`) |
| **outbred, phased individuals**: HPRC-style assemblies | 2 | split the sample into `sample.h1` and `sample.h2` (`--het split`) |

100 inbred accessions give 100 haplotypes; 100 phased outbred individuals give 200.
kMate estimates a frequency for each. See
[Building a panel](Building-a-panel#choosing---het).

---

## Over what span is the mixture estimated?

kMate estimates one founder-haplotype mixture **`h` per unit**, a unit being a stretch of
genome (`--unit`). The question it answers: **over how long a stretch is ancestry
constant?**

```
   one mixture per chromosome            one mixture per window
   ──────────────────────────            ──┬───┬───┬───┬───┬───
   ancestry constant along it              └───┴───┴───┴───┴──
                                          ancestry changes along it
```

### It depends on recombination

- **Little or no recombination**: selfing species, inbred lines, or an uncrossed founder
  (F0) mixture. Pool members carry whole founder chromosomes, so ancestry is constant
  along a chromosome and one mixture fitted from **all** its k-mers is best determined.
  → `--unit chrom`

- **Recombination**: outcrossing species, or several generations of crossing. Each
  chromosome is a **mosaic** of founder segments, which a chromosome-wide mixture would
  average away.
  → `--unit bp --window-bp 10000`

On a selfing benchmark, fitting per LD-block instead of per
chromosome moved AF error from **0.0033 to 0.0080** and outliers from **0.001% to
0.625%**: low-diversity regions such as centromeres cannot distinguish founder
haplotypes from local k-mers alone, so the local fit drifts. On a recombinant pool, a
chromosome-wide fit reports one average ancestry for a chromosome that has none.

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
unit. **Check the fallback rate**: on a real 8-founder MAGIC run at 10 kb windows,
94.3% of windows were fitted locally and 5.7% fell back, and those were concentrated
15.8× at centromeres, where the panel simply has no distinguishing k-mers. That is
structural, not a depth problem. A fallback rate that is high *everywhere* means
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
   Widespread fallback means the windows are too small for your depth; widen them.
3. **Check the founder-haplotype mixture.** Collapse onto a handful of haplotypes where
   many were expected means the unit is too small, or coverage too low, for the local fit
   to be identifiable.

Longer discussion of why local fits fail in low-diversity regions:
[`docs/EM_UNIT_CHOICE_AND_NONIDENTIFIABILITY.md`](https://github.com/Tatianabellagio/kMate/blob/master/docs/EM_UNIT_CHOICE_AND_NONIDENTIFIABILITY.md).

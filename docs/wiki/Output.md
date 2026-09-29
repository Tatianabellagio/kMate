# Output

`kmate run` writes a tab-separated table with a header, **one row per variant in the
panel**:

```
chrom  pos  ref_len  alt_len  alt_freq  info  n_called  se
```

| column | meaning |
|---|---|
| `chrom`, `pos` | position, 1-based |
| `ref_len`, `alt_len` | length of the REF and ALT alleles |
| `alt_freq` | **the estimate** — frequency of the ALT allele in the pool, 0–1 |
| `info` | how much founder mass had a genotype call here |
| `n_called` | how many founders had a genotype call here |
| `se` | standard error of `alt_freq` |

Alongside it, `*.h_per_chrom.npz` holds the estimated **founder mixture** — one value
per founder, summing to 1.

---

## There is no variant-type column

Work it out from the lengths:

| type | rule |
|---|---|
| SNP | `ref_len == 1 and alt_len == 1` |
| indel | `ref_len != alt_len` and `max(ref_len, alt_len) < 50` |
| structural variant | `max(ref_len, alt_len) >= 50` |

```python
import pandas as pd
df = pd.read_csv("MYPOOL_Chr1.tsv", sep="\t")
snp = (df.ref_len == 1) & (df.alt_len == 1)
sv  = df[["ref_len", "alt_len"]].max(axis=1) >= 50
```

---

## Comparing samples

Every sample run against the same panel produces **the same rows in the same order**.
So you can stack them directly:

```python
import pandas as pd
a = pd.read_csv("POOL_A_Chr1.tsv", sep="\t")
b = pd.read_csv("POOL_B_Chr1.tsv", sep="\t")
delta = b.alt_freq - a.alt_freq          # row-wise, no merge needed
```

> **Do not merge on `chrom` and `pos`.** A position can carry more than one variant
> (for example a SNP and an indel starting at the same base), so a position key will
> match the wrong allele. Join row-wise, or on all four of
> `chrom, pos, ref_len, alt_len`.

---

## Judging an estimate

- **`n_called`** is the honest support: a frequency derived from few founders is
  weaker, whatever `se` says. Filtering on it is usually wise.
- **`alt_freq` is a frequency among *called* founders.** Founders that were missing at
  a variant are excluded rather than counted as reference.
- **Check the founder mixture** before trusting per-variant numbers. If it has
  collapsed onto a handful of founders when you expected many, something upstream is
  wrong — coverage, the panel, or the wrong `--unit`.

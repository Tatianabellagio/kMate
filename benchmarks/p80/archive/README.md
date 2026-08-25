# p80/archive/

Superseded p80 simulation sets. Not consumed by anything — the live sims are in
`../sims/`.

## `sims_buggy_preFix_2026-06-18/`

Six regimes archived on the same date as the p231 set, under the same name — but the
diff (run 2026-08-25, since nothing recorded it at the time) shows **p80 was affected
differently, and its truth was never wrong.**

| artifact | archived vs live |
|---|---|
| `recomb_truth.tsv.gz` | **byte-identical uncompressed** (2,048,843 rows both; only the gzip container differs) |
| `ancestry.tsv`, `pool_weights.tsv`, `source_weights.tsv` | identical |
| `reads/r1.fq`, `r2.fq` | **differ** — 366,372,675 → 365,738,551 bytes (−0.17%) |
| `reads/` layout | 54 entries (per-clone subdirs) → 4 (merged FASTQs + BAM) |

So on p80 the regeneration was **read-side only**. This matters because the p231 fix
was the opposite — there the truth record set changed (462,947 monomorphic records
dropped by the 2026-06-02 segregating-only filter) while its atomized truth and pool
composition were untouched. See [`../../p231/archive/README.md`](../../p231/archive/README.md)
for that analysis.

The shared `sims_buggy_preFix_2026-06-18` name is therefore misleading: it labels one
event but two different changes, and on p80 nothing about the *truth* was buggy.

The read-side change is **not characterised** — file size alone does not identify it,
and the pool-composition inputs are byte-identical, so the founder draw is not the
cause.

### Consequence for old p80 numbers

Since p80's truth is unchanged, a pre-06-18 p80 number is scored against the same
yardstick as a post-06-18 one. Any difference comes from the reads, not the metric —
the opposite of the p231 situation, where the yardstick itself moved.

## `sims_raw/`

2026-05-22 `*_raw_p80_chr1` sims, superseded by the `hotspots_*` recombination-map
regimes the benchmark uses throughout.

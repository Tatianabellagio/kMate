# p231/archive/

Superseded p231 simulation sets. Not consumed by anything — the live sims are in
`../sims/`.

## `sims_buggy_preFix_2026-06-18/`

Six regimes (`cov10_n{50,231}_g{0,1,3}` ± `dom500`, `hotspots_p231_chr1`) set aside
2026-06-18 and regenerated.

**What actually changed — derived by diff on 2026-08-25**, since neither the git
history nor any doc recorded it. The dir name was the only provenance; this is the
re-derivation.

### 1. `recomb_truth_raw.tsv.gz` — 462,947 records dropped (17.7%)

Identical in all six regimes: **2,617,371 → 2,154,424 rows**. The retained count is
exactly the arch3 Chr1 `var_pa` record count (2,154,423 + header).

The dropped records are **monomorphic in the panel**:

| dropped rows | count |
|---|---|
| `truth_af == 0` | 420,941 |
| `truth_af == 1` | 1,573 |
| `0 < truth_af < 1` | **0** |
| NaN | 4 |

Not one polymorphic record among them, and the fixed table is a **strict subset**
(0 keys present in fixed but absent from buggy). That is the exact signature of the
**segregating-only panel filter applied 2026-06-02** — `bcftools view -e 'INFO/AC=0 ||
INFO/AC=INFO/AN'`, which dropped 17.8% of records genome-wide (`docs/PIPELINE_STATE.md`
§0). The truth tables had been built against the *pre-filter* panel and still carried
records `var_pa` no longer contained.

**So it was a record-set mismatch, not a value error.** On the 1,938,162 unique-key
rows shared by both, `truth_af` and `info` are identical — max |diff| = 0.

(Compare by *unique* keys only: `(chrom,pos,ref_len,alt_len)` is not unique for SVs, so
a naive merge cross-joins and manufactures ~281k spurious "differences". Same trap as
`../../SCORING_RULES.md` RULE 1.)

### 2. `recomb_truth_atomized.tsv.gz` — unchanged

Byte-identical uncompressed; only the gzip container differs. The atomization path was
not affected.

### 3. `reads/` — regenerated

`r1.fq`/`r2.fq` differ (366,389,200 → 364,053,706 bytes, −0.6%), and the per-clone
intermediate subdirs (54 entries) were dropped in favour of just the merged FASTQs +
BAM (4 entries). The read-side change is **not characterised** — size alone does not
identify it, and the pool composition inputs (`ancestry.tsv`, `pool_weights.tsv`,
`source_weights.tsv`, `region.bed`) are all byte-identical, so the founder draw is not
the cause.

## What this means for old numbers

Any p231 benchmark number computed before 2026-06-18 was scored against a truth table
carrying ~463k extra monomorphic records. Since those records are ~all `af=0`, including
them inflates any all-records metric — a tool that correctly predicts ~0 there gets a
large block of near-free agreement. **Pre-06-18 all-record R²/MAE are not comparable to
post-06-18 ones.** Fully-called / info-filtered subsets are far less affected, since the
dropped records are precisely the low-information tail.

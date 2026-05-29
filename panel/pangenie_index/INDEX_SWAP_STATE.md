# Custom k-mer index — state & validation plan

_Last updated: 2026-05-29_

## TL;DR

We have our own standalone k-mer-index builder
(`scripts/build_kmers_tsv.py`) that reproduces PanGenie-index's
`kmers.tsv.gz` without PanGenie's 9 GB graph. It has already produced an index
for **all 5 chromosomes** (`pang_135_haploid/ours_Chr*_kmers.tsv.gz`). The goal
is to replace the PanGenie-built index in the production cn_full pipeline with
our own. The format is a drop-in; what is **not yet validated** is that our
index yields the *same downstream results* as the PanGenie one. We are keeping
the PanGenie index in production for now and treating the swap + validation as a
follow-up task.

## Why the swap is sound in principle

The production PanGenie index and our index are built from the **same inputs**
(verified 2026-05-29):

| | Production PG index | Our index |
|---|---|---|
| builder | `panel/pangenie_genotyping/scripts/build_pangenie_index.sh` | `panel/pangenie_index/scripts/slurm_pang135_haploid.sh` |
| VCF | `pang_1001gplus_all.vcf.gz` (pang_135) | `pang_1001gplus_all.vcf.gz` (pang_135) |
| REF | `TAIR10.chr.fa` (non-iupacN) | `TAIR10.chr.fa` (non-iupacN) |
| k | 31 | 31 |
| GT handling | diploidize `X`→`X\|X`, then `PanGenie-index` | `--haploid` on raw vcfbub |
| output | `pang_135_pangenie_index_<chrom>_kmers.tsv.gz` | `ours_<chrom>_kmers.tsv.gz` |

**The index is built from pang_135, the genotype VCF is the 231-founder
v3qc_v3 panel — and that pairing is by design, not a mismatch.** The index only
supplies the *bubble windows + candidate unique k-mers*; `build_kmer_cn.py` uses
the v3qc_v3 VCF to fill in *which founders carry which allele* inside those
windows. pang_135 (78 cactus assemblies + the graph the 153 short-read ecotypes
were genotyped against) is a **superset**: PanGenie genotyping only assigns
existing alleles, it never creates new bubbles, so every k-mer the 231 panel
needs already exists in the pang_135 index.

## How cn_full is wired today (unchanged, still PanGenie index)

```
pang_135_pangenie_index_Chr1_kmers.tsv.gz   (PanGenie-built, ~97 MB)
        │  scripts/build_cn_full_v3qc_v3_chr1.sh  →  src/build_kmer_cn.py  (--treat-missing-as-n)
        │  VCF: founders_231_v3qc_v3.haploid.vcf.gz   REF: TAIR10.chr.iupacN.fa
        ▼
data/cn_full_231_v3qc_v3/cn_Chr1.cn.npz  (+ meta)   →  per_sample_per_chrom.py
```

The only line that changes for the swap is `KMERS=` in
`scripts/build_cn_full_v3qc_v3_chr1.sh` (line 17), then rebuild cn_full.

> Note on REF: the index was built with non-iupacN `TAIR10.chr.fa`; cn_full is
> currently built with `TAIR10.chr.iupacN.fa`. `build_kmer_cn.py` only uses REF
> to fetch flank sequence around bubbles, but keep this discrepancy in mind when
> diffing — ideally hold REF fixed across the PG-vs-ours comparison.

## What still needs validation (the 3 open doubts)

1. **haploid (ours) vs diploid (PG) equivalence.** Production PG index was
   diploidized (`X`→`X|X`) before indexing; ours uses `--haploid` on the raw
   vcfbub VCF. The slurm header asserts these *should* be byte-identical (same
   set of distinct path-tuples per bubble) but this is **not yet confirmed**.

2. **Unique-k-mer selection parity at full scale.** `kmer_index_construction.md`
   §7 reports 99.2% byte-identical bubbles on a *2-accession, 2 Mb* test only.
   Full-panel validation vs PanGenie's production index is **in progress**, and
   the per-allele cap selection *order* is not yet a strict PanGenie match
   (`select_unique_kmers` has a "we'll address strict ordering in validation
   phase" note). Consequence: our index may keep a slightly different *subset*
   of unique k-mers per bubble → cn_full density (and thus AF estimates) could
   shift. Both would be valid; they just need to be shown equivalent.

3. **`nan` empty-bubble sentinel.** `build_kmers_tsv.py` writes the literal
   string `nan` for bubbles with no surviving k-mers. `build_kmer_cn.py:87`
   (`kmers = p[3].split(",") if p[3] else []`) treats `"nan"` as truthy and
   would carry it as a bogus 3-char "k-mer" column. Harmless (never matches a
   real 31-mer, stays zero) but pollutes `kmer_index`/`bubble_id`. **Check
   whether PanGenie's own `kmers.tsv.gz` also emits `nan`** — if so this is
   pre-existing and cosmetic. One-line fix: have `parse_kmer_file` skip `"nan"`.

## Next step — validate that PanGenie and our index produce the same thing

Goal: show that swapping in our index does not change the production result.
Two levels, cheap → expensive:

### Level A — index-level diff (fast, no cn rebuild)
Compare the two `kmers.tsv.gz` for the **production panel**, bubble by bubble:
same bubble coordinates, same `unique_kmers` set (order-insensitive), same
overhang set. Quantify: % bubbles identical, and for differing bubbles whether
ours is a subset/superset of PG's k-mers.

> Subtlety: the existing `ours_Chr*` were built from `pang_1001gplus_all.vcf.gz`
> (pang_135), the **same** VCF PG indexed — so a *direct* PG-index-vs-ours diff
> is already possible today without rebuilding anything. That is the true
> apples-to-apples test of the two indexers. (Rebuilding ours from the v3qc_v3
> VCF would be a *different* comparison and is NOT needed for indexer parity.)

### Level B — downstream cn_full / AF diff (the result that matters)
1. Build a second cn_full from our index, holding VCF + REF identical to the
   current production build (only `KMERS=` changes):
   `cn_full_231_v3qc_v3_OURSIDX/cn_Chr1.cn.npz`.
2. Compare against `data/cn_full_231_v3qc_v3/cn_Chr1.cn.npz`: shapes, nnz,
   per-founder k-mer counts, and a column-wise (per-bubble) cn diff.
3. Run `per_sample_per_chrom.py` (global + window) on a fixed sample set (e.g.
   SEEDMIX reps) with each cn_full and diff the per-record AF TSVs. Pass
   criterion: AF differences within numerical noise / explainable by the
   k-mer-subset differences from Level A.

If Level A is clean, Level B should be a formality; if Level A shows
selection-order/subset drift, Level B quantifies whether it actually moves AF.

## Decision on record (2026-05-29)

Keep the PanGenie-built index in production for now. The superset argument makes
it safe; the swap to our own indexer is deferred to a future iteration, gated on
the validation above. Owner will pick up doubts 1–3 later.

## Pointers

- Builder: `panel/pangenie_index/scripts/build_kmers_tsv.py`
- Method writeup: `panel/pangenie_index/kmer_index_construction.md`
- Built indexes: `panel/pangenie_index/{pang_135_haploid,pang_135_diploid,pang_135_haploid_cap2x}/ours_Chr*_kmers.tsv.gz`
- PG reference index for comparison: `panel/pangenie_index/pg_reference/pgref_Chr*_kmers.tsv.gz`
- Production PG index in use: `panel/pangenie_genotyping/data/pang_135_pangenie_index_Chr*_kmers.tsv.gz`
- cn_full consumer: `src/build_kmer_cn.py`; wiring: `scripts/build_cn_full_v3qc_v3_chr1.sh`

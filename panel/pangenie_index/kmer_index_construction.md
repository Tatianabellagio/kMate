# k-mer index construction

This document describes the construction of the per-bubble unique-k-mer table
emitted by `panel/pangenie_index/scripts/build_kmers_tsv.py` (run at scale on the
pang_135 panel via `scripts/slurm_pang135_haploid.sh`). The algorithm reproduces
the `kmers.tsv.gz` output produced by `PanGenie-index` (PanGenie v4.2.1; Ebler
et al. 2022, Nat. Genet.), but as a standalone Python tool that avoids building
the 9 GB ancillary PanGenie graph state and accepts haploid cactus genotypes
natively.

> **Three distinct things, easy to conflate:**
> 1. **This in-house index** (`build_kmers_tsv.py` → `ours_<chrom>_kmers.tsv.gz`):
>    per-bubble windows + candidate unique k-mers, our reimplementation.
> 2. **The PanGenie genotyping index** (`panel/pangenie_genotyping/data/pang_135_pangenie_index_*`):
>    built by `PanGenie-index` to *genotype* the 153 short-read founders. It can
>    also serve as a `kmers.tsv.gz` source for `kmer_pa`, but is now kept only as
>    a **reference comparator** — production `kmer_pa` is built from the in-house
>    index (1); see §7.
> 3. **`kmer_pa`** (`src/build_kmer_pa.py` → founder × k-mer matrix): the
>    downstream product that consumes one `kmers.tsv.gz` index **plus** the panel
>    VCF. This doc is about (1)/(2); `kmer_pa` is built elsewhere.

## 1. Inputs

| Input | Description |
|---|---|
| `--vcf` | Pangenome VCF. Sites are expected to be vcfbub-filtered (top-level bubbles, no nested snarls). Both `.vcf` and `.vcf.gz`/`.bgz` are auto-detected. Diploid (`X|Y`) or haploid (`X`, with `--haploid`) genotypes accepted. Missing genotypes (`.`, `./.`, `.|.`) are tolerated and treated as N-allele paths (see §3.3). |
| `--ref` | Reference FASTA, indexed (`.fai`). Used to (a) validate REF, (b) extract bubble flanking sequence, (c) compute genome-wide k-mer counts via `jellyfish`. IUPAC ambiguity codes in REF positions are tolerated; ACGT-vs-ACGT mismatches are a hard error. |
| `-k`, `--kmer-size` | k-mer length (default 31, matching PanGenie default). |
| `--cap-biallelic`, `--cap-multiallelic` | Per-allele unique k-mer caps (default 16 / 32, matching PanGenie). |
| `--overhang-cap` | Per-side overhang k-mer cap (default 12, matching PanGenie). |
| `--no-add-reference` | Omit the synthetic all-REF path. PanGenie includes it by default. |
| `--jellyfish-threads`, `--jellyfish-hash` | `jellyfish count` thread count and hash size (defaults 4 / 100M; the production pang_135 run uses 8 / 3×10⁹). |
| `--keep-tempfiles` | Retain the intermediate `path_segments.fasta` + `.jf` hash for debugging (default: removed on exit). |

## 2. Outputs

Per chromosome: `<prefix>_<chrom>_kmers.tsv.gz`, a 5-column tab-separated table
with one row per bubble:

```
#chromosome   start   end   unique_kmers   unique_kmers_overhang
```

* `start`, `end`: 0-based half-open bubble coordinates on the reference.
* `unique_kmers`: comma-separated list of forward k-mer strings selected as
  uniquely identifying one allele within the bubble (see §4).
* `unique_kmers_overhang`: comma-separated list of forward k-mer strings drawn
  from the ±2k flanking regions, unique genome-wide (see §5).
* Both columns emit `nan` when no k-mers survive selection.

## 3. Stage 1 — VCF parsing and bubble construction

### 3.1 Per-site filters

For each VCF record, the following filters are applied (mirroring
`PanGenie-index/GraphBuilder`):

1. **Non-ACGT ALT.** If any ALT contains a non-ACGT base, the entire record is
   dropped (PanGenie prints a warning and skips). Variants with an IUPAC code
   in any ALT are therefore excluded from the index.
2. **Chromosome ends.** Records within `2·k` bp of either chromosome end are
   dropped, because flanking k-mers (which require `k-1` bases on either side
   of the bubble and an extra `k-1` for unique-genome-wide overhang queries)
   cannot be constructed safely.
3. **REF validation.** The VCF REF allele is compared base-by-base to the
   reference FASTA. A strict ACGT-vs-ACGT mismatch raises an error; positions
   where either side is non-ACGT (e.g., cactus normalizes IUPAC → `N` in REF
   while the unmodified TAIR10 reference carries the original IUPAC code) are
   tolerated, matching PanGenie's behavior of not validating REF at all.

### 3.2 Path representation

After filtering, each record's genotype field is parsed into an integer
allele-index per haplotype path. For diploid input, a sample contributes two
paths (`X|Y` → `[X, Y]`); for haploid input (with `--haploid`, used for native
cactus output), a sample contributes one path. Missing alleles (`.`) are
encoded as `-1`. If `--add-reference` is enabled (default), a synthetic
all-REF path (index 0 at every variant) is prepended to the path list.

### 3.3 Bubble clustering

Variants are scanned in (chrom, position) order. Consecutive variants whose
gap on the reference is less than `k-1` bp are merged into the same bubble;
otherwise a new bubble begins. A bubble is therefore a maximal cluster of
variants such that no k-mer of length k can span two adjacent bubbles
without crossing a flanking region.

### 3.4 Allele enumeration within a bubble

For each bubble of `m` variants, every haplotype path induces a tuple
`(a_1, ..., a_m)` where `a_i` is that path's allele index at variant `i`. The
**unique alleles** of the bubble are the unique such tuples, in order of first
appearance. Tuples containing any `-1` (missing) are excluded at this stage:
PanGenie materializes them as a synthetic N-allele that yields zero downstream
k-mers, so dropping them avoids redundant FASTA writes without affecting
k-mer output.

## 4. Stage 2 — Allele sequence construction and unique-k-mer selection

### 4.1 Allele sequence

For each unique allele `(a_1, ..., a_m)` of a bubble, the corresponding
nucleotide sequence is reconstructed as

```
left_flank + ALLELE(v_1, a_1) + INNER(v_1, v_2) + ALLELE(v_2, a_2) + ... + right_flank
```

where:
* `left_flank`, `right_flank` are the `k-1` reference bases immediately
  upstream / downstream of the bubble,
* `ALLELE(v_i, a_i)` is the REF sequence of variant `v_i` if `a_i == 0`,
  else the `a_i`-th ALT,
* `INNER(v_i, v_{i+1})` is the reference sequence strictly between
  `v_i.end` and `v_{i+1}.start`.

All unique allele sequences from all bubbles are concatenated into a single
`path_segments.fasta`, with reference unitigs between consecutive bubbles
preserved as their own FASTA entries so genome-wide k-mer counts on this
FASTA equal the panel-wide canonical-k-mer multiplicity.

### 4.2 Genome-wide canonical k-mer count

`jellyfish count -m k -C -s <hash_size>` is invoked on
`path_segments.fasta` to produce a canonical-k-mer multiplicity table. The
default hash size is 100M for small panels and is bumped to 3×10⁹ for
production runs on the full 135-assembly pangenome (matching PanGenie's
default). Throughout stage 3, the resulting `.jf` file is queried in-memory
via the official `dna_jellyfish` Python bindings; no `jellyfish dump` is ever
materialized.

### 4.3 Per-bubble unique k-mer selection

For each bubble, define the **local k-mer count** of a k-mer `m` as the number
of times `m` appears (forward strand) across all unique allele sequences of
the bubble, and the **global k-mer count** as the canonical multiplicity
returned by jellyfish on the path-segments FASTA.

A k-mer is **uniquely identifying** for allele `a` in bubble `B` if:
1. `m` is present in `a`'s sequence,
2. `local_count(m) == 1` (it occurs exactly once across all unique alleles of B),
3. `global_count(canonical(m)) == 1` (it occurs exactly once across the entire panel).

Within each bubble, k-mers are emitted by iterating alleles in ascending
allele-index order; for each allele, a round-robin queue selects up to
`cap_biallelic` (biallelic bubble) or `cap_multiallelic` (multi-allelic
bubble) k-mers per allele, ordered as PanGenie emits them (queue-pop order
from the per-allele std::map of unique forward k-mers).

A per-bubble total cap of `max(n_paths, 301)` k-mers is applied (matching
PanGenie). Bubbles with zero surviving k-mers emit `nan`.

## 5. Stage 3 — Overhang k-mer selection

For each bubble, a left and right **overhang region** of `2·k` bp is taken
immediately outside the bubble boundary, clipped to the previous/next
bubble's boundary on the same chromosome (so the same reference window is
never counted in two bubbles' overhangs). Within each side independently:

1. All forward k-mers of the region are enumerated, skipping any window
   containing a non-ACGT base.
2. K-mers are sorted lexicographically (PanGenie's `std::map` iteration
   order).
3. A k-mer is retained if its local count in the region is `1` and its
   global canonical count in the panel-wide jellyfish hash is also `1`.
4. Up to `overhang_cap` k-mers per side are emitted; the two sides are
   concatenated into the `unique_kmers_overhang` field.

Bubbles with zero surviving overhang k-mers emit `nan`.

## 6. Implementation notes

* **Language and dependencies.** Python 3.11; `pysam` for FASTA access;
  external `jellyfish` binary for k-mer counting and the `dna_jellyfish`
  bindings for in-memory queries. No PanGenie binary is required.
* **Coordinate system.** All internal coordinates are 0-based half-open;
  VCF positions are converted on input.
* **Canonical vs forward k-mers.** Allele-local k-mers and overhang k-mers
  are emitted as forward-strand sequences (matching PanGenie's
  `unique_kmers` storage); the genome-wide jellyfish hash is canonical, and
  the forward k-mer is canonicalized at query time only.
* **Memory.** The dominant memory consumers are (i) the jellyfish hash
  (3×10⁹ entries × ~8 B ≈ 24 GB on the 135-asm panel) and (ii) the in-memory
  list of all unique allele sequences across all bubbles. Peak observed RSS
  on the full pang_135 diploid run is ~18 GB.
* **Runtime optimization.** A C-language reimplementation of the per-k-mer
  inner loops in stage 3 (currently the dominant Python overhead) would
  close the remaining ~3× speed gap to native PanGenie; this is left as a
  future optimization and is not on the critical path for current use.

## 7. Validation

The implementation was validated against PanGenie-index v4.2.1 on a small
diploid test panel (2 accessions, Chr1:1–2 Mb, 3252 input variants → 1985
bubbles): 1969/1985 bubbles (99.2%) emit byte-identical
`unique_kmers` and `unique_kmers_overhang` fields. The remaining 16 bubbles
differ only in PanGenie occasionally including an overhang k-mer that
straddles the bubble-start position by one base (an off-by-one artefact of
PanGenie's left-flank window construction); no bubble emits a unique-k-mer
that is not also emitted by PanGenie.

> The small-panel figure above and the `build_kmers_tsv.py` module docstring
> report slightly different counts/mechanisms for the handful of differing
> bubbles (99.2% / 16 bubbles / overhang off-by-one here vs 99.4% / 12 bubbles /
> PanGenie phantom-`A` in the docstring). These predate the full-scale diff and
> have not been reconciled; treat the full-scale Level-A result below as the
> authoritative parity check.

**Full-scale validation: done (Level A, 2026-05-29).** The in-house index was
diffed against PanGenie's production index on the full 135-assembly pangenome
across **all five chromosomes** with `scripts/diff_index_vs_pg.py`. Because both
indexes are built from the **same inputs** (`pang_1001gplus_all.vcf.gz` +
`TAIR10.chr.fa`, k=31), any difference is attributable to the indexers alone;
the diff confirmed parity (ours is at worst a subset of PanGenie's k-mers per
bubble, never a wrong k-mer).

**Production now uses the in-house index.** After Level-A parity (above), the
swap is live: `scripts/build_kmer_pa_arch3.sh` defaults to the in-house
`ours_${CHR}_kmers.tsv.gz` (`INDEX=ours` → `data/kmer_pa_231_arch3_filt2inv`).
The pang_135 superset guarantees every k-mer the panel needs is present, so the
in-house index is a safe drop-in. The PanGenie-built index is kept only as a
**reference comparator** — `INDEX=pg sbatch scripts/build_kmer_pa_arch3.sh`
rebuilds the same matrix off the PG index into `data/kmer_pa_231_arch3_pgidx_filt2inv`.
The retired exploration (diploid byte-check, cap2x, the 2 Mb `pg_reference`
panel, the comparison test scripts, and the original `INDEX_SWAP_STATE.md`
working note) lives under `panel/pangenie_index/archive/`.

**Known issue (carried from the index swap, not yet fixed):** `build_kmers_tsv.py`
writes the literal `nan` for empty bubbles; `src/build_kmer_pa.py` parses a
`kmers.tsv.gz` field with `p[3].split(",") if p[3] else []`, which keeps `"nan"`
as a bogus 3-char "k-mer". Harmless (never matches a real 31-mer) but pollutes
the k-mer set; the one-line fix is to skip `"nan"` in the field parser. Fixing
it changes the `kmer_pa` build output, so it is gated behind a rebuild.

## 8. Default parameter table for paper methods

| Parameter | Value | Source |
|---|---|---|
| k-mer length | 31 | PanGenie default |
| Per-allele unique-k-mer cap (biallelic) | 16 | PanGenie default |
| Per-allele unique-k-mer cap (multi-allelic) | 32 | PanGenie default |
| Per-bubble total cap | `max(n_paths, 301)` | PanGenie default |
| Overhang window (each side) | 2·k bp | PanGenie default |
| Per-side overhang k-mer cap | 12 | PanGenie default |
| Jellyfish hash size | 3×10⁹ | PanGenie default |

## 9. References

* Ebler, J., Ebert, P., Clarke, W. E., Rausch, T., Audano, P. A., Houwaart,
  T., Mao, Y., Korbel, J. O., Eichler, E. E., Zody, M. C., Dilthey, A. T., &
  Marschall, T. (2022). Pangenome-based genome inference allows efficient
  and accurate genotyping across a wide spectrum of variant classes.
  *Nature Genetics*, 54, 518–525.
* Marçais, G., & Kingsford, C. (2011). A fast, lock-free approach for
  efficient parallel counting of occurrences of k-mers. *Bioinformatics*,
  27, 764–770. (`jellyfish`)
* Hickey, G., Heller, D., Monlong, J., Sibbesen, J. A., Sirén, J., Eizenga,
  J., Dawson, E. T., Garrison, E., Novak, A. M., & Paten, B. (2020).
  Genotyping structural variants in pangenome graphs using the vg toolkit.
  *Genome Biology*, 21, 35. (`vcfbub`)

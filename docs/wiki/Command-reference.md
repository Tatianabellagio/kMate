# Command reference

Every option, with its default. **Required** options have no default and must be given.

Run `kmate <command> --help` for the same information at the terminal.

> This page describes kMate 0.1.3. In 0.1.2, `run` still accepted `--kmer-weight` and
> `--emit-af-se` and spelled `--smooth-windows` as `--no-local-only`, and `decompose` had
> no `--gfa`. `kmate --help` is authoritative for the version you have installed.

| command | what it does |
|---|---|
| [`run`](#kmate-run) | estimate founder-haplotype frequencies and project to allele frequencies |
| [`build-index`](#kmate-build-index) | build the per-bubble k-mer index from a panel VCF |
| [`build-kmer-pa`](#kmate-build-kmer-pa) | build the founder-haplotype × k-mer matrix |
| [`build-var-pa`](#kmate-build-var-pa) | build the founder-haplotype × variant matrices |
| [`build-kmer-db`](#kmate-build-kmer-db) | count a pool's k-mers once for reuse |
| [`decompose`](#kmate-decompose) | multi-allelic → biallelic VCF |
| [`transfer-id`](#kmate-transfer-id) | copy `INFO/ID` onto a genotyped VCF |
| [`filter-pa`](#kmate-filter-pa) | re-filter an existing `kmer_pa` matrix |
| `selftest` | run the bundled fixture end to end |

---

## `kmate run`

### Inputs and outputs

| option | default | meaning |
|---|---|---|
| `--kmer-pa-prefix` | **required** | prefix of the k-mer matrices; `_<CHROM>.kmer_pa.npz` and `_<CHROM>.meta.npz` are appended |
| `--var-pa-prefix` | — | prefix of the variant matrices; `_<CHROM>.{var_pa,var_called,meta}.npz` appended. Use instead of the three options below |
| `--var-pa` | — | explicit path to `var_pa.npz` |
| `--var-called` | — | explicit path to `var_called.npz`. Guessed from `--var-pa` if omitted |
| `--var-meta` | — | explicit path to `meta.npz` |
| `--reads` | **required** | FASTQ(s) for one pool |
| `--kmer-db` | — | pre-counted Jellyfish DB from `build-kmer-db`; skips re-scanning the reads |
| `--sample` | **required** | sample name written into the output |
| `--out` | **required** | output TSV |
| `--chroms` | `Chr1 Chr2 Chr3 Chr4 Chr5` | chromosomes to process |
| `--h-only` | off | estimate the founder-haplotype mixture only, skip the allele-frequency projection |

### Estimation

| option | default | meaning |
|---|---|---|
| `--unit` | **`chrom`** | where the mixture is fitted: `chrom`, `bp`, `ld`, `tsv`. See [Haplotypes and windows](Haplotypes-and-windows) |
| `--normalize` | **`per_founder`** | M-step normalisation. `per_founder` divides each haplotype's update by its own k-mer content; `global` is the legacy multinomial form, under which k-mer-poor haplotypes collapse toward zero |
| `--haploblock-eps` | `0.0` | merge haplotypes whose k-mer presence differs by at most this fraction of the unit's k-mers. `0` merges only exactly identical ones |
| `--max-kmer-cov-mult` | `5.0` | repeat guard: drop k-mers seen more than this multiple of the estimated coverage (never below 1 read). `0` disables. See [Running kMate](Running-kMate#the-repeat-guard---max-kmer-cov-mult-on-by-default) |
| `--threads` | `4` | threads for k-mer counting |
| `--hash-size` | `3G` | Jellyfish hash size. Lower it on memory-capped jobs |

### Window options

Used when `--unit bp`, `ld` or `tsv`; ignored for `--unit chrom`.

| option | default | meaning |
|---|---|---|
| `--window-bp` | `10000` | window width for `--unit bp` |
| `--ld-r2` | `0.1` | r² cutoff for `--unit ld` |
| `--ld-window` | `100` | variants considered per LD block |
| `--ld-blocks` | — | precomputed LD blocks |
| `--blocks-tsv` | — | explicit block partition, required by `--unit tsv` |
| `--min-kmers-per-block` | `200` | a window with fewer *observed* k-mers gets no local fit |
| `--smooth-windows` | **off** | tie neighbouring windows together: enables the anchor and smoothing below, and fills low-support windows with the chromosome-wide mixture instead of NaN. Off by default, so each window uses only the information its own k-mers carry |
| `--global-anchor-weight` | `0.3` | anchor strength toward the chromosome-wide mixture. Applies only with `--smooth-windows` |
| `--hmm-smooth-passes` | `5` | cross-window smoothing passes. Applies only with `--smooth-windows` |
| `--hmm-smooth-alpha` | `0.5` | smoothing strength; smaller smooths more |
| `--hmm-smooth-recomb-rate` | `4e-08` | recombination rate used by the smoother |

### Deprecated

| option | meaning |
|---|---|
| `--block-mode` | alias for `--unit`: `global`→`chrom`, `window`→`bp`. `--unit` wins if both are given |
| `--local-only` / `--no-local-only` | old spelling; `--no-local-only` is `--smooth-windows` |

---

## `kmate build-index`

| option | default | meaning |
|---|---|---|
| `--vcf` | **required** | the **graph** VCF (multi-allelic, before `decompose`) |
| `--ref` | **required** | reference FASTA, indexed |
| `--out` | **required** | output prefix; writes `<prefix>_<CHROM>_kmers.tsv.gz` |
| `-k`, `--kmer-size` | `31` | k-mer length |
| `--haploid` | off | accept haploid genotypes. **Required for a haploid panel VCF** |
| `--cap-biallelic` | `16` | k-mers kept per allele at biallelic bubbles |
| `--cap-multiallelic` | `32` | k-mers kept per allele at multi-allelic bubbles |
| `--overhang-cap` | `12` | k-mers kept per bubble side |
| `--no-caps` | off | disable the caps above |
| `--no-add-reference` | off | omit the synthetic all-reference path |
| `--emit-overhang` | off | also write overhang k-mers |
| `--jellyfish-threads` | `4` | threads for the reference k-mer count |
| `--jellyfish-hash` | `100000000` | hash size for it |
| `--keep-tempfiles` | off | keep intermediates for debugging |
| `--tmp-dir` | next to `--out` | where the temporary files go (~2× genome size per assembly set; ~10 GB for 135 *Arabidopsis* assemblies). Use **node-local disk**: kMate warns when it is a network filesystem |

---

## `kmate build-kmer-pa`

| option | default | meaning |
|---|---|---|
| `--kmers` | **required** | `kmers.tsv.gz` from `build-index` |
| `--vcf` | **required** | the panel VCF (decomposed, biallelic) |
| `--ref` | **required** | reference FASTA. IUPAC codes are normalised to N internally |
| `--chrom` | **required** | chromosome to build |
| `--out` | **required** | output prefix |
| `--treat-missing-as-n` | off | mask a missing genotype with N over the record, dropping the overlapping k-mers. Without it, missing is treated as reference |
| `--filter-production` | off | apply the column filter below at build time |
| `--min-ac` | `2` | keep k-mers carried by at least this many haplotypes. **Use `1` on small panels**: on 8 founders, `2` discarded ~46% of k-mers |
| `--invariant-margin` | `1` | drop k-mers carried by more than `F - margin` haplotypes |
| `--max-bubbles` | — | stop after this many bubbles (testing) |

---

## `kmate build-var-pa`

| option | default | meaning |
|---|---|---|
| `--vcf` | **required** | the same panel VCF. Must be biallelic and haploid |
| `--out` | **required** | output prefix; writes `.var_pa.npz`, `.var_called.npz`, `.meta.npz` |
| `--chrom` | — | restrict to one chromosome |

---

## `kmate build-kmer-db`

| option | default | meaning |
|---|---|---|
| `--reads` | **required** | FASTQ(s) for one pool |
| `--out` | **required** | output Jellyfish DB |
| `--threads` | `4` | counting threads |
| `--hash-size` | `3G` | Jellyfish hash size |

---

## `kmate decompose`

Runs two bundled PanGenie scripts; see [Building a panel](Building-a-panel#acknowledgement).
Give **either** `--gfa` **or** both `--annotated-vcf` and `--biallelic-catalog`.

| option | default | meaning |
|---|---|---|
| `--genotyped-vcf` | **required** | the VCF to decompose; with `--gfa`, the Minigraph-Cactus VCF of that graph |
| `--gfa` | — | the graph's gzipped GFA; annotates `--genotyped-vcf` against it |
| `--annotated-vcf` | — | annotated multi-allelic VCF carrying `INFO/ID` (instead of `--gfa`) |
| `--biallelic-catalog` | — | biallelic catalog with the same symbolic IDs (instead of `--gfa`) |
| `--out` | **required** | output VCF, sorted; indexed when it ends in `.gz` |
| `--annotate-vcf-script` | bundled | use another copy of `annotate_vcf.py` |
| `--convert-to-biallelic` | bundled | use another copy of `convert-to-biallelic.py` |
| `--haploidize` | off | also reduce genotypes to one allele |
| `--het` | **`missing`** | with `--haploidize`: `missing` for inbred founders, `split` for outbred phased ones |
| `--threads` | `4` | threads for `bcftools` |
| `--keep-temp` | off | keep intermediates |
| `--citation` | — | print attribution and exit |

---

## `kmate transfer-id`

| option | default | meaning |
|---|---|---|
| `--cactus` | **required** | annotated VCF carrying `INFO/ID` |
| `--pg` | **required** | genotyped VCF from the same graph, lacking the IDs |
| `--out` | **required** | output VCF |

---

## `kmate filter-pa`

Applies the `build-kmer-pa` column filter to a matrix that was built unfiltered.

| option | default | meaning |
|---|---|---|
| `--in-prefix` | **required** | existing `kmer_pa` prefix |
| `--out-prefix` | **required** | filtered output prefix |
| `--min-ac` | `2` | as in `build-kmer-pa` |
| `--invariant-margin` | `1` | as in `build-kmer-pa` |

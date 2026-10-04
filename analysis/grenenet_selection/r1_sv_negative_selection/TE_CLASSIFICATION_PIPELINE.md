# TE / repeat classification of the arch3 SVs

How every structural variant (≥ 50 bp) in the arch3 panel got a repeat / TE class: the steps that were
run, what each one does, what went wrong on the way and how it was handled. Built 2026-10-03/04.

**Outputs** (`results/sv_adaptive/graffite/`)

| file | rows | one row per |
|---|---|---|
| `sv_repeat_class_records.csv` | 227,851 | SV record (panel allele), insertions and deletions |
| `sv_repeat_class_events.csv` | 89,940 | SV event (Truvari collapse), keyed to `sv_event_catalog.csv` |

---

## 0 · Input: the arch3 panel

`panel/arch3/chr{N}/merged_231_chr{N}_final.vcf.gz`: 231 founders, biallelic, 8.49 M records. The
biallelic records come from the HPRC Minigraph-Cactus decomposition (`prepare-vcf-MC/annotate_vcf.py`
+ `convert-to-biallelic.py`), **not** from `vcfwave`. The HPRC pipeline says outright that its input must
not be vcfwave-decomposed, and it never calls vcfwave. (`panel/arch3/README.md` line 42 says the
opposite; that line is wrong.)

`annotate_vcf.py` writes a type into every record ID (`Chr5-2-INS-…`), decided by allele length
(`define_id`): `INS` if REF is 1 bp, `DEL` if ALT is 1 bp, otherwise **`COMPLEX`**. `COMPLEX` means the
ALT path runs through *different* graph nodes from REF, so it is a replacement, not a pure gain or loss.

## 1 · Sites-only SV VCF for GraffiTE

`_make_graffite_vcf.py` (via `_submit_graffite_vcf.sh`, job 39593148)
→ `results/sv_adaptive/graffite/sv_sites.vcf`

- keeps biallelic, sequence-resolved records with |len(ALT) − len(REF)| ≥ 50: **227,851 records**
  (173,010 insertions, 54,841 deletions);
- drops genotypes (`--genotype false` doesn't need them);
- rewrites IDs as `{chrom}_{rec}`, where `rec` is the record's position in its chromosome VCF. This is the
  same key as `insertion_context.csv`, `rm_cargo_allins.csv` and `arch3_event_map.npz` (as `Chr|rec`).
  Native IDs are full bubble paths of thousands of characters; GraffiTE caps IDs at 50.

## 2 · GraffiTE

`_submit_graffite.sh`, job 39593195, 2 h 58 m, 448/448 tasks.

- GraffiTE `f37674b` (upstream HEAD, 2026-01-15), `--vcf sv_sites.vcf --genotype false`, standard
  profile, container `cgroza/collection/graffite:latest` (pre-seeded; byte-identical to
  `telib/graffite.sif`).
- Library `telib/combined_tair10_repbase.rm.lib`: TAIR10 TE families (320) + RepBase 31.09 `athrep.ref`
  (1,136), with headers converted to `>name#Class/Superfamily`.
- Its masking step (`bin/repmask_vcf.sh`), per chromosome:
  1. builds the query with `bcftools view --types indels` (ALT for ILEN > 0, REF for ILEN < 0);
  2. runs `RepeatMasker -lib … -s` (slow / sensitive), with low-complexity masking **on**;
  3. runs `build_dictionary.pl` and OneCodeToFindThemAll, which stitch fragments into copies;
  4. computes `total_match_span` as the union of OneCode copies over the allele length;
  5. keeps `total_match_span > 0.80` → `pangenome.vcf` (55,430 records), then runs the TSD search.

### What GraffiTE did not do, and why

**It skipped 71,865 of the 227,851 records (31.5%).** bcftools types a record as an indel only if,
after trimming shared bases, the shorter allele is empty. Everything else becomes `OTHER`, which
`--types indels` leaves out, with no warning:

| group | records | shorter allele vs longer |
|---|---|---|
| `INS`/`DEL`, one reference base also swapped (`C` → `TAAAGCG…TTTT`) | 17,682 | 1 bp unshared |
| `COMPLEX`, ≤ 5 bp unshared | 11,293 | near-clean gain/loss plus a SNP or small indel at the breakpoint |
| `COMPLEX`, 6–50 bp unshared | 6,157 | |
| `COMPLEX`, > 50 bp unshared | 36,733 | real replacement: both alleles carry their own sequence |

Of the 54,231 `COMPLEX` records GraffiTE masked only 48. The skipped records are biased towards long
SVs: median 1,355 bp against 617 bp for those masked, and about half of all insertions > 5 kb.

This is by design, not a fault in our input. GraffiTE's own discovery routes (svim-asm `--types INS,DEL`;
Sniffles2 filtered to `SVTYPE == INS|DEL`) only ever produce clean insertions and deletions. On
[issue #20](https://github.com/cgroza/GraffiTE/issues/20) its authors say they ran HPRC from a
decomposed Minigraph-Cactus VCF, and that complex SVs are out of scope: *"designed … for de novo
transposition events (which are not complex SVs)"*.

## 3 · The skipped records, masked with GraffiTE's own steps

`_graffite_prep_skipped.py`
- copies GraffiTE's raw per-chromosome output out of the Nextflow work dir (`telib/nxf_work`) to
  `graffite/raw/{chrom}/` (query FASTA, raw RepeatMasker `.out`, `.onecode.out`, `.tbl`, OneCode log,
  unfiltered annotated VCF);
- writes the 71,865 skipped records to `graffite/skipped/{chrom}/indels.fa`, in GraffiTE's convention:
  the whole ALT for insertions, the whole REF for deletions;
- checks per chromosome that masked + skipped = input.

`_submit_graffite_skipped.sh` (array 1–5; job 39624270, 56 min – 1 h 27 m per chromosome)
- runs GraffiTE's `repmask_vcf.sh` steps 2–3 verbatim: same container, library, `RepeatMasker -s`,
  `build_dictionary.pl`, OneCode, and the same parse of its output. It does not run the span filter or
  the TSD search;
- **parity check:** run on 300 sequences GraffiTE had itself masked (Chr5), the output matched
  GraffiTE's row for row: RepeatMasker 580/580 rows, OneCode 214/214 rows.

Operational notes:
- The first submission (39624242) failed in under a second on **n0041.savio4**, whose `/global/scratch`
  mount is broken (`transport endpoint shutdown`; it also caused the `trim` task failures on 2026-10-02).
  Resubmitted with `--exclude=n0041.savio4`.
- RepeatMasker's BLAST cache is redirected to scratch (`--home telib/rmhome_*`), because `$HOME` locking
  breaks `makeblastdb` under apptainer.
- The script `cd`s to an absolute path, not `$SLURM_SUBMIT_DIR`, which inside an `salloc` shell points at
  wherever salloc was started.

## 4 · Classification

`_graffite_classify.py` (env `basic`) → the two output tables.

All 227,851 records now have the same masking protocol. Everything is computed from the **raw
RepeatMasker `.out`**, not from GraffiTE's OneCode-based span (§5), as the union of hit intervals over
the classified allele (`seq_len`):

| column | definition |
|---|---|
| `span_te` | union of TE-class hits / allele length. TE = every class except `Simple_repeat`, `Low_complexity`, `Satellite` |
| `span_simple` | union of `Simple_repeat` + `Low_complexity` + `Satellite` hits / allele length |
| `span_any` | union of all hits / allele length |
| `span_graffite` | GraffiTE's own `total_match_span` (union of OneCode copies); for reference only, Helitron-blind |
| `cargo` | on `span_te`: `te_full` ≥ 0.80 · `te_most` ≥ 0.50 · `te_partial` ≥ 0.20 · `te_trace` > 0 · `non_te` = 0 (same cuts as `_rm_cargo_classify.py`; 0.80 is GraffiTE's own pME threshold) |
| `dominant` | what covers most of the allele: `te`, `simple`, or `unmasked` |
| `top_family`, `top_superfamily`, `top_family_frac`, `top_div` | the TE family covering the most of the allele (union bp per family), its library class, its coverage, and the lowest divergence among its hits |
| `helitron_bp`, `n_te_families`, `onecode_copies` | |
| `hprc_type` | `INS` / `DEL` / `COMPLEX`, read from the panel record ID |
| `unshared_bp` | bases of the shorter allele not covered by the prefix/suffix it shares with the longer one; > 0 means a replacement |
| `annotated_by` | `graffite` (155,986) or `graffite_steps_rerun` (71,865) |
| `event_key`, `is_rep` | event from `arch3_event_map.npz` (Truvari collapse defaults); representative = first record by position (`--keep first`) |

Event table: the representative record's class, plus `n_records`, `n_records_te_full`,
`frac_records_same_cargo`, `frac_records_same_superfamily` and `any_replacement`. The 89,940
representatives are exactly the keys of `sv_event_catalog.csv`. The 1,383 records of exactly 50 bp are
outside the event map, which uses > 50 bp, and have no event.

### Counts

Records:

| | non_te | te_full | te_most | te_partial | te_trace |
|---|---|---|---|---|---|
| insertions, masked by GraffiTE | 38,847 | 59,357 | 9,768 | 9,558 | 5,663 |
| insertions, skipped records | 12,020 | 17,923 | 7,959 | 6,924 | 4,991 |
| deletions, masked by GraffiTE | 18,715 | 10,308 | 1,522 | 1,354 | 894 |
| deletions, skipped records | 5,557 | 7,308 | 3,403 | 3,641 | 2,139 |

Events:

| | non_te | te_full | te_most | te_partial | te_trace | total |
|---|---|---|---|---|---|---|
| insertions | 17,213 | 21,750 | 4,774 | 4,036 | 2,458 | 50,231 |
| deletions | 21,638 | 12,295 | 2,235 | 2,118 | 1,423 | 39,709 |

Most common superfamilies among `te_full` events:
- insertions: RC/Helitron 4,886 · LTR/Copia 3,849 · DNA/MULE-MuDR 2,663 · LTR/Gypsy 2,554 · LINE/L1 2,323
- deletions: LTR/Gypsy 3,579 · RC/Helitron 1,885 · DNA/MULE-MuDR 1,617

Only 871 of the 25,532 multi-record events have records in different cargo classes.

## 5 · Problems found in GraffiTE (upstream, not ours)

**OneCode drops Helitrons.** OneCodeToFindThemAll emits no RC/Helitron copies at all. Its script
relabels `RC/Helitron` as `DNA/RC` (lines 303–304), but none appear in its output. On the
GraffiTE-masked records there are 41,555 raw Helitron hits (the most common class) and 0 in
`.onecode.out`. So GraffiTE's `total_match_span` and its `pangenome.vcf` are Helitron-blind: of the
21,703 records that are entire Helitrons (`te_full`, top superfamily RC/Helitron), **6** are in
`pangenome.vcf`. This is why the classification works from the raw `.out`. The cause hasn't been
traced.

**The TSD search runs without reference flanks.** Upstream commit `3f062c1` (2025-03-26, "Replace
bedtools getfasta with samtools faidx") changed `prepTSD.sh` so the flank FASTA is named by region
(`Chr2:33047-33077`). `TSD_Match_v2.sh` still looks flanks up as `{id}__L` / `{id}__R`, which never
match. Every TSD search since that commit compares only the unmasked edges of the SV allele itself,
with no flanking reference. Here: 395 PASS of 42,837 tested, and 31,771 of the failures are tagged
`poly-A`/`poly-T`, which is really the label for an empty comparison window. No upstream issue reports
it on `main` (checked 2026-10-04).

**Fixed and re-run (§5b).**

## 5b · TSDs, re-run with the flank fix

`graffite_prepTSD_flank_names.patch`, applied to the local GraffiTE checkout (`telib/graffite_src`),
restores the pre-`3f062c1` flank step: BED windows + `bedtools getfasta -nameOnly`, so flanks are named
`{id}__L` / `{id}__R` again. `_submit_graffite_tsd.sh` (job 39626019, 25 min) re-runs only the TSD stage,
entering at GraffiTE's documented `--RM_dir` point with the published `run/2_Repeat_Filtering/{1..5}`, so
nothing is re-masked. Output: `graffite/run_tsd/3_TSD_search/`. The script refuses to run if the patch
is not applied.

- Test on 200 Chr2 records before launching: PASS 3 → 114.
- Full run: **17,124 PASS of 42,837 tested** (40%), against 395 before; 15,549 insertions, 1,575
  deletions.
- TSD lengths match the superfamilies broadly: DNA/MULE-MuDR peaks at 8–9 bp, DNA/hAT at 7–8,
  LTR/Copia and LTR/Gypsy at 4–5; LINE/L1 is spread out, as expected.

Columns: `tsd_status` (`PASS` / `FAIL` / `not_tested`), `tsd_seq`, `tsd_len` per record; per event,
`n_records_tsd_tested`, `n_records_tsd_pass`, `rep_tsd_status`, `rep_tsd_len`.

Limits:
- **Scope is GraffiTE's own**: single-copy (`n_hits=1`) records whose OneCode span is > 0.80. Everything
  else is `not_tested` (185,014 records): replacements, partial TEs, multi-copy records, and Helitrons
  (OneCode drops them, and they make no TSD anyway). `not_tested` is not evidence of absence, and neither
  is `FAIL` for an old element, whose TSDs erode.
- **Short calls are suspect.** 27% of PASS calls are exactly 4 bp, the shortest the BLAST step
  (`word_size 4`) can report; PIF-Harbinger (expected 3 bp) is 63% 4 bp. Some lengths also look a base
  short (MuDR 8 rather than 9). Upstream has fixes for both, **merged only into the `v1.1dev` branch**:
  #117 (junction offsets taken from `WIN` instead of `WIN + 1`) and #119 (a TSD under 6 bp counts only
  if it sits on the junction). `v1.1dev` is 439 commits ahead of `main` with a rewritten TSD module, and
  has not been adopted here. Until it is, treat 4–5 bp PASS calls as weak evidence.

## 6 · Differences from the earlier cargo table (`rm_cargo_allins.csv`)

Same library and container. Different settings:

| | `rm_cargo_allins.csv` | this classification |
|---|---|---|
| records | insertions > 50 bp (172,220) | insertions and deletions ≥ 50 bp (227,851) |
| RepeatMasker | default speed, `-nolow -norna -no_is` | `-s`, low-complexity and RNA masking on (GraffiTE's setting) |
| TE coverage | all hits, satellites included | TE classes only; simple / low-complexity / satellite kept separate |
| denominator | net insertion size | length of the classified allele |

On the 122,516 insertions both cover (GraffiTE-masked set): cargo class agrees for **95.4%**,
r(`span_te`, old `te_frac`) = 0.934, and where both call `te_full` the top family agrees 99.6%. The
largest disagreement is **2,390 insertions that were `te_full` and are now `non_te`**. They are short
(median 120 bp), entirely simple repeat or satellite (`span_simple` median 1.0), and under `-nolow` were
claimed by TE consensus sequences (ATCOPIA68, Gypsy-26, …) or counted satellites (ATCLUST1, SAT-5S,
ATMSAT1) as TE. GraffiTE's README warns against `-nolow` for this reason.

## 7 · Caveats

- **Library-limited.** `non_te` / `unmasked` means "not in TAIR10 + RepBase", not "not a TE". Families
  absent from both libraries are invisible; a pangenome-derived library (panEDTA on our assemblies) is the
  way to close that gap.
- **Replacement records** (`COMPLEX`, or `unshared_bp > 0`; flagged per event by `any_replacement`): the
  class describes the whole longer allele, including sequence that replaces reference. Above 300 bp they
  are much less often whole TEs than clean records (insertions 300 bp–1 kb: 38% vs 53% `te_full`;
  > 5 kb: 37% vs 68%). Analyse them as their own group, or separately.
- **Superfamily labels are the library's.** A few TAIR10 Helitron entries carry plain `DNA`, including
  the common `ATREP18_AT1TE50795`, so Helitron counts are slightly low.
- **Detection power falls with length.** Below ~100 bp an alignment search finds little, so TE content
  and size are confounded for both biological and technical reasons. Stratify on size.
- **TSDs cover a minority of records** (§5b): 42,837 of 227,851 tested, short calls are weak, and
  everything untested has no polarity: "insertion" and "deletion" are relative to Col-0, not ancestral
  state.
- **Units.** Records are panel alleles. Per the project rule, SV analyses should use the event table.

## Re-running

```bash
cd analysis/grenenet_selection/r1_sv_negative_selection
sbatch _submit_graffite_vcf.sh                          # 1 → sv_sites.vcf
sbatch _submit_graffite.sh                              # 2 → graffite/run/, nxf_work
python _graffite_prep_skipped.py                        # 3a (env basic) → graffite/raw/, graffite/skipped/
sbatch --exclude=n0041.savio4 _submit_graffite_skipped.sh   # 3b
git -C ../../../telib/graffite_src apply $PWD/graffite_prepTSD_flank_names.patch   # 5b, once
sbatch _submit_graffite_tsd.sh                          # 5b → graffite/run_tsd/
python _graffite_classify.py                            # 4 (env basic) → sv_repeat_class_*.csv (reads run_tsd if present)
# parity check of 3b: SLURM_ARRAY_TASK_ID=5 TEST=300 TEST_FA=$PWD/results/sv_adaptive/graffite/raw/Chr5/indels.fa bash _submit_graffite_skipped.sh
```

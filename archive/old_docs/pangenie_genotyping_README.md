# pangenie_genotyping — direct genotyping of 151 missing GrENE-Net founders

**Goal**: replace Beagle imputation with **PanGenie direct genotyping**. The 151
GrENE-Net founders that lack long-read assemblies (and thus aren't in the
cactus pangenome) get genotyped at the cactus catalog from their public 1001
Genomes short reads. No imputation step → no 1.43× scale bias → cleaner per-
record alt_freq predictions in the production pool-seq pipeline.

## Why PanGenie over vg giraffe

PanGenie is k-mer based, purpose-built for re-genotyping pangenome variants
from short-read samples. ~30 min – 1h per Arabidopsis sample (vs 3-5h for
vg giraffe). Already integrated in our pipeline (we use its k-mer indices in
`build_kmer_cn.py`). Better SV accuracy than mapping-based genotypers per
Ebler 2022 (Nat Genet).

vg giraffe could supplement later for novel-variant discovery (SNPs / small
indels not in the catalog), but per-base it costs ~10× more compute and the
catalogue is already ~95% complete on this panel.

## Architecture

```
INPUT 1: cactus pang_69 (135 founders × all variants)
         /home/tbellagio/scratch/pang/pang_1001gplus/pang_all/output/
         (currently building, ETA ~3-4 days from 2026-04-29)

INPUT 2: 1001G short-read fastqs for 148 of 151 missing GrENE-Net founders
         from PRJNA273563 — see data/ena_manifest.tsv
         Download size: ~169 GB

PROCESS: For each of the 148 short-read samples:
         1. jellyfish count k-mers from fastqs
         2. PanGenie genotype against pang_69 → per-sample VCF with GT+DS+GQ
         3. Restrict to the 231 GrENE-Net founders' final catalog

OUTPUT:  151-founder VCF (148 PanGenie-genotyped + 3 ./. entries for unavailable),
         then merge with cactus_135 raw VCF →
         286-founder catalog (135 cactus + 151 short-read).
         Filter to 231 GrENE-Net founders → final cn_var/cn_kmer source.
```

## Files

```
data/
├── missing_151_ecotypes.txt  — IDs of the 151 GrENE-Net founders not in cactus_82
├── ena_manifest.tsv          — 148 ENA runs with FTP URLs, MD5s, base counts
└── ena_missing.tsv           — 3 ecotypes for which no suitable WGS run found

scripts/
├── build_ena_manifest.py     — query ENA REST API to populate ena_manifest.tsv
└── (TBD: download_reads.sh, run_pangenie.sh, merge_vcfs.sh)
```

## Coverage of public data

| Source | n | Notes |
|---|---|---|
| ENA: PRJNA273563 (Atwell 2015, 1001G main) | 148 | sample_alias = ecotype_id |
| ENA: PRJNA30811 (Cao 2011, 1001G earlier batch) | 1 | 9977 → SRR095751 / ICE226 (mapping from xwu's worklog) |
| xwu's BAMs (`/carnegie/nobackup/scratch/xwu/GrENE_net/vcf/dedup_<id>.bam`) | 2 | 100001 = `80acc_Cao_S89`, 100002 = `80acc_Cao_S95` — internal Cao 2011 sequencing, never publicly deposited under those IDs |
| **Total covered** | **151 / 151** | 100% coverage |

### Are the SRA fastqs raw or BAM-derived?

**They're raw.** Verified via ENA fields:
- `submitted_format` is empty for both PRJNA273563 and PRJNA30811 entries
- `submitted_ftp` is empty (no separate "submitted format" file beside the fastq)
- `fastq_ftp` IS the original raw FASTQ submission (predates widespread BAM-only submissions)

This matters because BAM-derived fastqs would have unmapped reads dropped → losing reads in SV regions, repetitive areas, novel insertions. Our raw fastqs preserve all reads.

### Are xwu's BAMs OK to use?

Yes, with one caveat: they're already aligned to TAIR10 + dedup'd, but they
**preserve unmapped reads** (verified: 268,212 unmapped of 11.4M total in
`dedup_100001.bam` = 2.3%). `samtools fastq -F 0` extracts ALL reads
(mapped + unmapped) → equivalent to having raw fastqs. The dedup step removed
PCR duplicates which is fine for genotyping.

### Where are the truly raw fastqs for 100001/100002?

The BAM @PG trail shows the source was `trimmed/80acc_Cao_S89-N.{1,2}.fastq.gz`
(and S95 for 100002), processed on the Caltech HPC at:
```
/central/groups/carnegie_poc/lczech/grenephase1/...
```
That path is not mounted from Carnegie. **The raw fastqs aren't on Carnegie**
— xwu transferred only the BAMs, not the raw fastqs, when migrating data here.

To get truly raw fastqs for these two, reach out to Lucas Czech (lczech) /
xwu for a Caltech-side transfer. Otherwise, the practical "same processing"
path is:

```
ENA path (149 ecotypes):
  raw fastq → Trimmomatic → Clumpify dedup → jellyfish count → PanGenie

xwu BAM path (100001, 100002):
  BAM → samtools fastq -F 0 → jellyfish count → PanGenie
       (BAM is already Trimmomatic+MarkDuplicates output;
        functionally equivalent to "post-trim + post-dedup" stage)
```

Both paths converge at the same logical stage ("trimmed + deduplicated
paired-end reads") before PanGenie ingests them. PanGenie is k-mer based →
doesn't care about pre-existing alignment.

### Naming-convention caveat

The 1001 Genomes Project deposited data across multiple ENA studies with
**different sample alias conventions**:
- **PRJNA273563** (Atwell et al. 2015): `sample_alias = ecotype_id` (e.g. "9542")
- **PRJNA30811** (Cao et al. 2011): `sample_alias = ICE_name` (e.g. "ICE226" for ecotype 9977)

The 9977 → SRR095751 / ICE226 mapping was found in xwu's GrENE-Net worklog:
`bwa mem -t 8 -R '@RG\tID:SRR095751\tSM:9977' ...`

### Note on naming conventions

The 1001 Genomes Project deposited its data across multiple ENA studies with
**different sample alias conventions**:
- **PRJNA273563** (Atwell et al. / 1001G main batch, 2015): `sample_alias = ecotype_id`
  (e.g., "9542"). Most of our 148 hits.
- **PRJNA30811** (Cao et al. 2011, earlier 1001G batch): `sample_alias = ICE_name`
  (e.g., "ICE226" for ecotype 9977). Manual override required.

The mapping from ecotype 9977 → SRR095751 was found in xwu's GrENE-Net worklog:
`bwa mem -t 8 -R '@RG\tID:SRR095751\tSM:9977' ...`
That's the only ecotype with this alias-mismatch issue documented in xwu's
notes; the other 7 single-fastq runs (9941, 9966, 9978, 9985, 10011, 10013,
10014) are in PRJNA273563 with matching aliases but stored as interleaved
single-files (paired layout, no `_1`/`_2` split).

For the 3 missing: we'd treat them as `./.` in the final cn_var (they exist
in the GrENE-Net VCF as SNPs only — we just can't genotype them at SVs).
This is a small mass loss (3/231 ≈ 1.3%) and can be handled the same way
the original 82-founder pipeline handled "missing mass".

## Storage / compute estimate

| Stage | Resource | Wall time |
|---|---|---|
| Download 148 fastqs | ~169 GB disk, ENA FTP | ~6-12 h (single-threaded), <1h with aria2c -x16 |
| Jellyfish k-mer count (148 × ~10 GB fastqs) | 8 cores, 16 GB RAM each | ~10-20 min/sample → SLURM array |
| PanGenie genotype (148 × pang_69 graph) | 8 cores, 50 GB RAM each | ~30-60 min/sample → SLURM array |
| Merge VCFs | bcftools merge | ~10 min |

Total: ~3-4 days wall clock with reasonable parallelism (24-task SLURM array).

## Sources

- ENA REST API: `https://www.ebi.ac.uk/ena/portal/api/search`
- 1001 Genomes Project (Alonso-Blanco et al. 2016 Cell, PRJNA273563)
- PanGenie (Ebler et al. 2022 Nat Genet, github.com/eblerjana/pangenie)

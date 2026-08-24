# preprocess provenance — trim + dedup parameters

This file documents the EXACT source for every parameter in `preprocess_one.sh`
and `preprocess_loo_one.sh`, so the trace doesn't live only in commit
messages.

## Data source

Both scripts process **1001 Genomes Project individual-accession short-read
fastqs** downloaded from ENA — primarily PRJNA273563 (Atwell 2015, 1001G main
batch, ~25× HiSeq 2000 PE) plus 1 PRJNA30811 (Cao 2011, ICE alias) and 2
xwu-internal Cao 2011 reruns (100001/100002, lane-concat'd locally).

This is **NOT** GrENE-Net Pool-seq data (which has its own canonical pipeline
at `/home/tbellagio/scratch/pang/grenenet_reads/run_trimmomatic.sh` +
`run_clumpify.sh`, with `SLIDINGWINDOW:4:20` and `minAdapterLength=8`).
The 1001G individual-accession data is older, lower-quality at the read ends,
and is processed with a deliberately lighter trim.

## Source of truth

**`/carnegie/nobackup/scratch/xwu/GrENE_net/vcf/sra/commands.sh`** — Xing Wu's
1001G short-read processing pipeline. Worklog dated 2022-12-13; used for
ecotypes 9940, 9977, 9992, 100001, 100002 (the GrENE-Net founders that needed
to be added to the SNP VCF from public 1001G data).

That file is the canonical reference for processing 1001G individual-accession
data on this machine.

## Trimmomatic

```
java -jar trimmomatic.jar PE -phred33 -threads 5 \
  R1.fastq.gz R2.fastq.gz \
  R1_P.fq.gz R1_U.fq.gz R2_P.fq.gz R2_U.fq.gz \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:True \
  LEADING:5 TRAILING:5 MINLEN:36
```

Parameter-by-parameter:

| param | value | rationale |
|---|---|---|
| Adapter file | TruSeq3-PE-2.fa | TruSeq Y-adapter + reverse-complement variants; correct for HiSeq 2000 / Genome Analyzer II era 1001G data |
| ILLUMINACLIP `seedMismatches` | 2 | Trimmomatic default, allows 2 mismatches in seed |
| `palindromeClipThreshold` | 30 | match score ≥30 in palindrome mode for paired-end overlap detection |
| `simpleClipThreshold` | 10 | min match score for non-palindrome adapter detection |
| `minAdapterLength` | **2** | catch even 2-bp adapter overhangs at read end |
| `keepBothReads` | True | keep R2 even when palindrome detection trims R1's adapter readthrough |
| LEADING | 5 | trim leading bases with Q<5 |
| TRAILING | 5 | trim trailing bases with Q<5 |
| MINLEN | 36 | discard reads shorter than 36bp after trimming |
| **SLIDINGWINDOW** | **(omitted)** | xwu deliberately omits this for 1001G data; the older read-end quality means SLIDINGWINDOW:4:20 would discard too much |

The `2:30:10:2:True` ILLUMINACLIP is more aggressive on tiny adapter remnants
than the GrENE Pool-seq pipeline (`2:30:10:8:TRUE`). Both are reasonable; for
1001G short reads we use xwu's choice for direct provenance.

## Dedup

xwu used **Picard MarkDuplicates** (BAM-stage; requires alignment first):

```
java -jar picard.jar MarkDuplicates I=merged_<eco>.bam O=dedup_<eco>.bam \
    M=marked_dup_metrics_<eco>.txt
```

Source: `/carnegie/nobackup/scratch/xwu/GrENE_net/vcf/sra/submission.txt` and
the `dedup_100001.bam` / `dedup_100002.bam` files at
`/carnegie/nobackup/scratch/xwu/GrENE_net/vcf/`.

We use **Clumpify (BBTools)** instead — fastq-stage, sequence-identity-based:

```
clumpify.sh in1=R1_P.fq.gz in2=R2_P.fq.gz out1=R1_P_dedup.fq.gz out2=R2_P_dedup.fq.gz \
    dedupe=t dupesubs=0 optical=f -Xmx30g
```

Why the swap: PanGenie is k-mer based and does not need an alignment. Picard
requires an aligned BAM — that's an extra alignment step we'd just throw away.
Clumpify is the standard fastq-stage analog that other GrENE-Net pipelines use
(`/home/tbellagio/scratch/pang/grenenet_reads/run_clumpify.sh`); the sequence
identity rule is functionally equivalent to Picard's outer-coordinate rule for
the purpose of "keep one copy per unique read pair".

## Diff between the main and LOO scripts

`preprocess_one.sh` and `preprocess_loo_one.sh` use the same trim+dedup params.
They differ only in:

- input manifest path (`ena_manifest.tsv` vs `loo_ena_manifest.tsv`)
- output directory (`data/preprocessed/` vs `data/loo_preprocessed/`)
- output filename suffix (`_1.dedup.fq.gz` vs `_1P_dedup.fq.gz`, the latter
  matches the canonical GrENE-Net naming)
- preprocess_loo_one.sh activates a single `pang` conda env (provides both
  trimmomatic and clumpify); preprocess_one.sh uses hard-coded jar paths from
  `sequencing_pipeline` env for trim and `pang` env for clumpify.

## In-flight inconsistency

SLURM array job **58045** (main preprocess of the 151 missing accessions) was
submitted before the threshold patch landed. Its already-completed tasks used
the GrENE Pool-seq trim (`minAdapterLength=8`, `SLIDINGWINDOW:4:20`). Per
SLURM semantics, sbatch captures the script content at submit time, so those
tasks are locked into the older parameters.

The over-trim is conservative — it removes some valid bases at low-quality
read ends but does not introduce errors. PanGenie's k-mer matching tolerates
the slightly lower per-read coverage that results.

Re-running 58045 with the new params would cost ~11h wall-clock (180 tasks ×
~30 min ÷ 8 concurrent) and was judged not worth it. The LOO pipeline (jobs
58131 + 58132) and any future re-run of `preprocess_one.sh` will use the
patched parameters.

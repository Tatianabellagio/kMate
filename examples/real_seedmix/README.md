# Real-data example — a sequenced SEEDMIX pool

This is the **Tier 2** example: kMate on a *real* pooled sample against the *real*
production panel. (For a 1-second offline check that your install works at all,
run `kmate selftest` first — that uses a tiny bundled fixture and needs no
external data.)

## What it does

[`run_seedmix.sh`](run_seedmix.sh) estimates per-record allele frequencies in a
real **GrENE-Net SEEDMIX** pool — a sequenced mixture of the 231 *Arabidopsis
thaliana* founders — against the production 231-founder `arch3` panel, on Chr1,
in **global** mode (SEEDMIX is a founder/F0 seed pool, so one founder mixture per
chromosome). SEEDMIX is the canonical validated pool: kMate's recovered founder
mixture agrees with the independent hapFIRE truth at mean *r* ≈ 0.98.

## Inputs (not in the git repo)

Unlike the bundled self-test, this example uses large files that live on the
cluster, not in version control:

| input | path | size |
|---|---|---|
| pool reads (R1/R2) | `…/pang/grenenet_reads/seed_mix/S1-1.{1,2}_P.fq.gz` | ~2.7 GB each |
| panel `kmer_pa` (Chr1) | `data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1.*` | ~2 GB |
| panel `var_pa` (Chr1) | `panel/arch3/chr1/var_pa_231_arch3_chr1.*` | ~0.9 GB |

These paths are the **template** to adapt to your own panel and reads. Build your
own panel matrices with `kmate build-var-pa` / `kmate build-kmer-pa` (see the top-level
README → *Building a panel*).

## Run it

```bash
cd examples/real_seedmix
sbatch run_seedmix.sh             # SEEDMIX sample S1
SAMPLE=S2 sbatch run_seedmix.sh   # a different sample (S1..S8)
```

The dense Chr1 `kmer_pa` for 231 founders needs ~80 GB RAM (the script requests
`--mem=80G`). **Submit with `sbatch`** — do not run on a login node. Chr1 takes a
few minutes once the reads are scanned. Output:

- `results/SEEDMIX_S1_chr1.tsv` — per-record `chrom pos ref_len alt_len alt_freq info n_called se`
- `results/SEEDMIX_S1_chr1.h_per_chrom.npz` — the recovered founder mixture `h`

## Count once, query per chromosome (all 5 chroms)

For a whole-genome run, scanning the multi-GB reads once and querying each
chromosome is much faster than re-scanning per chrom:

```bash
kmate build-kmer-db --reads R1.fq.gz R2.fq.gz --out pool.jf --threads 4
kmate run --kmer-pa-prefix … --var-pa … --var-meta … \
          --reads R1.fq.gz R2.fq.gz --kmer-db pool.jf \
          --sample S1 --out S1.tsv --chroms Chr1 Chr2 Chr3 Chr4 Chr5 \
          --block-mode global --kmer-weight inv_mb
```

Counts from the DB are byte-identical to the per-chrom count path.

## Validate against truth

The recovered founder mixture can be checked against the independent hapFIRE
SEEDMIX truth:

```bash
python ../../scripts/validate_seedmix_vs_hapfire.py
```

(See that script for the truth-file location; hapFIRE truth lives at the
`fc_moilab` project path, symlinked at `~/grenenet-phase1`.)

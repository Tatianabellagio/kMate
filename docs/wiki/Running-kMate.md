# Running kMate

One pool at a time.

```bash
kmate run \
    --kmer-pa-prefix panel/kmer_pa/kmer_pa \
    --var-pa-prefix  panel/var_pa/var_pa \
    --reads R1.fq.gz R2.fq.gz \
    --sample MYPOOL --out MYPOOL_Chr1.tsv \
    --threads 8 --chroms Chr1 --unit chrom
```

Both `--*-prefix` options are **prefixes, not files**. kMate appends the chromosome
and suffix itself, so with `--chroms Chr1` it looks for:

```
panel/kmer_pa/kmer_pa_Chr1.kmer_pa.npz      the k-mer evidence
panel/kmer_pa/kmer_pa_Chr1.meta.npz
panel/var_pa/var_pa_Chr1.var_pa.npz         the projection target
panel/var_pa/var_pa_Chr1.var_called.npz     where each haplotype had a call
panel/var_pa/var_pa_Chr1.meta.npz           the variants themselves
```

If a file is missing, kMate says which one. (The individual `--var-pa`,
`--var-called` and `--var-meta` options still exist if your files are not named to a
common stem.)

This writes the frequency table — see [Output](Output) — plus a `.h_per_chrom.npz`
holding the estimated founder-haplotype mixture.

---

## The one setting you must think about: `--unit`

`--unit` decides **over what span** the founder-haplotype mixture is estimated — a whole
chromosome, or a window. Choose it from your species' biology, not from the default.

| your pools are… | use | why |
|---|---|---|
| **selfing, inbred, or a founder mix** (F0 seed pools) | `--unit chrom` | ancestry is essentially constant along a chromosome, so using all its k-mers gives the best-determined mixture |
| **recombinant** — a few generations of outcrossing | `--unit bp --window-bp 10000 --kmer-weight inv_mb` | ancestry is a mosaic along the chromosome and must be fitted locally |

This matters. On a selfing benchmark, fitting per LD-block instead of per chromosome
moved the error from 0.0033 to 0.0080 and outliers from 0.001% to 0.625%, because
low-diversity regions such as centromeres cannot tell founder haplotypes apart from
local k-mers alone. In the other direction, a chromosome-wide fit on a recombinant pool
averages away the very mosaic you are trying to measure.

**Window size is also limited by sequencing depth** — a window needs enough *observed*
k-mers to be fittable at all. → **[Haplotypes and windows](Haplotypes-and-windows)**
explains both constraints and how to check you chose well.

Leave `--normalize per_founder` on. It corrects for founder haplotypes differing in how
many k-mers they contribute; without it, k-mer-poor haplotypes collapse toward zero.

---

## Whole genome, counting reads once

Counting k-mers is the slow part. Do it once per pool and reuse it for every
chromosome:

```bash
kmate build-kmer-db --reads R1.fq.gz R2.fq.gz --out pool.jf --threads 8

for CHR in Chr1 Chr2 Chr3 Chr4 Chr5; do
  kmate run \
      --kmer-pa-prefix panel/kmer_pa/kmer_pa \
      --var-pa-prefix  panel/var_pa/var_pa \
      --reads R1.fq.gz R2.fq.gz --kmer-db pool.jf \
      --sample MYPOOL --out MYPOOL_${CHR}.tsv \
      --threads 8 --chroms $CHR --unit chrom
done
```

Results are identical to re-scanning the reads each time, and about twice as fast.

---

## Many samples

One job per sample. Here is a complete SLURM array script — adapt the header to your
cluster:

```bash
#!/bin/bash
#SBATCH --job-name=kmate
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --array=1-100          # one task per sample
set -uo pipefail

PANEL=/path/to/panel
MANIFEST=samples.tsv           # header + columns: sample_id, r1, r2
OUT_DIR=results
CHROMS="Chr1 Chr2 Chr3 Chr4 Chr5"

read -r SAMPLE R1 R2 < <(awk -v n=$((SLURM_ARRAY_TASK_ID + 1)) 'NR==n' "$MANIFEST")
mkdir -p "$OUT_DIR"

# k-mer DB on NODE-LOCAL storage -- see the warning below
JF_DIR=${SLURM_TMPDIR:-/dev/shm}
JF=$JF_DIR/$SAMPLE.jf
trap 'rm -f "$JF"' EXIT

kmate build-kmer-db --reads "$R1" "$R2" --out "$JF" --threads 8

for CHR in $CHROMS; do
    OUT=$OUT_DIR/${SAMPLE}_${CHR}.tsv
    [ -s "$OUT" ] && continue          # resumable: skip what is already done
    kmate run \
        --kmer-pa-prefix $PANEL/kmer_pa/kmer_pa \
        --var-pa-prefix  $PANEL/var_pa/var_pa \
        --reads "$R1" "$R2" --kmer-db "$JF" \
        --sample "$SAMPLE" --out "$OUT" \
        --threads 8 --chroms $CHR --unit chrom || exit 1
done
```

Submit with `sbatch --array=1-<N> runner.sh`.

> ### Put the k-mer database on node-local storage
> Never on shared network storage. With many jobs reading their databases off the same
> shared filesystem at once, we measured queries up to **48× slower**. `$SLURM_TMPDIR`
> or `/dev/shm` as above. This is the single biggest performance factor at scale, and
> it is not specific to kMate.
>
> If you use `/dev/shm`, the database counts against your job's memory — request
> ~32 GB. The `trap` frees it even if the job is killed.

The repository also contains the production runner this is distilled from,
`grenenet/run_site_array_perchrom.sh`, which adds manifest chunking past the array-size
limit and finer resume logic.

---

## Memory

Roughly **20 GB** for a few hundred founder haplotypes on one chromosome, because the k-mer
matrix is loaded densely. Request ~32 GB. Fewer haplotypes need much less.

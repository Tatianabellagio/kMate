# Running kMate

One pool at a time, one chromosome at a time.

```bash
kmate run \
    --kmer-pa-prefix panel/kmer_pa/kmer_pa \
    --var-pa     panel/var_pa/var_pa_Chr1.var_pa.npz \
    --var-called panel/var_pa/var_pa_Chr1.var_called.npz \
    --var-meta   panel/var_pa/var_pa_Chr1.meta.npz \
    --reads R1.fq.gz R2.fq.gz \
    --sample MYPOOL --out MYPOOL_Chr1.tsv \
    --threads 8 --chroms Chr1 \
    --unit chrom --kmer-weight uniform --normalize per_founder
```

`--kmer-pa-prefix` is a **prefix**: kMate appends `_<CHR>.kmer_pa.npz` itself.

This writes the frequency table (see [Output](Output)) and a `.h_per_chrom.npz`
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

Counting k-mers is the slow part. Do it once and reuse it for every chromosome:

```bash
kmate build-kmer-db --reads R1.fq.gz R2.fq.gz --out pool.jf --threads 8

for CHR in Chr1 Chr2 Chr3 Chr4 Chr5; do
  kmate run \
      --kmer-pa-prefix panel/kmer_pa/kmer_pa \
      --var-pa     panel/var_pa/var_pa_${CHR}.var_pa.npz \
      --var-called panel/var_pa/var_pa_${CHR}.var_called.npz \
      --var-meta   panel/var_pa/var_pa_${CHR}.meta.npz \
      --reads R1.fq.gz R2.fq.gz --kmer-db pool.jf \
      --sample MYPOOL --out MYPOOL_${CHR}.tsv \
      --threads 8 --chroms $CHR --unit chrom
done
```

Results are identical to re-scanning the reads each time, about twice as fast.

---

## Many samples

Run one job per sample. Two things matter at scale:

- **Count each pool's k-mers once** (`--kmer-db`), as above.
- **Put the k-mer database on node-local or RAM storage**, never on shared network
  storage. With many jobs reading their databases off the same shared filesystem at
  once, we measured queries up to **48× slower**. This is the single biggest
  performance factor, and it is not specific to kMate.

The repository contains a SLURM array runner (`grenenet/run_site_array_perchrom.sh`)
that does both; point it at any panel with environment variables.

---

## Memory

Roughly **20 GB** for a few hundred founder haplotypes on one chromosome, because the k-mer
matrix is loaded densely. Request ~32 GB. Fewer haplotypes need much less.

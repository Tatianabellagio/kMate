# Functional tracks for candidate loci: eQTL, ATAC, TFBS

How the three functional tracks in the MOI-LAB zoom-Manhattan figures are built, and what
changes when we apply them to our non-SNP candidates. Written 2026-09-17 from the shared
Drive sources (`MOI-LAB/PROJECTS/syntheticevolution/`), for the kMate GrENE-Net candidates
in `results/screen_visual_review_round2.csv`.

Upstream files read (Drive):

| file | role |
|---|---|
| `scripts/zoom_manhattan/plot_zoom_manhattan.py` | the composite figure (gene model / GEA / eQTL / GPN) |
| `scripts/zoom_manhattan/extract_viewer_tracks.py` | harvests tracks out of the plotly browser HTMLs |
| `scripts/zoom_manhattan/README.md` | data sources and known gaps |
| `scripts/run_eqtl.py` + `pipelines/eqtl/run_eqtl_one_gene.R` | the cis-eQTL computation |
| `scripts/NARROW_MYB6_tfbs_fimo.py` | the FIMO ref-vs-alt TFBS turnover analysis |

⚠ `extract_viewer_tracks.py` does NOT compute anything: it scrapes per-SNP arrays out of
pre-built plotly HTMLs. The real recipes are the two scripts below plus one static file.

---

## 1. cis-eQTL ("exp. QTL" panel)

`run_eqtl.py --agi <AGI>` wraps `run_eqtl_one_gene.R`. Per gene:

- **Expression**: 1001T transcriptomes, `TG_data_20180606.RData`.
  `TG.genes$d_log2_batch` (batch-corrected log2) for the gene level, `TG.trans$d_log2`
  for isoforms, sample ids from `TG.meta$index`; samples with `batch_comb == "MU"` dropped.
- **Genotypes**: PLINK `1001gbi.{bed,bim,fam}` (1001 Genomes biallelic SNPs).
- **cis window**: gene body ±10 kb (`--window`, default 10000), cut with
  `plink2 --chr N --from-bp --to-bp --make-bed`.
- **Association**: phenotype written into column 6 of the temporary `.fam` (missing = -9),
  requires >= 50 matched samples, then `gemma -lm 4` (linear model, all test statistics).
  Output `<AGI>_<pheno>.assoc.txt` with `rs, ps, beta, se, af, p_wald`; the figure plots
  `-log10(p_wald)`.
- **Three phenotypes per gene**
  1. gene level;
  2. one per isoform (`AT1G01010.1`, `.2`, …);
  3. **isoform ratio** = logit(isoform / sum of isoforms), with the bottom 5% of total
     expression masked to NA — i.e. a splicing QTL.
- **Failure mode to expect**: a gene not expressed in 1001T has no phenotype and therefore
  no eQTL (their aquaporin AT1G52180: TPM > 0 in 10 of 874 samples).

**For us.** Nothing about this is SNP-specific: the phenotype is expression, the genotypes
are whatever PLINK set we hand it. Two options, in increasing order of work:
*(a)* run it as-is on the 1001g SNPs and ask whether our indel/SV sits inside a cis-eQTL
peak — cheap, but the eQTL is then SNP-tagged by construction and says nothing our
`best_r2_snp` does not; *(b)* add our non-SNP genotypes to the PLINK set so the indel is
itself tested for expression association — the informative version for a kMate candidate,
since the question is whether the indel, not its SNP proxy, changes expression.

## 2. ATAC ("open chromatin" spans)

Not computed. A static multi-tissue peak union, `colab_data/ATAC-seq_multitissue.csv`
(~2.2 MB; columns `chrom, start, end, number_tissue`), filtered to the plot window:

```python
d = df[(df.chrom == chrom) & (df.end >= lo) & (df.start <= hi)]
```

drawn as pale spans, with `number_tissue` = in how many tissues the peak was called (so it
can be shaded by breadth). Adopting it is one file copy plus an interval join.

**For us.** The useful statistic is not decoration: *does the candidate indel/SV overlap an
ATAC peak, and in how many tissues?* For promoter candidates (FUS3, GPX6, AT2G30000,
AT4G13200) that is a direct regulatory prior, and it is testable against a background of
matched non-candidate indels — i.e. it can be a number, not only a picture.

## 3. TFBS (FIMO, PlantTFDB) — an ALLELE CONTRAST, not a motif scan

`NARROW_MYB6_tfbs_fimo.py`, for a window around the peak:

1. reference sequence via `samtools faidx TAIR10.fa Chr4:5995300-5996050`;
2. a second sequence with the selected (adaptive) alleles substituted, from the `.bim`
   alleles;
3. `fimo --thresh 1e-4 Ath_TF_binding_motifs.meme locus.fa` on **both** sequences at once
   (PlantTFDB *A. thaliana* motifs; `Ath_TF_binding_motifs_information.txt` maps motif ->
   TF family);
4. hits keyed on `(motif, strand, start, end)` and compared between the two sequences:
   - ref only -> **LOST** (adaptive allele destroys the site)
   - alt only -> **GAINED**
   - both -> weakened / strengthened by FIMO score
5. house filters: keep `q < 0.005`, drop low-complexity matches — a homopolymer run >= 5
   or <= 2 distinct bases ("the CAM5 lesson").

Their own note on why the contrast rather than the significance: the SNPs sit in a short
window, so per-window q inflates; the robust signal is the ref-vs-alt **asymmetry** at the
same threshold. In the figure a TFBS box is red when it overlaps a variant with r² >= 0.5
to the focal variant, i.e. it sits on the focal haplotype.

**For us — two changes, one of them substantive.**

- *Coordinate mapping (must fix).* The script maps a FIMO hit back to the genome with
  `gstart = WIN_START + start - 1`, which holds only because a SNV substitution keeps ref
  and alt the same length. For an indel the alt sequence is shifted after the variant, so
  the alt-side hits need an explicit offset map (identity before the variant, ±|size|
  after) or the gain/loss pairing silently compares different places.
- *A deletion is a stronger test than a SNV.* A 2 bp promoter deletion (FUS3) or a 1.2 kb
  deletion (GPX6) removes motif sequence outright, so LOST/GAINED is a statement about
  sequence presence rather than a score shift at one base. This is the one place where our
  variant class is an advantage over the SNP-based original.

---

## What is needed to run any of this

| input | where |
|---|---|
| `Ath_TF_binding_motifs.meme` + `..._information.txt` | PlantTFDB (Drive: `arabidopsispangenomes/data/motifs/`) |
| `ATAC-seq_multitissue.csv` | Drive `colab_data/` |
| `TAIR10.fa` (+ `.fai`) | local (`lib`) |
| `1001gbi.{bed,bim,fam}` | Drive `arabidopsisgenomes/` |
| `TG_data_20180606.RData`, `gene_infoV2.RData` | Drive `data/eqtl/inputs/` |
| `fimo` (MEME suite), `plink2` | install via mamba (`meme`, `plink2`) |
| `gemma` | already here: `~/miniforge3/envs/gwas_tools/bin/gemma` |

Order of work for our candidates: ATAC overlap (cheapest, gives a number), then the
indel-aware TFBS turnover on the promoter candidates, then eQTL last — it is the only one
needing a second genotype panel and a real decision about whether to test the indel itself.

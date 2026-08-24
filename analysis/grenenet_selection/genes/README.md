# genes — which genes, and are they real

Two layers, deliberately kept separate because the second exists to check the
first.

```
genes/
    attribution/   block -> gene, from the per-site GWAS      (tracked)
    dissection/    per-locus validation, from the LFMM-WZA GEA (untracked WIP)
```

| | `attribution/` | `dissection/` |
|---|---|---|
| upstream scan | the **GWAS** (`r3_persite_gwas/`: class-split, multitrait, 31 per-site) | the **GEA** (`r2_gea_nonsnp/phase1_replication/multiaxis/` new-peaks) |
| how a gene is picked | clq0.9 blocks BH-significant in the non-SNP scan but **not** the SNP scan → block interval → overlapping TAIR10 genes | dissect the locus: GEA Manhattan + founder-haplotype panel + founder-LD triangle, sharing a genomic x-axis |
| unit | **block attribution** | **the variant's own position** |
| what it gives you | a candidate list | whether the candidate survives a look at where the signal actually sits |

> **Why they are not merged.** `dissection/` exists *because block attribution
> can be wrong. Its worked case: block Chr2_5291 was attributed to CARK8/CARK9
> (AT2G30740/30730), but the real Bonferroni SV is a 5.1 kb insertion ~29 kb
> away at Chr2:13,127,635, with founder r² ≤ 0.05 to those genes and r² = 0.88
> to a different local block (GASA12 / HVT1). A sparse-SV `merge_small_blocks`
> step had chained empty blocks together, ballooning the block footprint to
> ~40 kb. So an `attribution/` gene list is a **hypothesis**, not a result —
> treat it as one until `dissection/` has looked at the locus.

---

## attribution/

| script | what it does |
|---|---|
| `nonsnp_only_genes.py` | the core: union the non-SNP-only blocks across every contrast (multitrait JOINT/GLOBAL, CLIMATE × bio1–19 + PC1, 31 per-site scans) → clq0.9 interval → overlapping TAIR10 genes |
| `nonsnp_only_genes_describe.py` | annotates that list via public gene APIs |
| `go_enrichment_nonsnp.py` | GO / stress-term over-representation of the non-SNP-only genes |
| `genes_from_regions.py` | helper: genes overlapping arbitrary regions, via Ensembl Plants REST |
| `plot_gi_locus.py` | annotated Manhattan for the GIGANTEA (AT1G22770) non-SNP-only hit |
| `_build_nonsnp_only_genes_nb.py` | → `notebooks/nonsnp_only_genes.ipynb` |

Reads and writes `../../r3_persite_gwas/results/varexp/`, since that is where
the class-split GWAS output lives.

> ⚠ For gene *function*, use `r2_gea_nonsnp/phase1_replication/annotate_genes_tair_uniprot.py`
> (TAIR GO + UniProt), not the older mygene-only annotator — NCBI carries no
> free-text summaries for Arabidopsis loci, so a mygene-only run leaves the
> `summary` column empty and the functional categories become name-driven.

## dissection/

Currently **untracked** work in progress.

| script | what it does |
|---|---|
| `screen_sig_blocks.py` | genome-wide screen of per-class Bonferroni block leads → classify each variant CDS / UTR / promoter / intron / TE / intergenic against the full TAIR10 GFF, and flag `block_gene_mismatch` (the CARK-type re-attribution) |
| `augment_snp_cosig.py` | adds `snp_cosig_2kb` — is there a Bonferroni SNP within 2 kb? separates SNP-shadowed from SV/indel-unique |
| `theme_filter.py` | narrows to climate / flowering / circadian / stress themes |
| `plot_locus_combined.py` | the parametrized combined locus figure for any gene — reuse this for new candidates |
| `plot_cark_combined.py` | the original CARK8/9 figure that `plot_locus_combined.py` generalizes |
| `plot_sv_garden_dynamics.py` | per-garden AF trajectory of one variant over generations, cold vs warm |

Outputs → `dissection/results/` (`loci/`, `regions/`, screens as CSV) with
figures under `dissection/results/plots/`.

> **Caveat carried from the screen:** it runs on raw, uncalibrated LFMM p-values
> throughout. It is a candidate *net*, not a calibrated hit list. The triage
> rule that emerged is recurrence across axes + founder frequency + LD-to-gene —
> **not** p-value.

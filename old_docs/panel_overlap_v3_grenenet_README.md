# Panel overlap — v3 vs GrENE-Net 231 SNP catalog

**Question**: how much of the GrENE-Net 231-founder SNP catalog (`greneNet_final_v1.1.recode.vcf`, the panel hapFIRE consumes natively) is also represented in the v3 cactus_em panel (`founders_231_chr.haploid.vcf.gz`, the panel cactus_em projects through)?

**Why it matters**:
- hapFIRE estimates per-block h on its 3.24 M SNP positions.
- Projecting hapFIRE's h through `cn_var_v3` only produces meaningful per-record AFs at positions present in `cn_var_v3`.
- Positions in GrENE-Net but NOT in v3 → no projection target (silent loss).
- Positions in v3 but NOT in GrENE-Net → hapFIRE doesn't have h there (different gap).

This is the v3 analog of the older test at `/home/tbellagio/scratch/freqk_gr/panel_overlap_test/` (which compared the syri-derived 82-accession panel vs GrENE-Net). Here both panels cover the **same 231 founders**, so no founder-subset complication is needed — direct position overlap.

## Inputs

| panel | path | n SNPs (Chr1-5, biallelic) |
|---|---|---|
| **v3** (this project) | `pangenie_genotyping/data/merged/founders_231_chr.haploid.vcf.gz` | TBD (extract via `extract_v3_positions.sh`) |
| **GrENE-Net 231** | `/carnegie/nobackup/scratch/xwu/GrENE_net/greneNet_final_v1.1.recode.vcf` | 3,235,480 (already extracted at `data/grenenet231_snp_positions.tsv`) |

## Pipeline

1. Submit `extract_v3_positions.sh` (SLURM, ~5-10 min) → produces:
   - `data/v3_all_positions.tsv` — every (chrom, pos) in v3
   - `data/v3_snp_positions.tsv` — biallelic SNP only
   - `data/v3_snp_positions_dedup.tsv` — multi-allelic SNP positions, dedup'd
   - `.nochr.tsv` versions with 'Chr' prefix stripped to match GrENE-Net

2. Open `overlap_v3_vs_grenenet.ipynb` — Venn diagrams + counts + per-chrom breakdown.

## Chrom-naming caveat

GrENE-Net VCF uses chromosome names `1, 2, 3, 4, 5` (no prefix). v3 VCF uses `Chr1, Chr2, ...`. The extract script produces `.nochr.tsv` versions for direct comparison.

## Files

| File | Purpose |
|---|---|
| `extract_v3_positions.sh` | SLURM job to extract v3 position lists |
| `data/v3_*_positions.tsv` (+ `.nochr.tsv`) | v3 (chrom, pos) lists |
| `data/grenenet231_snp_positions.tsv` | copied from old test — 3.24M SNP positions |
| `overlap_v3_vs_grenenet.ipynb` | Venn diagrams + per-chrom breakdown |
| `results_summary.md` | Headline numbers (filled in after notebook runs) |

# benchmarks/vs_hapfire/ — hapFIRE's inputs

hapFIRE's pools and runs for the paper figures are in `../simulations/` (founder sets `greneNet231` / `greneNet78`, fits
`<pool>.hapfire`, scored with kMate's in `score.py`). Each tool runs on its own panel: hapFIRE needs a fully imputed SNP
panel, and a panel shared with kMate would have to call missing genotypes reference (`PANEL_MISMATCH_BUG.md`). Its reads
come from its own SNP-only genomes, drawn with the same founders as kMate's pool (checked by `score.py`).

| path | what | built by (`scripts/`) |
|---|---|---|
| `work/greneNet_chr1.vcf.gz` | the GrENE-Net VCF (231 founders, SNPs, fully imputed), Chr1 | `bcftools view -r 1` of `work/greneNet_final_v1.1.recode.vcf.gz` (recorded in its header) |
| `work/greneNet_chr1_p78.vcf.gz` | the same cut to the 78 long-read founders, sites still polymorphic among them | `build_greneNet_p78_vcf.sh` |
| `greneNet_fastas/<founder>.chr.fa` | each founder's Chr1 from the GrENE-Net VCF (TAIR10 coordinates, SNPs only) | `build_greneNet_fastas.sh` |
| `work/greneNet_gt_matrix.npz` | founder dosages of the GrENE-Net SNPs, for hapFIRE's allele-frequency truth | `extract_greneNet_gt_matrix.py` (`extract_gt_matrix.sbatch`) |
| bwa index of `ref_xing` TAIR10 | the reference hapFIRE's reads are aligned to | `build_bwa_index.sh` |

`ecotype_resolution/` is a separate analysis with its own README. The scripts that ran the comparison before, its old
results and the previous README are in `../archive/scripts_before_simulations_2026-10-17/vs_hapfire/`; the old pools and
runs (`sims_greneNet*`, `results/greneNet_fair*`) are listed in `../simulations/logs/old_folders_to_delete.txt`.

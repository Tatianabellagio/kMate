# benchmarks/p80/ — inputs of the long-read (78-founder) benchmarks

The simulations and fits that use these inputs run from `../simulations/` (run table `runs.py`, figures `longread_grid`
and `units_p78`). The scripts that ran them before, the old results, the coordinate-bug audit (`audit_sims/`) and the
old README are in `../archive/scripts_before_simulations_2026-10-17/p80/`.

| path | what | built by (`scripts/`) |
|---|---|---|
| `data/pangenome_p80_chr1.vcf.gz`, `data/var_pa_p80.*`, `data/kmer_pa_p80/` | the 80 long-read assemblies' Chr1 variant and k-mer panels | `01_build_biallelic_p80.sh`, `04_build_var_pa_p80.sh`, `03_build_kmer_pa_p80.sh` |
| `data/var_pa_p78.*`, `data/kmer_pa_p78_{raw,filt2inv}/` | the same cut to the 78 long-read founders of the 231 panel (`data/samples_78_acc.txt` = `longread_founders_78.txt`), records and k-mers no longer variable among them dropped; `filt2inv` = the private k-mer filter for 78 founders | `build_p78_panels.py` |
| `fastas_78n/<founder>.chr.fa`, `.chain` | each of the 78 founders' Chr1 from `data/pangenome_p80_chr1.vcf.gz`, `bcftools consensus -H 1 -M N`, with its TAIR10 → founder chain | `05n_build_fastas_p78n.sh` |
| `data/var_pa_p78.ld_blocks_r2_0.10_Chr1.tsv` | r² ≥ 0.1 LD blocks of `var_pa_p78` (`--unit ld`) | `../simulations/build_ld_blocks.py` |

Not used any more and kept only until deleted: `fastas_80/`, `data/kmer_pa_p80_{filt2,ours,ours_filt2inv}/`,
`data/var_pa_p80.ld_blocks_r2_0.10_Chr1.tsv`, and the old simulation and result folders (`sims*/`, `results/`; listed in
`../simulations/logs/old_folders_to_delete.txt`).

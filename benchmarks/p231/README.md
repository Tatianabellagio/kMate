# benchmarks/p231/ — inputs of the 231-founder benchmarks

The simulations and fits that use these inputs run from `../simulations/` (run table `runs.py`). The scripts that ran
them before, the old results and their README are in `../archive/scripts_before_simulations_2026-10-17/p231/`.

| path | what | built by (`scripts/`) |
|---|---|---|
| `fastas_231n/<founder>.chr.fa` | each founder's Chr1 as the panel VCF (`panel/arch3/chr1/merged_231_chr1_final.vcf.gz`) records it, `bcftools consensus -H 1 -M N` (missing genotypes as N, as the k-mer panels treat them) | `05n_build_fastas_p231n.sh` |
| `fastas_231n/<founder>.chain` | TAIR10 → founder coordinates, from the same `bcftools consensus` call (`.chaincheck`: regenerated genome identical) | `05c_build_chains_p231n.sh` |
| `data/kmer_pa_p231/` | k-mer panel of the 231 founders, every k-mer (founder-type figure, "all k-mers") | `03_build_kmer_pa_p231.sh` |
| `data/kmer_pa_p231_filt2inv/` | the same panel with the private k-mer filter (k-mers carried by 2 to 230 founders) | `03c_build_kmer_pa_filt2inv_p231.sh` |
| `data/ld_blocks_r2_0.10_Chr1.tsv` | r² ≥ 0.1 LD blocks of the production variant panel (`--unit ld`) | `../simulations/build_ld_blocks.py` |

The production panels the main runs use (`data/kmer_pa_231_arch3_filt2inv`, `panel/arch3/chr1/var_pa_231_arch3_chr1*`)
come from the arch3 panel pipeline (`panel/arch3/`). Not used any more and kept only until deleted: `fastas_231/`
(genomes with missing genotypes as reference), `data/kmer_pa_p231_filt2/`, and the old simulation and result folders
(`sims*/`, `results*/`, listed in `../simulations/logs/old_folders_to_delete.txt`).

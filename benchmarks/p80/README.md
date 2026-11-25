# benchmarks/p80/ — inputs of the long-read (78-founder) benchmarks

The simulations and fits that use these inputs run from `../simulations/` (run table `runs.py`, figures `longread_grid`
and `units_p78`). The scripts that ran them before, the old results, the coordinate-bug audit (`audit_sims/`) and the
old README are in `../archive/scripts_before_simulations_2026-10-17/p80/`.

| path | what | built by (`scripts/`) |
|---|---|---|
| `data/pangenome_p80_chr1.vcf.gz`, `data/var_pa_p80.*` | the 80 long-read assemblies' Chr1 variants and variant panel | `01_build_biallelic_p80.sh`, `04_build_var_pa_p80.sh` |
| `data/kmer_pa_p80_ours/` | their k-mer panel, built with the production (in-house) k-mer index; the source of the 78-founder k-mer panels | `03d_build_kmer_pa_ours_p80.sh` (also writes `kmer_pa_p80_ours_filt2inv/`) |
| `data/var_pa_p78.*`, `data/kmer_pa_p78_{raw,filt2inv}/` | `var_pa_p80` and `kmer_pa_p80_ours` cut to the 78 long-read founders of the 231 panel (`data/samples_78_acc.txt` = `longread_founders_78.txt`), records and k-mers no longer variable among them dropped; `filt2inv` = the private k-mer filter for 78 founders | `build_p78_panels.py` |
| `fastas_78n/<founder>.chr.fa`, `.chain` | each of the 78 founders' Chr1 from `data/pangenome_p80_chr1.vcf.gz`, `bcftools consensus -H 1 -M N`, with its TAIR10 → founder chain | `05n_build_fastas_p78n.sh` |
| `data/var_pa_p78.ld_blocks_r2_0.10_Chr1.tsv` | r² ≥ 0.1 LD blocks of `var_pa_p78` (`--unit ld`) | `../simulations/build_ld_blocks.py` |

Not used by the paper but read by the exploratory benchmarks (`../founder_h_uncertainty/`, `../unit_localonly_p231/`,
`../vs_vg_giraffe/`, `../vs_hapfire/ecotype_resolution/`): `data/kmer_pa_p80/` (`03_build_kmer_pa_p80.sh`),
`data/kmer_pa_p80_filt2/` (`03b_build_kmer_pa_filt2_p80.sh`), `data/kmer_pa_p80_ours_filt2inv/`, `data/var_pa_p80.ld_blocks_r2_0.10_Chr1.tsv`.

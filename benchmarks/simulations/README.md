# benchmarks/simulations

Every simulated pool and every kMate and hapFIRE fit behind the paper's kMate figures, from one simulator, one fit
runner and one scorer.
The figures are drawn by `analysis/grenenet_selection/notebooks/_build_kmate_benchmarks_nb.py` from `results/`.

| file | what it does |
|---|---|
| `config.py` | paths (founder genomes, k-mer and variant panels, tools), read settings, Slurm resources |
| `runs.py` | the run table: every pool and fit, grouped by figure (`python runs.py` prints the counts) |
| `simulate.py POOL` | one pool: mosaic genomes, seeded wgsim reads, per-record truth → `pools/POOL/` |
| `fit.py FIT` | one kMate or hapFIRE fit, timed → `fits/FIT/` |
| `submit.py [FIGURE ...]` | submits every pool and fit not done yet, then `score.py` |
| `score.py` | scores every fit, writes `results/*.tsv` |

Names carry the parameters: pool `p231_n50_g3_self97_s42_10x` = 50 individuals from the 231 founders, 3 generations of
97% selfing, seed 42, 10×; fit `p231_n50_g3_self97_s42_10x.p231.chrom` = that pool, production k-mer panel, one mixture
per chromosome. A pool or fit listed by two figures is run once. hapFIRE's pools (`greneNet231_*`, `greneNet78_*`) hold
the same founders as kMate's, simulated from its own GrENE-Net SNP-only genomes; its fits (`<pool>.hapfire`) align
with bwa and run hapFIRE on its GrENE-Net VCF (`../vs_hapfire/scripts/` builds those inputs).

## The simulation

1. **Mosaics** (`sims/scripts/make_recomb_mosaics.py`). n individuals drawn without replacement from the founder set's
   list with the pool seed (so a hapFIRE pool, drawn from the same list, holds the same founders as its kMate pool), then
   `gen` generations: each offspring a clone of one parent with probability `selfing`, otherwise a recombinant of two
   (Poisson crossovers, 4 cM/Mb). Each TAIR10 segment is cut from its founder at the founder's own coordinates (the
   chain written with the founder genome), so a whole founder keeps its whole sequence (hapFIRE's SNP-only genomes are
   in TAIR10 coordinates and need none). Founder genomes are each
   founder's Chr1 as the panel records it, missing genotypes written as N (`../p231/scripts/05n_build_fastas_p231n.sh`,
   `../p80/scripts/05n_build_fastas_p78n.sh`; chains `05c_build_chains_p231n.sh`).
2. **Reads** (wgsim, seeded per individual from the pool name). Individual i has pool share w_i (1/n, or `dominant`
   for individual 1 and the rest equal) and gets round(depth × w_i × L_i / 300) pairs of 150 bp, L_i its length without N:
   the same per-base depth along every individual. Error 0.001, mutation 0.001 (15% indels), insert 500 ± 50, as VISOR
   SHORtS ran wgsim for the earlier simulations.
3. **Truth** (`sims/scripts/compute_recomb_truth.py`): per record of each variant panel, Σ w_i × allele of the founder
   individual i carries there, over the individuals whose founder is called at the record.

Every step is deterministic: a pool is reproduced from its name. Mosaic FASTAs are deleted after the reads are drawn;
`shares.tsv` records each individual's share, length, read pairs and wgsim seed.

## Run, from a clean start

Inputs (built once; each script's header gives its usage and environment):

| input | built by |
|---|---|
| 231 founder genomes (`../p231/fastas_231n/`) and their chains | `../p231/scripts/05n_build_fastas_p231n.sh`, `05c_build_chains_p231n.sh` |
| 78 long-read founder genomes + chains (`../p80/fastas_78n/`) | `../p80/scripts/05n_build_fastas_p78n.sh` |
| 231 production panels (`data/kmer_pa_231_arch3_filt2inv`, `panel/arch3/chr1/var_pa_231_arch3_chr1*`) | the arch3 panel pipeline (`panel/arch3/`) |
| founder-type k-mer panels (`../p231/data/kmer_pa_p231{,_filt2inv}`) | `../p231/scripts/03_build_kmer_pa_p231.sh`, `03c_build_kmer_pa_filt2inv_p231.sh` |
| 78-founder panels (`../p80/data/var_pa_p78*`, `kmer_pa_p78_{raw,filt2inv}`) | `../p80/scripts/01_build_biallelic_p80.sh`, `03_build_kmer_pa_p80.sh`, `04_build_var_pa_p80.sh`, `build_p78_panels.py` |
| hapFIRE's genomes, VCFs, genotype matrix, bwa index | `../vs_hapfire/scripts/build_greneNet_fastas.sh`, `build_greneNet_p78_vcf.sh`, `extract_greneNet_gt_matrix.py`, `build_bwa_index.sh` |

Then (`kmate` env, compute node):

```bash
cd benchmarks/simulations
python runs.py                      # what would run
python submit.py --dry-run          # what is missing
python submit.py                    # submit it (skips what is done); score.py runs after every fit
```

and, once `score.py` has written `results/`, the figures (`basic` env, compute node):

```bash
cd analysis/grenenet_selection/notebooks && python _build_kmate_benchmarks_nb.py
```

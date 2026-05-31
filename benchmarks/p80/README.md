# `benchmarks/p80/` — 80-cactus-founder control experiment

> **Naming:** the AF-estimation method is now called **kMate** (see `../ALGORITHM.md`). The `cactus_em_*` result-dir names and `07*_run_cactus_em_*.sh` scripts here retain the legacy name `cactus_em`; a full path rename is deferred.

**Standalone project.** Successor to `control_p82/` after the 2026-05-16 panel-QC decisions dropped 2 cactus assemblies. Tests the same hypothesis: that the +41% cactus h-bias and the off-diagonal streaks in `FINAL_RESULTS_cov10_v3.ipynb` vanish (or shrink) when the panel is homogeneous (all founders are long-read cactus assemblies, no PG short-read genotypes).

See `../HANDOFF_P82_CONTROL.md` for the full hypothesis-test rationale (still applies; only the panel size and decomposition method changed). See `../SESSION_2026-05-16.md` for the 2-assembly exclusion decision.

## What changed vs `control_p82/`

### Panel: dropped 2 assemblies

Per `../data/exclude_list.txt` and `../data/flag_list.tsv`:

| Assembly_ID | Accession_ID | Name  | Reason |
|---|---|---|---|
| 101003 | 5772 | Set-1 | Cactus FASTA mislabeled (k-mer J=0.964 to T980 dup); 1001G J=0.300 |
| 100852 | 9947 | Ped-0 | ~10× ONT; LOO PG-vs-cactus disagreement 30.34% (rank 4/75) |

Net cactus panel: **80 founders.** The 4 untriaged-tail founders (Hau-0/7164, Mt-0/6939, Aitba-1/9606, Brö1-6/8231) and 4 ambiguous-pair founders (Ei-2/St-0, Stw-0/Ct-1) are **kept** — they're `flag-*` not `exclude-*` in `flag_list.tsv`.

### Decomposition: arch (annotate_vcf + convert-to-biallelic), not vcfbub + norm

p82 used `vcfbub -l 0 -r 100000 | bcftools norm -m -any`. p80 follows the production arch3 chr1 cactus-side procedure:

**Shared prerequisite (already built in `../panel/arch3/chr1/`):**
- `chr1_135_annotated.sorted.vcf.gz` — full 135-assembly pangenome, Chr1, with symbolic `INFO/ID` from `annotate_vcf.py -gfa pang_1001gplus_all.gfa.gz` (the production A1 step).
- `chr1_135_annotated_biallelic.sorted.vcf.gz` — biallelic catalog from the same annotate step.

**p80-specific cactus-side step (`01_build_biallelic_p80.sh`):**
1. Subset the 135-asm annotated VCF to the 80 Asm_IDs we want (= 82 from `pang_1001gplus_82acc.vcf.gz` minus `data/exclude_list.txt`). All 80 are confirmed present in the 135-asm catalog.
2. `convert-to-biallelic.py <biallelic catalog>` (the PanGenie helper) → one biallelic row per atomic variant ID; carrier GT derived from each sample's GT against `INFO/ID`'s per-ALT IDs.
3. sort + bgzip + fill-tags AC/AN/F_MISSING + drop AC=0 + reheader Asm_ID → Acc_ID.

No `transfer_id` step (unlike `panel/arch3/chr1/jobA3_*.sh`): production needed it because `cactus_78.vcf.gz` came from a different source (`cactus_pang69_1001g.vcf.gz`) than the 135-asm annotation catalog. Here we subset the annotated 135-asm VCF directly, so `INFO/ID` is already present.

### PanGenie-index: pang_135 (production), not p82's

p80 uses the production **`pang_135_pangenie_index_Chr1_kmers.tsv.gz`** for kmer_pa, NOT `control_p82/data/kmers_p82/`. The reason is graph identity:

- The arch3 A1 catalog used here was annotated against `pang_1001gplus_all.gfa.gz` (135-asm production graph). Our p80 VCF inherits that graph's bubble structure.
- `kmers_p82/` was built from the standalone 82-acc minigraph-cactus run — a different cactus invocation with different bubble topology and k-mer dictionary.

Mixing them would have `build_kmer_pa.py` reconstruct haplotypes from 135-asm records inside 82-acc-defined bubble regions — a silent inconsistency that exactly mirrors the v3 sim-FASTA bug (`memory/project_v3_singleton_kmer_bug`). The pang_135 PG-index is the topology-matching one; production `kmer_pa` is built on the same pang_135 graph (but from the in-house index, not the PG-index — `docs/PIPELINE_STATE.md` §0). This benchmark uses the pang_135 PG-index, which shares that graph topology.

## Auditor caveats to call out in any writeup

1. **Numbers from `../HANDOFF_P82_CONTROL.md` do NOT carry over.** The p82 vs v3 apples-to-apples table (775,944 shared SNPs; MAE/RMSE/R²) was tied to the p82 build. p80 differs on both axes (sample drop AND decomposition method). Re-run every number.
2. **5772 drop is a mixed effect.** With 5772 + T980 both present in p82, EM mass-split across the duplicate columns (cactus-FASTA-mislabel artifact). Dropping 5772 *also* removes that simplex degeneracy. So any p80 vs p82 lift mixes (a) drop-noisy-9947, (b) remove-5772-T980-duplicate, (c) arch decomposition. Don't conflate.
3. **Decomposition is different.** The arch path uses graph-topology symbolic IDs (`convert-to-biallelic.py`) rather than reference-aligned `bcftools norm -m -any`. The variant ID space and exact (chrom, pos, ref, alt) atomization may differ for MNPs and complex bubbles. Intersect on (chrom, pos, ref, alt) when comparing against p82 or v3.
4. **PanGenie-index is `pang_135` (production), not p82's.** Required for graph-topology consistency with the arch3 A1 catalog (both come from `pang_1001gplus_all`). k-mers private to 5772/9947 in pang_135 land as all-zero kmer_pa columns after p80 sample-subset; not wrong, just inefficient. If pruning is needed later, drop columns where `kmer_pa.sum(axis=0) == 0`.

## Layout

```
benchmarks/p80/
├── README.md                                  # this file
├── data/                                      # all artifacts derive from one source VCF
│   ├── samples_80_asm.txt                           # 80 Asm_IDs to keep (= 82-acc − exclude_list)
│   ├── samples_80_rename.tsv                        # Assembly_ID -> Accession_ID (80 rows)
│   ├── samples_80_acc_order.txt                     # canonical 80-Acc_ID order (post-rename)
│   ├── pangenome_p80_chr1.vcf.gz                    # CANONICAL biallelic VCF (80 samples, Acc_ID)
│   ├── pangenome_p80_chr1.vcf.gz.tbi
│   ├── kmer_pa_p80/
│   │   ├── kmer_pa_Chr1.kmer_pa.npz                           # founder x k-mer (80 x ~6M)
│   │   └── kmer_pa_Chr1.meta.npz
│   ├── var_pa_p80.var_pa.npz                        # founder x variant
│   └── var_pa_p80.meta.npz
├── fastas_80/                                 # bcftools-consensus founder FASTAs (Chr1)
│   └── <Accession_ID>.chr.fa[.fai]            # 80 of these
├── sims/                                      # VISOR sim dirs (one per regime)
│   ├── cov10_n50_g1_s42_hotspots_p80_chr1/
│   ├── cov10_n50_g3_s42_hotspots_p80_chr1/
│   └── ...
├── results/
│   ├── cactus_em_global/<regime>/<sample>.tsv
│   ├── cactus_em_star2/<regime>/<sample>.tsv
│   ├── v3_panel_on_p80_reads/...              # for v3-vs-p80 apples-to-apples (optional)
│   └── FINAL_RESULTS_cov10_p80.ipynb          # TODO (Phase D)
├── scripts/
│   ├── 01_build_biallelic_p80.sh              # subset 135-asm annotated to 80 + convert-to-biallelic + fill-tags + reheader
│   ├── 03_build_kmer_pa_p80.sh                # kmer_pa_p80 (uses the pang_135 PG-index, NOT control_p82's)
│   ├── 04_build_var_pa_p80.sh                 # var_pa_p80
│   ├── 05_build_fastas_p80.sh                 # 80-founder consensus FASTAs (array 1-80)
│   ├── 06_run_sim_p80.sh                      # recomb sim per regime
│   ├── 07_run_cactus_em_p80.sh                # cactus_em per (regime, method)
│   ├── 08_run_v3panel_on_p80reads.sh          # v3 panel on p80 reads (optional)
│   ├── make_recomb_mosaics_p80.py             # mosaic builder
│   └── submit_all_p80.sh                      # chained-dependency submission
└── logs/                                      # SLURM stdout/stderr
```

## How to run

```bash
cd /carnegie/nobackup/scratch/tbellagio/kmate/benchmarks/p80
bash scripts/submit_all_p80.sh
```

Or step-by-step:

```bash
# Phase A: panel artifacts (arch decomposition; A1 already built in panel/arch3/chr1/)
A1=$(sbatch --parsable scripts/01_build_biallelic_p80.sh)
A3=$(sbatch --dependency=afterok:$A1 --parsable scripts/03_build_kmer_pa_p80.sh)
A4=$(sbatch --dependency=afterok:$A1 --parsable scripts/04_build_var_pa_p80.sh)
A5=$(sbatch --dependency=afterok:$A1 --parsable scripts/05_build_fastas_p80.sh)

# Phase B: sims
B1=$(sbatch --dependency=afterok:$A3:$A4:$A5 --parsable scripts/06_run_sim_p80.sh 50 1)
B3=$(sbatch --dependency=afterok:$A3:$A4:$A5 --parsable scripts/06_run_sim_p80.sh 50 3)

# Phase C: methods
sbatch --dependency=afterok:$B1 scripts/07_run_cactus_em_p80.sh n50_g1 global
sbatch --dependency=afterok:$B1 scripts/07_run_cactus_em_p80.sh n50_g1 star2
sbatch --dependency=afterok:$B3 scripts/07_run_cactus_em_p80.sh n50_g3 global
sbatch --dependency=afterok:$B3 scripts/07_run_cactus_em_p80.sh n50_g3 star2
```

## External tools / paths

- **convert-to-biallelic**: `../external/pangenie-tools/pipelines/run-from-callset/scripts/convert-to-biallelic.py`
- **A1 shared catalog (annotated 135-asm Chr1)**: `../panel/arch3/chr1/chr1_135_annotated.sorted.vcf.gz`
- **A1 shared biallelic catalog**: `../panel/arch3/chr1/chr1_135_annotated_biallelic.sorted.vcf.gz`
- **82-acc source VCF (Asm_ID list reference)**: `/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.vcf.gz`
- **Reference FASTA**: `/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa`
- **PG-index (reused)**: `../panel/pangenie_genotyping/data/pang_135_pangenie_index_Chr1_kmers.tsv.gz` (production, 135-asm)

## Consistency invariant

Every p80-specific artifact (`kmer_pa_p80`, `var_pa_p80`, `fastas_80/`) derives from the **same canonical biallelic VCF** (`data/pangenome_p80_chr1.vcf.gz`). Do not symlink to p82's FASTAs or to v3 artifacts — variant-set mismatch would recreate the v3 simulation bug (`memory/project_v3_singleton_kmer_bug`).

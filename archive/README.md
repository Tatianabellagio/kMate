# archive/ — superseded methods, lost benchmarks, and the original deliverable

**Created 2026-05-07.** Things in here are kept for reproducibility / postmortem reference but are NOT part of the current production pipeline. None of these directories are imported or read by any active script (verified by grep before the move).

If you're wondering whether something was tried before, check `../METHODS_TRIED.md` first. This README explains *what* lives here and *why each piece is here* (rather than at the top level).

---

## Directory map

| Path | What | Why archived | Postmortem |
|---|---|---|---|
| `route3_kallisto_em/` | kallisto-style read-level pseudoalignment + equivalence-class EM | LOST — 6 pp WORSE SNP R² than `window_10kb` per-k-mer Poisson at the same resolution | `RESULTS_LOG.md` 2026-05-06 "kallisto-EM" entry; `METHODS_TRIED.md` §2 Route 3 |
| `bench_kmc/` | KMC 3 + rapidgzip vs jellyfish k-mer-counting bench | NO WIN — KMC: +80% memory regression. rapidgzip: ~7% speedup once the per-file-CLI bug was fixed | `HANDOFF.md` 2026-05-02 entry "k-mer counting alternatives benchmarked" |
| `hapfire_projection/` | The original Chr1 hapFIRE-projection deliverable | SUPERSEDED by the cactus_em end-to-end recipe (locked in 2026-05-07) | `SUMMARY.md` (kept at root) for the original writeup |

---

## `route3_kallisto_em/`

The architecture (k-mer hash → per-read compatibility-set intersection → equivalence-class EM) is canonical from kallisto (Bray 2016). We applied it with founders-as-transcripts. Implementation worked, scaled (40 min/sample on 1 core, 261 K unique ECs at full cov50), but produced R²=0.842 SNP at 10 kb on n50_g3 vs `window_10kb`'s 0.902. The architectural reason: EC-multinomial loses the absence signal that per-k-mer Poisson uses (k-mers expected to be present but observed at low count constrain `h^T cn[:,k]` downward; ECs only see observed reads).

| File | Provenance |
|---|---|
| `per_sample_kallisto_em.py` | from `poolfreq/src/` |
| `precompute_panel_u64_index.py` | from `poolfreq/scripts/` |
| `eval_kallisto_em.py` | from `poolfreq/scripts/` |
| `append_kallisto_to_session_dataset.py` | from `poolfreq/scripts/` (appended kallisto runs to `af_long`/`af_summary`) |
| `run_anchor_sweep.sh` | from `poolfreq/scripts/` (anchor-weight sweep launcher for kallisto-EM) |
| `run_kallisto_other_regimes.sh` | from `sims/visor_freqk/scripts/` (Route 3 launcher for n200_g1 / n50_g1 / skewed) |

**When to revisit.** If we need to push to whole-genome 5-chrom 8-kb-window scale on cluster-wide concurrency and per-iter wall is the bottleneck, kallisto-style EC EM scales O(unique_ec) per iter (~10⁵-10⁶) instead of O(K)~10⁷. Possible rescue paths in `METHODS_TRIED.md` §6 (Mora set-cover regularization, Route 1 anchor stacked on EC EM, KAGE 2 variant-correlation prior).

---

## `bench_kmc/`

36 GB total (mostly `rgz_59320/`, `rgz2_59551/`, `work_59298/` — rapidgzip and KMC intermediate output). Could be deleted entirely; everything is reproducible from the scripts.

| File | Notes |
|---|---|
| `run_bench.sh` | KMC 3 vs jellyfish bench |
| `run_rapidgzip_bench.sh` / `run_rapidgzip_bench2.sh` | rapidgzip vs zcat front-end |
| `compare_outputs.py` | Output identity check (jellyfish reference) |
| `build_query_fa.py`, `query_chr1.fa` | Query input |
| `bench_*.{out,err}`, `rgz_bench*.{out,err}` | SLURM logs |
| `work_59298/` | KMC intermediate |
| `rgz_59320/`, `rgz2_59551/` | rapidgzip intermediate |

**Conclusion (HANDOFF.md 2026-05-02).** jellyfish at 8 threads is already saturating CPU (618%) with zcat as front-end. Gzip decompression only bottlenecks below ~4 threads. The next real lever is SSHash or LP-MPHF, ~1-2 days to integrate; deferred until per-chrom + memex concurrency proves insufficient.

---

## `hapfire_projection/`

The original deliverable that the project pivoted away from. Method: run hapFIRE on SNPs, project the 231-vector ecotype-frequency output onto SVs via a precomputed founder × SV genotype matrix (`f_SV(v) = Σ_e G_SV(e, v) · f_ecotype(e)`). Worked very well on the no-recomb sim (MAE 0.0006-0.0013 across cov{10,20,50}, 75-100× better than freqk) and on the 8-sample SEEDMIX (r=0.99 vs panel-truth). But:

1. The no-recomb sim is best-case for hapFIRE's intact-founder assumption — recombination violates it.
2. Panel ascertainment ceiling: only 80 of 231 founders had assemblies → SVs polymorphic only among the un-genotyped 151 are invisible.
3. Throughput: ~3 h/sample for hapFIRE end-to-end (HARP serial across coarse blocks), vs ~25 min/sample for the cactus_em recipe.

The pivot to cactus_em (2026-04-27 onward) addressed all three.

### Layout

```
hapfire_projection/
├── SIMULATION_RESULTS.ipynb     # original no-recomb sim notebook
├── scripts/                     # 12 files: builders, projectors, hapFIRE wrappers
│   ├── build_founder_sv_matrix.py     # SV panel VCF → 50K×82 0/1 dosage
│   ├── project_seedmix_svs.py         # F @ G.T projection + freqk comparison
│   ├── seedmix_truth_analysis.py      # recipe-based panel-truth
│   ├── plot_seedmix.py / plot_simulation.py / posterior_combine.py
│   ├── project_sv_freq.py             # generic per-sample projector
│   ├── compare_methods.py / diagnostic_freqk_vs_hapfire_disagreement.py
│   ├── run_hapfire_one.sh / run_all_analysis.sh
│   └── test_projection_truth.py       # math sanity check
├── data/
│   ├── founder_sv_matrix.parquet      # 50,446 SVs × 82 ecotypes, 0/1 dosage
│   ├── founder_sv_meta.parquet        # chrom/pos/alt_idx/var_type/sv_size/AC
│   ├── seedmix_recipe.tsv             # un-normalized recipe (the normalized version
│   │                                    `seedmix_recipe_normalized.tsv` STAYS at
│   │                                    `../../data/` — used by imputation/)
│   └── vcf/                           # greneNet_Chr1_only.vcf (798 MB) +
│                                        greneNet_final_v1.1.recode.vcf.gz (92 MB)
└── results/
    ├── rep{1,9,10,11}_cov{10,20,50}_1kb_f{10,30,50,70,90}/    # 25 hapFIRE outputs
    ├── simulation_summary.tsv / simulation_combined.tsv      # aggregated truth-vs-method
    ├── simulation_truth_vs_estimate_by_cov.png               # key figure
    ├── simulation_mae_by_cov.png                             # key figure
    ├── seedmix_3way_comparison.tsv.gz                        # 8-sample × SV table
    ├── seedmix_S1_comparison.png / seedmix_correlation_bars.png
    ├── seedmix_per_sample_summary.tsv
    ├── hapfire_proj_svs_seedmix.parquet                      # projection output
    ├── hapfire_vs_freqk_svs_seedmix_AF.csv.gz
    ├── hapfire_vs_freqk_svs_seedmix_per_sample_corr.tsv
    └── sanity_check_truth_projection.tsv
```

### What stayed at the root

These files are still used by the **active** imputation pipeline at `../../imputation/`:

- `../../data/sv_panel_to_accession_id.tsv` (Assembly_ID → 1001G ID mapping)
- `../../data/vcf_samples_231.txt` (231-founder VCF header order)
- `../../data/seedmix_recipe_normalized.tsv` (sum=1 recipe — used by imputation validation)

`../../SUMMARY.md` (the original writeup) is also kept at the root; it now opens with a status note pointing to FINAL_RESULTS.

### When to revisit

- If a real evolved-sample test ever shows hapFIRE-proj outperforming cactus_em at sites where `cn_var_231_v2` ascertainment is poor (the 151 imputed founders), the projection scripts in `scripts/` are ready to re-run on the existing xwu hapFIRE outputs.
- The posterior-combine framework (`scripts/posterior_combine.py`) is the half-day version of the "LD-empowered k-mer" idea (HANDOFF "open next steps" #4) — could plug into a freqk + cactus_em combination if we ever want one.

---

## How to restore something from archive

```bash
# Single file (untracked):
mv archive/route3_kallisto_em/per_sample_kallisto_em.py poolfreq/src/

# Tracked file (use git mv to keep history):
git mv archive/hapfire_projection/scripts/posterior_combine.py scripts/

# Restore an output dir for re-analysis:
mv archive/hapfire_projection/results/rep9_cov50_1kb_f30 results/
```

If you restore something, update `METHODS_TRIED.md` and remove the corresponding row from this README.

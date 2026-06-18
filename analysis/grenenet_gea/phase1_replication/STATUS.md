# Status — phase-1 GEA replication (updated 2026-06-17)

Re-running the GrENE-Net phase-1 SNP GEA pipeline on kMate AF time-series, for
three variant classes, to reproduce the CAM5 headline and ask whether SVs/indels
carry signal SNPs miss. Original pipeline: `/global/scratch/users/tbellg/gea_grene-net/`.

This folder is **self-contained**: it holds ONLY the phase-1 replication (3 models ×
last-gen × WZA → datasets → figure). The PC1/PC1-LFMM robustness branch lives in
`../wza_investigation/`; one-block deep-dives live in `followups/`.

## ★ FINAL DELIVERABLES ★ (all under `results/grenenet_gea/phase1_replication/`)
- **gene table (per gene):** `gene_summary_gen9_bio1_deg7cap2000.csv` — one row/gene:
  symbol, function, region (Chr:start-end), categories, n_models, best_q, q_floored.
- **gene table (per block,gene):** `significant_genes_annotated_gen9_bio1_deg7cap2000.csv`
  (+ `significant_genes_…csv` without function; `gene_function_…csv`; mapping cross-check
  `gene_mapping_verification_…csv`).
- **model comparison:** `compare_all3models_gen9_bio1_deg7.csv`.
- **figure:** `manhattan_3models_deg7cap2000.png` ← notebook
  `notebooks/manhattan_3models_deg7cap2000.ipynb` (the 3×3 class×model Manhattan panel).
- **per-record + block stats:** `{kendall,lfmm,binomial}/`, `wza/`, inputs in `class_matrices/`.

## The branches (dimensions)

| dimension | values | notes |
|---|---|---|
| **variant class** | snp · smallindel · sv | filter: MAF≥0.05 + finite in ≥50% pools; every kept record carries an LD block |
| **generation** | gen1 · gen3 · **gen9** | gen9 = "last-gen" (flower-weighted timepoint merge) = the phase-1 analysis unit |
| **model** | kendall · lfmm · binomial | three independent per-record GEA tests |
| **WZA correction** | **deg7-cap2000 (PRIMARY)** · deg7-nocap · deg2nocap | deg-7 fits the SNP-count→SD support curve best (RMSE 0.94 vs deg-2 2.06); cap removes tail-extrapolation. See `../wza_investigation/RESULTS.md` |
| **climate** | bio1 | bio2–19 not yet run |
| **honest null** | PC1 / PC1-after-LFMM + site-permutation | separate branch (`../wza_investigation/`, `run_block_pc1*`) |

### Three models
- **kendall** (`run_kendall.py`) — Kendall-τ between per-pool AF and per-pool bio1.
  Raw, structure-uncorrected. Run for all 3 gens.
- **lfmm** (`run_lfmm_gea.sh` → `run_lfmm_lastgen.R`, K=16) — latent-factor mixed model,
  absorbs site structure as latent factors. gen9 only (the phase-1 published model).
- **binomial** (`run_binomial.py`) — per-SNP binomial GLM `[alt,ref] ~ const+z(bio1)`,
  allele counts = **AF×total_flowers×2** (diploid genomes; phase-1's flower-count scheme,
  from `create_partitions_allele_freq.py`). gen9 only. Env-coef Wald p → WZA.

All three emit `{chrom,pos,ref_len,alt_len,MAF,block,(stat),pval}` so `run_wza.py`
consumes them identically.

## Coverage (what exists)

Per-record model outputs (`results/.../{kendall,lfmm,binomial}/`):

| model | gen1 | gen3 | gen9 |
|---|---|---|---|
| kendall  | ✓ | ✓ | ✓ |
| lfmm     | — | — | ✓ |
| binomial | — | — | ✓ |

WZA block-level (`results/.../wza/`, gen9 last-gen comparison is the complete one):

| run | deg2nocap | deg7nocap | deg7cap2000 |
|---|---|---|---|
| kendall gen1/gen3 | ✓ | — | — |
| kendall gen9      | snp only | ✓ | ✓ |
| lfmm gen9         | — | ✓ | ✓ |
| binomial gen9     | ✓ | ✓ | ✓ |

## ★ Result: CAM5 (block 2_1265) reproduces in all 3 models, SNP gen9, deg-7 ★
(`compare_all3models_gen9_bio1_deg7.csv`; CAM5 honest poly-free p ≈ 1.6e-4 from wza_investigation §4)

| model | regime | CAM5 rank /16,447 | p | q | BH-sig | tot BH |
|---|---|---|---|---|---|---|
| binomial | deg7cap2000 | 33 | **2.0e-4** | 0.092 | top-0.2% | 9 |
| kendall  | deg7cap2000 | 5–8 | 8.5e-6 | 0.017 | ✅ | 19 |
| lfmm     | deg7cap2000 | 15  | 1.5e-5 | 0.017 | ✅ | 28 |

- **binomial's deg-7 CAM5 p (2.0e-4) lands on the poly-free ground truth (1.6e-4)** — the
  best-calibrated of the three; kendall/lfmm are more anti-conservative.
- deg-2 → deg-7 moves binomial CAM5 rank 169 → 33 (deg-7 lifts CAM5-sized ~13-SNP windows).
- Binomial raw p<0.05 ~84% (pool over-precision: N=flowers×2 genomes treated independent),
  but WZA's relative SNP-number correction absorbs the ~uniform inflation → *fewer* BH hits
  (9) than kendall (19)/lfmm (28). cap vs nocap barely moves anything under deg-7.
- **SV is null at CAM5 in all three models** (CAM5 was a SNP signal; only block 2_1264 — a
  null arm block — overlaps the gene in the SV set).
- Cross-class recurrence (deg7): kendall shares **4_2519** (Chr4 CRK defense cluster) across
  all 3 classes — the dominant raw peak, structure-suspect (see wza_investigation).

## Key scripts (all in this folder)
- `build_class_matrices.py` — split/filter pool matrices → per-class/gen `{af.npy, records.csv, pools.csv}`.
- `run_kendall.py` / `run_lfmm_gea.sh`(+`run_lfmm_lastgen.R`) / `run_binomial.py`(+`.sh`) — the 3 models.
- `run_wza.py` — WZA wrapper; `--poly-deg/--sample-snps/--min-entries/--regime` for any correction.
  Default = canonical deg-2; pass `--poly-deg 7 --min-entries 10 --sample-snps 2000 --regime deg7cap2000`.
- `run_wza_binomial.sh` — SLURM driver for binomial WZA.
- `compare_blocks.py` — BH-FDR/Bonferroni per output, cross-class/gen recurrence, CAM5 rank
  in every WZA output. `--deg deg7|deg2|both`. Output: `compare_all3models_*`.
- `build_significant_genes.py` — BH-sig blocks (all 9 combos) → genomic span → overlapping
  TAIR10 genes (local GFF) → symbol+description via Ensembl Plants REST API → gene×combo
  table. Output: `significant_genes_gen9_bio1_deg7cap2000.csv` (93 blocks, 220 genes;
  CRK cluster block 4_2519 = 7/9 combos, CAM5 = 3, SRO2 = 4).
- `verify_gene_mapping.py` — independent cross-check of block→gene overlap vs Ensembl
  `overlap/region` (Araport11). 219/222 (98.6%) confirmed; ~1.4% = TAIR10-vs-Araport11
  annotation version. Output: `gene_mapping_verification_*.csv`.
- `annotate_gene_function.py` — per-gene NCBI name+summary+GO via mygene.info, curated
  climate/stress/flowering classification (evidence recorded, 67/77 GO-backed). Outputs:
  `gene_function_*.csv` + `significant_genes_annotated_*.csv` (sig table + function). 58/220
  climate/stress/flowering (CRK cluster = defense/oxidative; CAM5 = calcium).
- `build_gene_summary.py` — GENE-level table (1 row/gene): function + how many of the 3
  models/classes/9 combos flagged it. Output: `gene_summary_*.csv` (17 genes hit by all 3
  models = CRK cluster + AT1G33040–33080; CAM5 = 2 models).
- `_build_manhattan_3x3_nb.py` → `notebooks/manhattan_3models_deg7cap2000.ipynb` — the 3×3
  (class×model) last-gen WZA Manhattan grid + `manhattan_3models_deg7cap2000.png`.
- `followups/crk_block_af_support.py` — deep-dive on the top block (CRK cluster 4_2519):
  raw AF vs climate (up-in-cold, site τ=−0.42) + kMate panel support (n_called 226/231) + error.
- `logs/` — archived SLURM `.out` job logs.
- **PC1 / PC1-LFMM robustness branch MOVED to `../wza_investigation/`** (`run_block_pc1*`,
  `build_block_pc1_matrix.py`, `manhattan_pc1.ipynb`, results in `wza_investigation/pc1/`):
  it's the LD-collapse + structure-correction check, not the phase-1 method.

## Environments
- **statsmodels (binomial)** lives in the **`basic`** env, NOT `kmate`. WZA/compare use `kmate`.
- Always **sbatch to a clean node** (savio4_htc / co_moilab / savio_lowprio). The OOD Jupyter
  node stalls in I/O (D-state) on big-npy loads.

## NEXT (open)
1. Extend binomial + lfmm to **gen1/gen3** for full Kendall-parity (currently gen9 only).
2. **bio2–19** — only bio1 run so far (the whole climate dimension is open).
3. Optional: deg-7 WZA for kendall gen1/gen3 (currently only deg-2 there).
4. Fold the honest PC1/site-permutation null into the cross-model story.

## Target
Phase-1 headline **CAM5 / AT2G27030, Chr2 ~11.53 Mb** = LD blocks 2_1264 (null, gene midpoint)
+ **2_1265** (the real signal, gene 3' end, ~13 SNPs).

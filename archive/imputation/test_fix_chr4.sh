#!/bin/bash
#SBATCH --job-name=test_fix_chr4
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/test_fix_chr4_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/test_fix_chr4_%j.err

# =============================================================================
# Test the proper fix for the 1.43× slope bias.
#
# Bug: cn_var_231 was built from `cactus_svs + imputed_151`. cactus_svs has
#      no SNP records, so 80 cactus founders end up with ./. (=> 0) at all
#      3.24M SNP records. EM compensates by over-attributing to imputed
#      founders. Predicted alt_freq scales 1.43× over recipe truth.
#
# Fix:  Build from `ref_80 + imputed_151` instead. ref_80 has the 80
#       cactus founders' SNP genotypes (from grene_80 used as Beagle's
#       reference panel input).
#
# Test on Chr4 (smallest chromosome) end-to-end:
#   1. bcftools merge ref_80 + imputed_151, restrict to chrom 4
#   2. Rename chrom 4 → Chr4
#   3. Build cn_kmer_v2 + cn_var_v2 for Chr4
#   4. Run per_sample_driver on SEEDMIX_S1 reads, Chr4-only
#   5. Compute slope: should be ~1.0 if fix works
# =============================================================================
set -euo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work

BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix
PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python
POOLFREQ=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
mkdir -p test_fix

echo "[$(date)] === Step 1: merge ref_80 + imputed_151 → corrected 231-founder VCF (Chr4 only) ==="
$BCF merge ref_80.vcf.gz imputed_151.vcf.gz -r 4 \
    -Oz -o test_fix/merged_v2_chr4.vcf.gz --threads 8
$TABIX -p vcf test_fix/merged_v2_chr4.vcf.gz
echo "  records: $($BCF view -H test_fix/merged_v2_chr4.vcf.gz | wc -l)"
echo "  samples: $($BCF query -l test_fix/merged_v2_chr4.vcf.gz | wc -l)"

echo ""
echo "[$(date)] === Step 2: rename chrom 4 → Chr4 ==="
echo -e "4\tChr4" > test_fix/rename_chr4.txt
$BCF annotate --rename-chrs test_fix/rename_chr4.txt test_fix/merged_v2_chr4.vcf.gz \
    -Oz -o test_fix/merged_v2_Chr4.vcf.gz --threads 8
$TABIX -p vcf test_fix/merged_v2_Chr4.vcf.gz

echo ""
echo "[$(date)] === Step 3a: build cn_var_v2 for Chr4 ==="
$PYTHON ${POOLFREQ}/src/build_cn_var.py \
    --vcf test_fix/merged_v2_Chr4.vcf.gz \
    --out test_fix/cn_var_v2_Chr4

echo ""
echo "[$(date)] === Step 3b: build cn_kmer_v2 for Chr4 (this is the slow step ~1.5h) ==="
$PYTHON ${POOLFREQ}/src/build_kmer_cn.py \
    --kmers /carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_test/pangenie_idx_rawv2_Chr4_kmers.tsv.gz \
    --vcf   test_fix/merged_v2_Chr4.vcf.gz \
    --ref   /home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa \
    --chrom Chr4 \
    --out   test_fix/cn_v2_Chr4

echo ""
echo "[$(date)] === Step 4: density check on cn_var_v2 vs original cn_var_231 (Chr4 records only) ==="
$PYTHON -u <<'PYEOF'
import numpy as np, pandas as pd
from scipy.sparse import load_npz

# Original cn_var_231, restricted to Chr4
cn_var_old = load_npz("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231.cn_var.npz")
meta_old = np.load("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231.meta.npz", allow_pickle=True)
chr4_old = meta_old["chrom"] == "Chr4"
cv_old = cn_var_old[:, chr4_old]
founders_old = list(meta_old["founders"])

# New cn_var_v2 for Chr4
cn_var_new = load_npz("test_fix/cn_var_v2_Chr4.cn_var.npz")
meta_new = np.load("test_fix/cn_var_v2_Chr4.meta.npz", allow_pickle=True)
founders_new = list(meta_new["founders"])

# Identify cactus vs imputed founders
panel_map = pd.read_csv("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/sv_panel_to_accession_id.tsv", sep="\t")
cactus_set = set(panel_map.Accession_ID.astype(str))
in_cactus_old = np.array([f in cactus_set for f in founders_old])
in_cactus_new = np.array([f in cactus_set for f in founders_new])

print(f"OLD cn_var_231 (Chr4): shape {cv_old.shape}, cactus founders: {in_cactus_old.sum()}, imputed: {(~in_cactus_old).sum()}")
print(f"NEW cn_var_v2_Chr4:    shape {cn_var_new.shape}, cactus founders: {in_cactus_new.sum()}, imputed: {(~in_cactus_new).sum()}")

den_old = np.asarray(cv_old.mean(axis=1)).flatten()
den_new = np.asarray(cn_var_new.mean(axis=1)).flatten()
print(f"\nPer-founder density (Chr4 records):")
print(f"  OLD cactus  mean: {den_old[in_cactus_old].mean():.4f} (was: ~0)")
print(f"  OLD imputed mean: {den_old[~in_cactus_old].mean():.4f}")
print(f"  NEW cactus  mean: {den_new[in_cactus_new].mean():.4f}  ← should match imputed if fix works")
print(f"  NEW imputed mean: {den_new[~in_cactus_new].mean():.4f}")
print(f"  ratio NEW imputed/cactus: {den_new[~in_cactus_new].mean()/den_new[in_cactus_new].mean():.2f}× (was ~50× in OLD)")
PYEOF

echo ""
echo "[$(date)] === Step 5: run per_sample_driver on SEEDMIX_S1 reads, Chr4 cn ==="
mkdir -p test_fix/results
$PYTHON -u ${POOLFREQ}/src/per_sample_driver.py \
    --cn-kmer-prefix test_fix/cn_v2 \
    --cn-var test_fix/cn_var_v2_Chr4.cn_var.npz \
    --cn-var-meta test_fix/cn_var_v2_Chr4.meta.npz \
    --reads /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz \
            /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz \
    --sample SEEDMIX_S1_v2 \
    --out test_fix/results/SEEDMIX_S1_v2.tsv \
    --threads 8 \
    --block-mode global

echo ""
echo "[$(date)] === Step 6: slope diagnostic ==="
$PYTHON -u <<'PYEOF'
import numpy as np, pandas as pd
from scipy.sparse import load_npz

cn_var_v2 = load_npz("test_fix/cn_var_v2_Chr4.cn_var.npz")
meta_v2 = np.load("test_fix/cn_var_v2_Chr4.meta.npz", allow_pickle=True)
founders = list(meta_v2["founders"])
recipe = pd.read_csv("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/seedmix_recipe_normalized.tsv", sep="\t")
rd = dict(zip(recipe.ID.astype(str), recipe.seed_prop))
h_truth = np.array([rd.get(str(f), 0.0) for f in founders])
print(f"Recipe panel mass on cn_var_v2 founders: {h_truth.sum():.4f}")
h_truth = h_truth / h_truth.sum()
af_truth = (h_truth.astype(np.float32) @ cn_var_v2.toarray()).astype(np.float64)

df = pd.read_csv("test_fix/results/SEEDMIX_S1_v2.tsv", sep="\t")
af_pred = df["alt_freq"].to_numpy()
print(f"af_truth: n={len(af_truth):,}, mean={af_truth.mean():.4f}")
print(f"af_pred:  n={len(af_pred):,}, mean={af_pred.mean():.4f}")
slope, intercept = np.polyfit(af_truth, af_pred, 1)
r = np.corrcoef(af_pred, af_truth)[0,1]
r2 = 1 - ((af_pred - af_truth)**2).sum() / ((af_truth - af_truth.mean())**2).sum()
rmse = float(np.sqrt(np.mean((af_pred - af_truth)**2)))
print(f"\n=== HEADLINE: slope = {slope:.4f}  (was 1.43 with the bug)")
print(f"  intercept = {intercept:.5f}")
print(f"  R² (raw)  = {r2:.4f}  (was 0.68)")
print(f"  Pearson r = {r:.4f}")
print(f"  RMSE      = {rmse:.4f}")
print(f"\n  → if slope ≈ 1.00, the root-cause fix works → trigger full rebuild")
print(f"  → if slope still > 1.2, something else is contributing")
PYEOF
echo ""
echo "[$(date)] DONE"

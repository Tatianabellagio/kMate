#!/bin/bash
#SBATCH --job-name=loo_imp_concord
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/loo_imp_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/loo_imp_%A_%a.err

# =============================================================================
# loo_concord_imputed_one.sh — re-run LOO concordance on imputed merged VCF
# as the truth source. Compares PanGenie LOO short-read call (per-sample VCF)
# against the imputed 231-founder VCF (where cactus-side `.` cells have been
# Beagle-imputed). Tells us whether imputation improved or degraded the
# LOO concordance vs the original cactus assembly truth.
#
# Output: data/loo_concordance_imputed/<eco>_summary.tsv (same schema as the
# pre-imputation loo_concordance/ TSVs).
# =============================================================================
set -eo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping
MANIFEST=$BASE/data/loo_ena_manifest.tsv
GT_DIR=$BASE/data/loo_genotyped
TRUTH_VCF=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work_merged/founders_231_imputed_multiallelic.vcf.gz
OUT_DIR=$BASE/data/loo_concordance_imputed
mkdir -p $OUT_DIR

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $MANIFEST)
[ -n "$LINE" ] || { echo "ERROR: no row at $IDX" >&2; exit 1; }
ECOTYPE=$(echo "$LINE" | cut -f1)

PG_VCF=$GT_DIR/${ECOTYPE}_genotyping.vcf.gz
if [ ! -s "$PG_VCF" ]; then
    echo "ERROR: PanGenie LOO VCF missing for $ECOTYPE — skipping (pre-imputation runs also failed for the 3 truncated R2 ecotypes)" >&2
    exit 0
fi
if [ -f "$OUT_DIR/${ECOTYPE}_summary.tsv" ]; then
    echo "[$(date)] $ECOTYPE: already done"; exit 0
fi

# In the imputed merged VCF, sample names are 1001G ecotype IDs (renamed in
# merge_vcfs.sh Step A). So truth-sample == sample.
PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python
echo "[$(date)] $ECOTYPE: loo_concordance vs imputed truth"
$PYTHON $BASE/scripts/loo_concordance.py \
    --pangenie-vcf $PG_VCF \
    --truth-vcf    $TRUTH_VCF \
    --sample       $ECOTYPE \
    --truth-sample $ECOTYPE \
    --out          $OUT_DIR/${ECOTYPE}

echo "[$(date)] DONE"
ls -lh $OUT_DIR/${ECOTYPE}_summary.tsv

#!/bin/bash
#SBATCH --job-name=pg_loo
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --output=logs/loo_%A_%a.out
#SBATCH --error=logs/loo_%A_%a.err

# =============================================================================
# validate_loo.sh
# Stage 5: PanGenie LOO concordance check on a held-out cactus founder.
#
# For one cactus accession with both long-read assembly truth AND short reads
# (we'd need to identify which cactus founders also have public ENA reads),
# re-genotype it with PanGenie using its short reads, then compare to its
# cactus-derived ground truth genotypes.
#
# Target: ≥97% concordance, matching the 98.5% Beagle LOO benchmark.
#
# Indexed by SLURM_ARRAY_TASK_ID over a list of LOO targets.
# =============================================================================
mkdir -p logs
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
eval "$(conda shell.bash hook)"
conda activate pangenie

BASE=/global/scratch/users/tbellg/hapfire_sv/pangenie_genotyping

# Ecotypes to LOO-test: cactus founders that ALSO have ENA short reads.
# (Their cactus assembly = ground truth; their short reads = test input)
LOO_LIST=$BASE/data/loo_targets.txt   # e.g. 9542, 7186 etc.
mkdir -p $BASE/data/loo

IDX=${SLURM_ARRAY_TASK_ID:-1}
ECOTYPE=$(sed -n "${IDX}p" $LOO_LIST)
if [ -z "$ECOTYPE" ]; then echo "ERROR: no ecotype at $IDX"; exit 1; fi
echo "[$(date)] LOO test for cactus founder $ECOTYPE"

# Run PanGenie genotype using its short reads
PREP_R1=$BASE/data/preprocessed/${ECOTYPE}_1.dedup.fq.gz
PREP_R2=$BASE/data/preprocessed/${ECOTYPE}_2.dedup.fq.gz
PANG69_VCF=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa

OUT_VCF=$BASE/data/loo/${ECOTYPE}_pangenie.vcf
PanGenie -i "$PREP_R1 $PREP_R2" -r $REF -v $PANG69_VCF -o $OUT_VCF -s $ECOTYPE -t 8 -j 8

# Compare to cactus-derived truth
BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
TRUTH_VCF=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz
echo "[$(date)] concordance check"

bgzip -f $OUT_VCF; tabix -p vcf ${OUT_VCF}.gz
$BCF stats -s $ECOTYPE ${OUT_VCF}.gz $TRUTH_VCF \
    > $BASE/data/loo/${ECOTYPE}_concordance.txt

# Pull headline: % concordant
grep -E "non-Reference Discordance|nRDs" $BASE/data/loo/${ECOTYPE}_concordance.txt | head
echo "[$(date)] DONE — full report at $BASE/data/loo/${ECOTYPE}_concordance.txt"

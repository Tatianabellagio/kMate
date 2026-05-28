#!/bin/bash
# =============================================================================
# launch_loo_pangenie.sh
# Submit a SLURM array of LOO PanGenie genotype jobs (genotype + concordance).
# Pass an array spec as $1 to restrict (e.g. "1-5" for the smallest 5 LOO
# candidates). Default: all 78 manifest rows, 8 concurrent.
#
# Prerequisites:
#   1. cactus_all (job 56180) finished + pang_1001gplus_all.vcf.gz exists
#   2. build_pangenie_index.sh has run (pang_69 graph indexed)
#   3. loo_preprocessed/ populated (preprocess_loo_one.sh array done)
#   4. cactus truth VCF accessible at $PANG69_VCF for concordance step
# =============================================================================
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
PANG69_VCF=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz

# Sanity checks
if [ ! -s $PANG69_VCF ]; then
    echo "ERROR: pang_69 VCF not found ($PANG69_VCF). Wait for cactus_all (56180)." >&2
    exit 1
fi
if [ ! -d $BASE/data/loo_preprocessed ] || [ "$(ls -A $BASE/data/loo_preprocessed 2>/dev/null | wc -l)" -lt 50 ]; then
    echo "ERROR: loo_preprocessed/ is empty or sparse. Run launch_loo_preprocess.sh first." >&2
    exit 1
fi

N=$(awk 'NR>1' $BASE/data/loo_ena_manifest.tsv | wc -l)
ARRAY="${1:-1-${N}%8}"

JOB=$(sbatch --parsable --array=$ARRAY $BASE/scripts/pangenie_loo_one.sh)
echo "Submitted LOO PanGenie array $JOB ($ARRAY)"
echo ""
echo "Per-sample concordance: $BASE/data/loo_concordance/<ecotype>_summary.tsv"

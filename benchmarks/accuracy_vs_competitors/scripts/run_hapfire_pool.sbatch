#!/bin/bash
#SBATCH --job-name=hapfire_pool
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/hfpool_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/hfpool_%x_%j.out
set -eo pipefail
# hapFIRE (SNP-only) for ONE p80 pool, scored vs recomb_truth -> competitor rows file.
# align reads -> TAIR10 BAM -> hapFIRE on the shared SNP panel -> score SNP AF.
ROOT=/global/scratch/users/tbellg/kmate; cd "$ROOT"
POOL=${1:?pool}; PANEL=${2:-p80}; T=${SLURM_CPUS_PER_TASK:-8}
SC=benchmarks/accuracy_vs_competitors/scripts
WORK=benchmarks/accuracy_vs_competitors/work
RES=benchmarks/accuracy_vs_competitors/results
ROWS=benchmarks/benchmark_runs/rows_competitors/${POOL}_hapfire.tsv
mkdir -p benchmarks/benchmark_runs/rows_competitors
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
SIMD=benchmarks/$PANEL/sims/$POOL
TRUTH=$SIMD/recomb_truth.tsv.gz; [ -f "$TRUTH" ] || TRUTH=$SIMD/recomb_truth_raw.tsv.gz
BAM=$WORK/${POOL}.tair10.srt.bam
source ~/miniforge3/etc/profile.d/conda.sh

# 1. align (reuse if present)
if [ ! -s "$BAM" ]; then
  conda activate kmate
  echo "[$(date)] align -> $BAM"; bash $SC/align_pool.sh "$POOL" "$PANEL" "$T"
fi

# 2. hapFIRE (shared SNP panel)
conda activate hapfire
export PATH="$ROOT/external/HapFIRE/src:$PATH"
echo "[$(date)] hapFIRE on $(hostname): $POOL"
( cd $RES && python $ROOT/external/HapFIRE/hapFIRE.py \
    -vcf $ROOT/$WORK/shared_snps_${PANEL}_Chr1.vcf.gz \
    -bam $ROOT/$BAM -fa "$REF" -output "hapfire_${POOL}" )

# 3. score vs recomb_truth -> competitor rows
conda activate kmate
echo "[$(date)] score hapFIRE vs recomb_truth"
python3 benchmarks/score_competitor_snp.py --tool hapFIRE --pool "$POOL" \
  --truth "$TRUTH" --est $RES/hapfire_${POOL}_snp_frequency.txt --table "$ROWS"
echo "[$(date)] DONE $POOL"

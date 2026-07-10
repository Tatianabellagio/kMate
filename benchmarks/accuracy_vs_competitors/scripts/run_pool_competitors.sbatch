#!/bin/bash
#SBATCH --job-name=acc_compet
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=05:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/compet_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/compet_%x_%j.err
set -eo pipefail
# End-to-end competitor comparison for ONE g0 p80 pool:
#   align reads -> hapFIRE (shared SNP panel) -> vg giraffe (graph) -> truth -> score.
# All tools estimate AF for the SAME shared SNP panel; scored vs the same MAR truth.
# Usage: sbatch run_pool_competitors.sbatch <POOL> [PANEL]
ROOT=/global/scratch/users/tbellg/kmate
cd "$ROOT"
POOL=${1:?need POOL e.g. cov10_n50_g0_s42_hotspots_p80_chr1}
PANEL=${2:-p80}
SC=benchmarks/accuracy_vs_competitors/scripts
WORK=benchmarks/accuracy_vs_competitors/work
RES=benchmarks/accuracy_vs_competitors/results
T=${SLURM_CPUS_PER_TASK:-8}
source ~/miniforge3/etc/profile.d/conda.sh

# regime = strip covNN_ prefix and _s42_..._chr1 suffix  (e.g. n50_g0); keep coverage
COV=$(echo "$POOL" | sed -E 's/^cov([0-9]+)_.*/\1/')
REG=$(echo "$POOL" | sed -E 's/^cov[0-9]+_//; s/_s42_.*$//')
KM=$(ls benchmarks/$PANEL/results/cactus_em_global_filt2_mb/$REG/*_cov${COV}_*.tsv 2>/dev/null | head -1)
echo "[$(date)] POOL=$POOL PANEL=$PANEL REG=$REG host=$(hostname) cpus=$T"
echo "  kMate est = ${KM:?no kMate est found for $REG}"

echo "[$(date)] === align reads -> TAIR10 BAM ==="
conda activate kmate
bash $SC/align_pool.sh "$POOL" "$PANEL" "$T"

echo "[$(date)] === build truth (MAR, g0) ==="
python $SC/build_truth.py \
  --panel-vcf benchmarks/$PANEL/data/pangenome_${PANEL}_chr1.vcf.gz \
  --sites-vcf $WORK/shared_snps_${PANEL}_Chr1.vcf.gz \
  --weights benchmarks/$PANEL/sims/$POOL/pool_weights.tsv \
  --out $WORK/truth_${POOL}.tsv

echo "[$(date)] === hapFIRE ==="
conda activate hapfire
export PATH="$ROOT/external/HapFIRE/src:$PATH"
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
( cd $RES && /usr/bin/time -v -o "hapfire_${POOL}_time.txt" \
    python $ROOT/external/HapFIRE/hapFIRE.py \
    -vcf $ROOT/$WORK/shared_snps_${PANEL}_Chr1.vcf.gz \
    -bam $ROOT/$WORK/${POOL}.tair10.srt.bam \
    -fa "$REF" -output "hapfire_${POOL}" )

echo "[$(date)] === vg giraffe ==="
conda activate vg
bash $SC/run_vg.sh "$POOL" "$PANEL" "$T"

echo "[$(date)] === score (no-filter + miss90) ==="
conda activate basic
F=$(bcftools query -l benchmarks/$PANEL/data/pangenome_${PANEL}_chr1.vcf.gz 2>/dev/null | wc -l)
for tag in nofilt miss90; do
  THR=""; [ "$tag" = miss90 ] && THR="--miss-thr 0.9"
  python $SC/score_competitors.py \
    --truth $WORK/truth_${POOL}.tsv \
    --kmate "$KM" \
    --hapfire $RES/hapfire_${POOL}_snp_frequency.txt \
    --vg $RES/vg_${POOL}_snp_frequency.txt \
    --out-prefix $RES/${POOL}_${tag} --F "$F" $THR
done
echo "[$(date)] DONE $POOL"

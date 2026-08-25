#!/bin/bash
#SBATCH --job-name=hf_ecotype
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=benchmarks/vs_hapfire_ecotype_resolution/logs/hfeco_%x_%j.out
#SBATCH --error=benchmarks/vs_hapfire_ecotype_resolution/logs/hfeco_%x_%j.out
set -eo pipefail
ROOT=/global/scratch/users/tbellg/kmate; cd "$ROOT"
POOL=${1:?pool}; PANEL=${2:-p80}; T=${SLURM_CPUS_PER_TASK:-8}
SC=benchmarks/vs_vg_giraffe/scripts; WORK=benchmarks/vs_vg_giraffe/work
RES=benchmarks/vs_vg_giraffe/results
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
BAM=$WORK/${POOL}.tair10.srt.bam
source ~/miniforge3/etc/profile.d/conda.sh
if [ ! -s "$BAM" ]; then conda activate kmate; echo "[$(date)] align"; bash $SC/align_pool.sh "$POOL" "$PANEL" "$T"; fi
conda activate hapfire; export PATH="$ROOT/external/HapFIRE/src:$PATH"
echo "[$(date)] hapFIRE: $POOL"
( cd $RES && python $ROOT/external/HapFIRE/hapFIRE.py \
    -vcf $ROOT/$WORK/shared_snps_${PANEL}_Chr1.vcf.gz -bam $ROOT/$BAM -fa "$REF" -output "hapfire_${POOL}" )
conda activate kmate
python3 benchmarks/vs_hapfire_ecotype_resolution/scripts/score_ecotype.py --tool hapFIRE \
  --hapfire-ecotype $RES/hapfire_${POOL}_ecotype_frequency.txt \
  --weights benchmarks/$PANEL/sims/$POOL/pool_weights.tsv --pool "$POOL" --panel "$PANEL" \
  --table benchmarks/vs_hapfire_ecotype_resolution/results/ecotype_rows/${POOL}_hapfire.tsv
echo "[$(date)] DONE $POOL"

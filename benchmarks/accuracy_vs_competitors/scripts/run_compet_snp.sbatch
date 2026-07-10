#!/bin/bash
#SBATCH --job-name=compet_snp
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=10:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/csnp_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/csnp_%x_%j.out
set -eo pipefail
# SNP-side competitors for ONE p80 pool: align reads -> hapFIRE (shared SNP panel)
# + vg giraffe (d2 graph) -> per-SNP AF. Idempotent: skips a step whose output exists.
# Usage: sbatch run_compet_snp.sbatch <POOL> [PANEL=p80]
cd /global/scratch/users/tbellg/kmate
ROOT=$(pwd)
POOL=${1:?need POOL}; PANEL=${2:-p80}; T=${SLURM_CPUS_PER_TASK:-8}
SC=benchmarks/accuracy_vs_competitors/scripts
WORK=benchmarks/accuracy_vs_competitors/work
RES=benchmarks/accuracy_vs_competitors/results
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
source ~/miniforge3/etc/profile.d/conda.sh
echo "[$(date)] POOL=$POOL PANEL=$PANEL host=$(hostname) cpus=$T"

BAM=$WORK/${POOL}.tair10.srt.bam
HF=$RES/hapfire_${POOL}_snp_frequency.txt
VG=$RES/vg_${POOL}_snp_frequency.txt

# 1) align -> TAIR10 BAM (hapFIRE input)
if [ ! -s "$BAM" ]; then
  conda activate kmate
  bash $SC/align_pool.sh "$POOL" "$PANEL" "$T"
else echo "[$(date)] skip align (have $BAM)"; fi

# 2) hapFIRE on the shared SNP panel
if [ ! -s "$HF" ]; then
  conda activate hapfire
  export PATH="$PWD/external/HapFIRE/src:$PATH"
  ( cd $RES && /usr/bin/time -v -o "hapfire_${POOL}_time.txt" \
      python $ROOT/external/HapFIRE/hapFIRE.py \
        -vcf $ROOT/$WORK/shared_snps_${PANEL}_Chr1.vcf.gz \
        -bam $ROOT/$BAM -fa "$REF" -output "hapfire_${POOL}" )
else echo "[$(date)] skip hapFIRE (have $HF)"; fi

# 3) vg giraffe (d2) -> AD-coverage per-SNP AF  (rescue-off: AF-equivalent, ~35x faster)
if [ ! -s "$VG" ]; then
  conda activate vg
  VG_RESCUE=0 bash $SC/run_vg.sh "$POOL" "$PANEL" "$T"   # writes $VG directly
else echo "[$(date)] skip vg-SNP (have $VG)"; fi

echo "[$(date)] DONE SNP-side $POOL"

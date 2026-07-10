#!/bin/bash
#SBATCH --job-name=acc_hapfire
#SBATCH --account=fc_moilab
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/hapfire_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/hapfire_%j.err
set -eo pipefail
# hapFIRE on the SHARED p80 SNP panel + a sim pool's TAIR10-aligned BAM.
# Phases 1-2 only (-haplotype omitted -> BigLD/Phase-4 skipped, no R needed):
# independent-LD partition + HARP like/freq -> per-SNP frequency.
ROOT=/global/scratch/users/tbellg/kmate
POOL=${1:-cov10_n231_g0_s42_hotspots_p80_chr1}
WORK=$ROOT/benchmarks/accuracy_vs_competitors/work
RES=$ROOT/benchmarks/accuracy_vs_competitors/results
HF=$ROOT/external/HapFIRE
PY=/global/home/users/tbellg/miniforge3/envs/hapfire/bin/python
export PATH="$HF/src:$PATH"

VCF=$WORK/shared_snps_p80_Chr1.vcf.gz
BAM=$WORK/${POOL}.tair10.srt.bam
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
cd "$RES"

echo "[$(date)] hapFIRE accuracy run on $(hostname): pool=$POOL"
echo "  VCF=$VCF  BAM=$BAM  harp=$(command -v harp)"
/usr/bin/time -v -o "$RES/hapfire_${POOL}_time.txt" $PY "$HF/hapFIRE.py" \
    -vcf "$VCF" -bam "$BAM" -fa "$REF" \
    -output "hapfire_${POOL}"
echo "[$(date)] done"; ls -la "$RES"/hapfire_${POOL}*

#!/bin/bash
#SBATCH --job-name=bench_hapfire
#SBATCH --account=fc_moilab
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=benchmarks/speed_vs_hapfire/logs/hapfire_%j.out
#SBATCH --error=benchmarks/speed_vs_hapfire/logs/hapfire_%j.err
set -eo pipefail

# Speed benchmark: hapFIRE (= HARP wrapper) on the SAME 231-founder Chr1 pool
# that kMate is timed on (visor_freqk del/rep1/cov50/var_del_1kb_n231_f30).
# Phases 1-2 only (-haplotype False): independent-LD partition (python) +
# HARP like/freq -> per-SNP frequency. BigLD (the R step) never runs.
# Reuse of the prebuilt partition (-block) is available for the full pipeline.

ROOT=/global/scratch/users/tbellg/kmate
WORK=$ROOT/benchmarks/speed_vs_hapfire/work
RES=$ROOT/benchmarks/speed_vs_hapfire/results
HF=$ROOT/external/HapFIRE
PY=/global/home/users/tbellg/miniforge3/envs/hapfire/bin/python

export PATH="$HF/src:$PATH"      # harp binary
mkdir -p "$RES"
cd "$RES"

VCF=$WORK/greneNet_chr1.vcf.gz
BAM=$WORK/sim.chr1named.bam
REF=/global/scratch/users/tbellg/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

echo "[$(date)] hapFIRE benchmark start on $(hostname), ${SLURM_CPUS_PER_TASK:-?} cpus"
echo "  VCF=$VCF"; echo "  BAM=$BAM"; echo "  REF=$REF"; echo "  harp=$(command -v harp)"

# NOTE: -haplotype is OMITTED on purpose. argparse type=bool makes "-haplotype False"
# evaluate to True (bool("False")==True), so passing it would wrongly enable Phase 4
# (BigLD, which needs R). Omitting it leaves the real default False -> Phase 4 skipped.
# We therefore time Phase 1 (independent-LD partition) + Phase 2 (HARP like/freq ->
# per-SNP frequency). BigLD never runs, so the slow partition-build step is skipped
# entirely (matching the "reuse / save that step" intent).
/usr/bin/time -v -o "$RES/hapfire_time.txt" $PY "$HF/hapFIRE.py" \
    -vcf "$VCF" \
    -bam "$BAM" \
    -fa  "$REF" \
    -output hapfire_n231_chr1_cov50_f30

echo "[$(date)] hapFIRE benchmark done"
ls -la "$RES"

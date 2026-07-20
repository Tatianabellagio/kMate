#!/bin/bash
#SBATCH --job-name=p231_gw_hf
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/gw_hf_%j.out
#SBATCH --error=logs/gw_hf_%j.out

# =============================================================================
# ONE-OFF diagnostic: genome-wide hapFIRE (all 5 chroms) vs the Chr1-only run,
# same reads-generation recipe, testing whether the h/ecotype R2 collapse at
# high N (-0.02 to -0.20 at N=150 Chr1-only) is explained by hapFIRE's
# block-averaging mechanism needing genome-wide scope (per
# docs/FOUNDER_NORMALIZATION_FIX.md sec 4.4 + the real production recipe in
# panel/arch3/chr1/jobA7_compare_vs_hapfire.sh's commands.sh) rather than a
# pipeline bug. Reads from run_sim_p231_genomewide_test.sh.
#
# Usage: sbatch run_hapfire_genomewide_test.sh N SEED
# =============================================================================
set -euo pipefail
N=${1:?Usage: N SEED}
SEED=${2:?Usage: N SEED}
COV=10

ROOT=/global/scratch/users/tbellg/kmate
WORK=$ROOT/benchmarks/speed_vs_hapfire/sims_genomewide_test/cov${COV}_n${N}_g0_s${SEED}_p231_genomewide
OUT=$ROOT/benchmarks/speed_vs_hapfire/results/genomewide_test/n${N}_s${SEED}
mkdir -p "$OUT" logs

BWA=/global/home/users/tbellg/miniforge3/envs/kmate/bin/bwa
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/kmate/bin/samtools
HFPY=/global/home/users/tbellg/miniforge3/envs/hapfire/bin/python
HF=$ROOT/external/HapFIRE
REF=/global/scratch/users/tbellg/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
VCF=$ROOT/benchmarks/speed_vs_hapfire/work/greneNet_final_v1.1.recode.vcf.gz

READS=$WORK/reads
[ -s "$READS/r1.fq" ] || { echo "ERROR: missing $READS/r1.fq" >&2; exit 1; }
[ -s "${REF}.bwt" ]    || { echo "ERROR: missing bwa index $REF.bwt" >&2; exit 1; }

SAMPLE=p231_hapfire_gw_n${N}_s${SEED}
cd "$OUT"

echo "[$(date)] ALIGN (genome-wide) N=$N s=$SEED  (host=$(hostname), cpus=${SLURM_CPUS_PER_TASK:-8})"
/usr/bin/time -v -o "$OUT/align_time.txt" bash -c "
    set -euo pipefail
    $BWA mem -t ${SLURM_CPUS_PER_TASK:-8} $REF $READS/r1.fq $READS/r2.fq 2> $OUT/bwa.log \
        | $SAMTOOLS sort -@ ${SLURM_CPUS_PER_TASK:-8} -o $OUT/sim.srt.bam -
    $SAMTOOLS index $OUT/sim.srt.bam
"

# reheader ALL 5 chroms Chr1-5 -> 1-5 to match the genome-wide greneNet VCF
$SAMTOOLS view -H "$OUT/sim.srt.bam" \
    | sed -E 's/SN:Chr([1-5])\b/SN:\1/' > "$OUT/newhdr.sam"
$SAMTOOLS reheader "$OUT/newhdr.sam" "$OUT/sim.srt.bam" > "$OUT/sim.renamed.bam"
$SAMTOOLS index "$OUT/sim.renamed.bam"
rm -f "$OUT/sim.srt.bam" "$OUT/sim.srt.bam.bai" "$OUT/newhdr.sam"

echo "[$(date)] HAPFIRE (genome-wide) N=$N s=$SEED"
export PATH="$HF/src:$PATH"
/usr/bin/time -v -o "$OUT/hapfire_time.txt" $HFPY "$HF/hapFIRE.py" \
    -vcf "$VCF" \
    -bam "$OUT/sim.renamed.bam" \
    -fa  "$REF" \
    -output $SAMPLE

echo "[$(date)] DONE N=$N s=$SEED"
ls -la "$OUT"

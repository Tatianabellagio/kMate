#!/bin/bash
#SBATCH --job-name=oldpanel_hf
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=logs/oldpanel_hf_%j.out
#SBATCH --error=logs/oldpanel_hf_%j.out

# ONE-OFF diagnostic: run hapFIRE on OUR simulated reads (39-founder-overlap
# pool) but using the OLD/prior "200_test_chr1.recode.vcf" panel instead of
# our 231-founder greneNet_chr1.vcf.gz. Isolates: does OUR sim/alignment
# pipeline work fine when paired with the panel type that gave good published
# -matching results before, or does it still collapse?
set -euo pipefail
N=${1:?Usage: N SEED}
SEED=${2:?Usage: N SEED}
COV=10

ROOT=/global/scratch/users/tbellg/kmate
WORK=$ROOT/benchmarks/speed_vs_hapfire/sims_oldpanel_test/cov${COV}_n${N}_g0_s${SEED}_p231_oldpanel39
OUT=$ROOT/benchmarks/speed_vs_hapfire/results/oldpanel_test/n${N}_s${SEED}
mkdir -p "$OUT" logs

BWA=/global/home/users/tbellg/miniforge3/envs/kmate/bin/bwa
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/kmate/bin/samtools
HFPY=/global/home/users/tbellg/miniforge3/envs/hapfire/bin/python
HF=$ROOT/external/HapFIRE
REF=/global/scratch/users/tbellg/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
VCF=/global/scratch/projects/fc_moilab/projects/grenenet-phase1/drive_zenodo/data-intermediate/hapfire_performance/200_test_chr1.recode.vcf

READS=$WORK/reads
[ -s "$READS/r1.fq" ] || { echo "ERROR: missing $READS/r1.fq" >&2; exit 1; }

SAMPLE=p231_hapfire_oldpanel_n${N}_s${SEED}
cd "$OUT"

echo "[$(date)] ALIGN (old-panel test) N=$N s=$SEED"
/usr/bin/time -v -o "$OUT/align_time.txt" bash -c "
    set -euo pipefail
    $BWA mem -t ${SLURM_CPUS_PER_TASK:-8} $REF $READS/r1.fq $READS/r2.fq 2> $OUT/bwa.log \
        | $SAMTOOLS sort -@ ${SLURM_CPUS_PER_TASK:-8} -o $OUT/sim.srt.bam -
    $SAMTOOLS index $OUT/sim.srt.bam
"

echo "[$(date)] HAPFIRE (old-panel test, 200_test_chr1.recode.vcf) N=$N s=$SEED"
export PATH="$HF/src:$PATH"
/usr/bin/time -v -o "$OUT/hapfire_time.txt" $HFPY "$HF/hapFIRE.py" \
    -vcf "$VCF" \
    -bam "$OUT/sim.srt.bam" \
    -fa  "$REF" \
    -output $SAMPLE

echo "[$(date)] DONE N=$N s=$SEED"
ls -la "$OUT"

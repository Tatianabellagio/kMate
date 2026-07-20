#!/bin/bash
#SBATCH --job-name=gren_hf
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=logs/gren_hf_%j.out
#SBATCH --error=logs/gren_hf_%j.out

# FAIR hapFIRE run: reads simulated from greneNet-derived fastas, aligned to
# TAIR10 (ref_xing, contig "1"), hapFIRE tested against the SAME greneNet VCF
# (contig "1"). Self-consistent panel -- reads and reference come from the same
# SNP catalog, so no arch3 indel/SV reads to trip HARP. Mirror of kMate-on-arch3.
#
# Usage: sbatch run_hapfire_greneNet.sh N COV SEED
set -euo pipefail
mkdir -p logs
N=${1:?Usage: N COV SEED}
COV=${2:?Usage: N COV SEED}
SEED=${3:?Usage: N COV SEED}

ROOT=/global/scratch/users/tbellg/kmate
BASE=$ROOT/benchmarks/speed_vs_hapfire
WORK=$BASE/sims_greneNet/cov${COV}_n${N}_g0_s${SEED}_greneNet_chr1
OUT=$BASE/results/greneNet_fair/n${N}_cov${COV}_s${SEED}
mkdir -p "$OUT" logs

BWA=/global/home/users/tbellg/miniforge3/envs/kmate/bin/bwa
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/kmate/bin/samtools
HFPY=/global/home/users/tbellg/miniforge3/envs/hapfire/bin/python
HF=$ROOT/external/HapFIRE
REF=/global/scratch/users/tbellg/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
VCF=$BASE/work/greneNet_chr1.vcf.gz          # contig "1", matches BAM

READS=$WORK/reads
[ -s "$READS/r1.fq" ] || { echo "ERROR: missing $READS/r1.fq" >&2; exit 1; }
[ -s "${REF}.bwt" ]    || { echo "ERROR: missing bwa index $REF.bwt" >&2; exit 1; }

SAMPLE=greneNet_hapfire_n${N}_cov${COV}_s${SEED}
cd "$OUT"

echo "[$(date)] ALIGN (greneNet fair) N=$N cov=$COV s=$SEED"
/usr/bin/time -v -o "$OUT/align_time.txt" bash -c "
    set -euo pipefail
    $BWA mem -t ${SLURM_CPUS_PER_TASK:-8} $REF $READS/r1.fq $READS/r2.fq 2> $OUT/bwa.log \
        | $SAMTOOLS sort -@ ${SLURM_CPUS_PER_TASK:-8} -o $OUT/sim.srt.bam -
    $SAMTOOLS index $OUT/sim.srt.bam
"

echo "[$(date)] HAPFIRE (greneNet fair) N=$N cov=$COV s=$SEED"
export PATH="$HF/src:$PATH"
/usr/bin/time -v -o "$OUT/hapfire_time.txt" $HFPY "$HF/hapFIRE.py" \
    -vcf "$VCF" \
    -bam "$OUT/sim.srt.bam" \
    -fa  "$REF" \
    -output $SAMPLE

echo "[$(date)] DONE N=$N cov=$COV s=$SEED"
ls -la "$OUT"
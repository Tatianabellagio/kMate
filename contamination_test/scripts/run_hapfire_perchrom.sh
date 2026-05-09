#!/bin/bash
#SBATCH --job-name=hapfire_1141
#SBATCH --output=/home/tbellagio/scratch/hapfire_sv/contamination_test/logs/hapfire_%x_%j.out
#SBATCH --error=/home/tbellagio/scratch/hapfire_sv/contamination_test/logs/hapfire_%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=600G
# Default --mem=600G is enough for chr1 (peak measured 503 GB at F=1141).
# submit_all.sh overrides per-chrom (chr2/chr4: 450G, chr3: 500G, chr5: 550G,
# chr1: 600G) so SLURM can pack 2-3 concurrent jobs on the bse-2021-001 1.5TB node.
# Only bse-2021-001 has enough memory; jobs at any of these sizes go there.

# Run hapFIRE on one (sample, chrom) pair against the 1141-panel.
#
# Args (positional):
#   $1 = sample tag, e.g. SEEDMIX_S1   -> uses xwu's filtered_seeds-${N}.bam
#   $2 = chrom number 1..5
#
# Output is written to:
#   results/<sample_tag>_chr<N>/<sample_tag>_chr<N>_*.txt

set -eo pipefail

SAMPLE="${1:?sample tag required, e.g. SEEDMIX_S1}"
CH="${2:?chrom number required, 1..5}"

# Sample tag -> BAM map: SEEDMIX_S{1..8} -> xwu's filtered_seeds-{1..8}.bam
N="${SAMPLE#SEEDMIX_S}"
BAM="/carnegie/nobackup/scratch/xwu/GrENE_net/seed_mix/filtered_bam/filtered_seeds-${N}.bam"

ROOT=/home/tbellagio/scratch/hapfire_sv/contamination_test
# rsynced from xwu's Savio (xingwu@savio:/global/scratch/users/xingwu/GrENE_net/vcf/imputation/).
# 1141 ecotypes, fully phased, MAC>=7, biallelic SNPs, all 6 GrENE extras present.
VCF="${ROOT}/vcf/1001G_80pilot_israel_regmap_overlapping_biallelic_chr${CH}_mac7.recode.vcf"
REF=/home/tbellagio/scratch/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

HAPFIRE_DIR=/carnegie/nobackup/scratch/xwu/haplotype_frequency_estimation/hapFIRE_sourcecode
HAPFIRE=${HAPFIRE_DIR}/hapFIRE.py
HARP_BIN=${HAPFIRE_DIR}/bin
export PATH="${HARP_BIN}:${PATH}"

PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python

OUTDIR="${ROOT}/results/${SAMPLE}_chr${CH}"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# pre-flight checks (cheap)
[[ -f "${VCF}" ]] || { echo "missing VCF: ${VCF}"; exit 2; }
[[ -f "${BAM}" ]] || { echo "missing BAM: ${BAM}"; exit 2; }
[[ -f "${BAM}.bai" ]] || { echo "missing BAM index: ${BAM}.bai"; exit 2; }
[[ -f "${REF}" ]] || { echo "missing REF: ${REF}"; exit 2; }
[[ -f "${REF}.fai" ]] || { echo "missing REF index: ${REF}.fai"; exit 2; }

OUT_PREFIX="${SAMPLE}_chr${CH}"

echo "[$(date)] hapFIRE start: sample=${SAMPLE} chr=${CH}"
echo "  VCF: ${VCF}"
echo "  BAM: ${BAM}"
echo "  REF: ${REF}"
echo "  OUT: ${OUTDIR}/${OUT_PREFIX}"

${PYTHON} ${HAPFIRE} \
    -v "${VCF}" \
    -b "${BAM}" \
    -f "${REF}" \
    -o "${OUT_PREFIX}"

echo "[$(date)] hapFIRE done"
ls -la "${OUTDIR}"

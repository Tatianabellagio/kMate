#!/bin/bash
#SBATCH --job-name=imp_231
#SBATCH --partition=bse
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/imp_merged_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/imp_merged_%j.err

# =============================================================================
# impute_merged_panel.sh
# Beagle 5.x imputation on the merged 231-founder VCF (founders_231_chr.vcf.gz)
# to fill missing GTs across both cohorts (mostly cactus haploid `.`).
#
# Following xwu's pattern (vcf/imputation/commands.sh):
#   - per-chromosome (memory-bounded; Beagle holds full chrom in RAM)
#   - window=10000  overlap=1000  ne=10000  (A. thaliana-tuned)
#   - 16 threads, -Xmx48g
#
# Differences from xwu's default and our 03_run_imputation.sh:
#   - No separate ref panel — use merged VCF as both gt and ref. Beagle 5.x's
#     gt= mode imputes internal missing using LD within the input.
#   - SVs included (not just biallelic SNPs). Multi-allelic bubbles split via
#     bcftools norm -m -any so Beagle gets per-ALT biallelic records.
#   - Pre-filter MAC>=1 (drop monomorphic records — Beagle can't impute
#     records with no carrier; saves time + avoids errors).
#   - Diploidize cactus-side haploid GTs first (Beagle requires diploid input).
#
# Outputs:
#   imputation/work/founders_231_imputed.vcf.gz   (full imputed VCF)
# =============================================================================
set -euo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
WORK=$BASE/imputation/work_merged
mkdir -p $WORK $BASE/imputation/logs
cd $WORK

BEAGLE=$BASE/imputation/beagle.jar
BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix

INPUT=$BASE/pangenie_genotyping/data/merged/founders_231_chr.vcf.gz

# ---- Step A: diploidize + decompose multi-allelic + filter monomorphic ----
PREP=$WORK/founders_231_dipl_split.vcf.gz
if [ ! -s $PREP ]; then
    echo "[$(date)] Step A: diploidize haploid GTs + split multi-allelic + filter MAC=0"
    # bcftools +setGT can pad haploid to diploid via -t a -n c:'major'-style operations,
    # but the simplest reliable way for inbred A. thaliana is awk in the body:
    #   X     -> X|X     (haploid call to diploid homozygous)
    #   .     -> .|.     (haploid missing to diploid missing)
    # Then bcftools norm -m -any to decompose multi-allelic records to biallelic.
    # Then bcftools view -e 'AC=0 || AC=AN' to drop monomorphic-after-decomposition.
    $BCF view $INPUT 2>/dev/null | \
        awk 'BEGIN{OFS="\t"} /^#/{print; next}
             {for(i=10;i<=NF;i++) {
                  n=split($i, fmt, ":")
                  gt=fmt[1]
                  if (gt !~ /[\/\|]/) {
                      # haploid: X -> X|X, . -> .|.
                      fmt[1] = gt "|" gt
                      $i = fmt[1]
                      for(j=2;j<=n;j++) $i = $i ":" fmt[j]
                  }
              }
              print
             }' | \
        $BCF norm -m -any -Ou 2>/dev/null | \
        $BCF view -e 'AC=0 || AC=AN' -Oz -o $PREP --threads 8
    $TABIX -p vcf $PREP
fi
echo "[$(date)] PREP records: $($BCF index -n $PREP)"
echo "[$(date)] PREP samples: $($BCF query -l $PREP | wc -l)"

# ---- Step B: per-chromosome Beagle ----
for i in 1 2 3 4 5; do
    OUT_PREFIX=$WORK/imputed_chr${i}
    if [ -s ${OUT_PREFIX}.vcf.gz ]; then
        echo "[$(date)] chr${i}: ${OUT_PREFIX}.vcf.gz exists, skipping"
        continue
    fi
    echo "[$(date)] Beagle on Chr${i}"
    # Beagle 5.x: window/overlap are in cM (NOT bp like 4.x). Defaults are
    # window=40, overlap=2 (cM) — fine for our 30 Mb A. thaliana chroms.
    # ne=10000 is reasonable for an inbred A. thaliana panel.
    java -Xmx48g -jar $BEAGLE \
        gt=$PREP \
        chrom=Chr${i} \
        out=$OUT_PREFIX \
        nthreads=16 \
        ne=10000 \
        seed=42
    $TABIX -p vcf ${OUT_PREFIX}.vcf.gz
    echo "[$(date)] Chr${i}: $($BCF index -n ${OUT_PREFIX}.vcf.gz) records"
done

# ---- Step C: concat per-chrom outputs ----
FINAL=$WORK/founders_231_imputed.vcf.gz
echo "[$(date)] Concat all chroms -> $FINAL"
$BCF concat \
    $WORK/imputed_chr1.vcf.gz $WORK/imputed_chr2.vcf.gz \
    $WORK/imputed_chr3.vcf.gz $WORK/imputed_chr4.vcf.gz \
    $WORK/imputed_chr5.vcf.gz \
    -Oz -o $FINAL --threads 8
$TABIX -p vcf $FINAL

echo "[$(date)] DONE"
echo "  records: $($BCF index -n $FINAL)"
echo "  samples: $($BCF query -l $FINAL | wc -l)"
ls -lh $FINAL

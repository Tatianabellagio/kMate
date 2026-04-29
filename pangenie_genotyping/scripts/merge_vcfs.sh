#!/bin/bash
#SBATCH --job-name=merge_pg
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/merge_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/merge_%j.err

# =============================================================================
# merge_vcfs.sh
# Stage 4: build the combined 231-founder VCF.
#   80 cactus accessions   → from cactus pang_69 VCF (long-read truth)
#   151 short-read founders → from per-sample PanGenie VCFs
# bcftools merge → 231-founder catalog → ready for cn_kmer / cn_var build.
# =============================================================================
set -euo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping
GT_DIR=$BASE/data/genotyped
OUT_DIR=$BASE/data/merged
mkdir -p $OUT_DIR

BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix

# Step A: subset pang_69 VCF to the 80 cactus accessions that are in our 231 panel
PANG69_VCF=/home/tbellagio/scratch/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz
PANEL_MAP=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/sv_panel_to_accession_id.tsv
KEEP_80=$OUT_DIR/cactus_overlap_80.txt
awk 'NR>1 {print $3}' $PANEL_MAP | sort -u > $KEEP_80   # 1001G IDs of cactus founders
echo "[$(date)] cactus founders to keep: $(wc -l < $KEEP_80)"

# Need to rename cactus pang assembly IDs → 1001G IDs first (same as imputation pipeline did)
# Use the existing rename file from the imputation step
RENAMED_CACTUS=$OUT_DIR/cactus_pang69_1001g.vcf.gz
if [ ! -s $RENAMED_CACTUS ]; then
    echo "[$(date)] rename cactus assembly IDs → 1001G IDs"
    $BCF reheader -s /carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work/sample_rename.txt \
        $PANG69_VCF | $BCF view -S $KEEP_80 --force-samples -Oz -o $RENAMED_CACTUS --threads 8
    $TABIX -p vcf $RENAMED_CACTUS
fi

# Step B: collect all PanGenie per-sample VCFs (151 of them) into a single multi-sample VCF
# bcftools merge can do this directly across files
echo "[$(date)] collect 151 PanGenie VCFs"
PG_VCFS=$(ls $GT_DIR/*.vcf.gz 2>/dev/null | sort)
N_PG=$(echo "$PG_VCFS" | wc -l)
echo "  found $N_PG PanGenie VCFs (expected 151)"

PG_MERGED=$OUT_DIR/pangenie_151.vcf.gz
if [ ! -s $PG_MERGED ]; then
    $BCF merge $PG_VCFS -Oz -o $PG_MERGED --threads 8
    $TABIX -p vcf $PG_MERGED
fi

# Step C: merge cactus_80 + pangenie_151 → 231-founder VCF
COMBINED=$OUT_DIR/founders_231.vcf.gz
echo "[$(date)] merge cactus_80 + pangenie_151 → $COMBINED"
$BCF merge $RENAMED_CACTUS $PG_MERGED -Oz -o $COMBINED --threads 8
$TABIX -p vcf $COMBINED

# Step D: rename chroms 1..5 → Chr1..Chr5 (cactus convention used by cn_build)
COMBINED_CHR=$OUT_DIR/founders_231_chr.vcf.gz
echo -e "1\tChr1\n2\tChr2\n3\tChr3\n4\tChr4\n5\tChr5" > $OUT_DIR/rename_chrs.txt
$BCF annotate --rename-chrs $OUT_DIR/rename_chrs.txt $COMBINED \
    -Oz -o $COMBINED_CHR --threads 8
$TABIX -p vcf $COMBINED_CHR

echo "[$(date)] DONE"
echo "  records: $($BCF view -H $COMBINED_CHR | wc -l)"
echo "  samples: $($BCF query -l $COMBINED_CHR | wc -l)"
ls -lh $COMBINED_CHR

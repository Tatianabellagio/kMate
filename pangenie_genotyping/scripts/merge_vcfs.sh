#!/bin/bash
#SBATCH --job-name=merge_pg
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/merge_%j.out
#SBATCH --error=logs/merge_%j.err

# =============================================================================
# merge_vcfs.sh
# Stage 4: build the combined 231-founder VCF.
#   80 cactus accessions   → from cactus pang_69 VCF (long-read truth)
#   151 short-read founders → from per-sample PanGenie VCFs
# bcftools merge → 231-founder catalog → ready for cn_kmer / cn_var build.
# =============================================================================
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv/pangenie_genotyping
GT_DIR=$BASE/data/genotyped
OUT_DIR=$BASE/data/merged
mkdir -p $OUT_DIR

BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix

# Optional GQ filter on PanGenie cells. PanGenie reports GQ per genotype call;
# cells below this are masked to ./. so they don't pollute the merged catalog.
# Set to 0 (or unset) to disable.
PG_GQ_MIN=${PG_GQ_MIN:-0}

# ---- Step A: subset cactus pang_all VCF to the 80 GrENE-overlap founders ----
# pang_all has 135 cactus founders. Of those, 80 map to GrENE-Net 231 ecotype IDs
# (sample_rename.txt). The other 55 are dropped — they're not in GrENE-Net.
# Chroms in pang are already Chr1..Chr5 (verified) — no chrom rename needed.
PANG69_VCF=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz
SAMPLE_RENAME=/global/scratch/users/tbellg/hapfire_sv/imputation/work/sample_rename.txt
KEEP_80=$OUT_DIR/cactus_overlap_80.txt
awk '{print $2}' $SAMPLE_RENAME | sort -u > $KEEP_80
echo "[$(date)] cactus founders to keep (1001G IDs): $(wc -l < $KEEP_80)"

RENAMED_CACTUS=$OUT_DIR/cactus_pang69_1001g.vcf.gz
if [ ! -s $RENAMED_CACTUS ]; then
    echo "[$(date)] rename cactus assembly IDs → 1001G IDs and subset to GrENE-overlap 80"
    # No --force-samples: every ID in KEEP_80 must exist after rename, else fail loud.
    $BCF reheader -s $SAMPLE_RENAME $PANG69_VCF | \
        $BCF view -S $KEEP_80 -Oz -o $RENAMED_CACTUS --threads 8
    $TABIX -p vcf $RENAMED_CACTUS
fi
N_CACTUS=$($BCF query -l $RENAMED_CACTUS | wc -l)
echo "  cactus subset samples: $N_CACTUS (expected 80)"
[ "$N_CACTUS" = "80" ] || { echo "ERROR: expected 80 cactus samples, got $N_CACTUS"; exit 1; }

# ---- Step B: collect & merge per-sample PanGenie VCFs ------------------------
echo "[$(date)] collect PanGenie per-sample VCFs"
PG_VCFS=$(ls $GT_DIR/*.vcf.gz 2>/dev/null | sort)
N_PG=$(echo "$PG_VCFS" | wc -l)
echo "  found $N_PG PanGenie VCFs (expected 151)"
[ "$N_PG" = "151" ] || { echo "ERROR: expected 151 PanGenie VCFs, got $N_PG"; exit 1; }

# Optional per-sample GQ filter — mask cells with GQ < $PG_GQ_MIN to ./.
if [ "$PG_GQ_MIN" -gt 0 ]; then
    echo "[$(date)] applying GQ < $PG_GQ_MIN → ./. on each PanGenie VCF"
    FILT_DIR=$OUT_DIR/pangenie_gqfilt
    mkdir -p $FILT_DIR
    NEW_PG_VCFS=""
    for v in $PG_VCFS; do
        out=$FILT_DIR/$(basename $v)
        if [ ! -s $out ]; then
            $BCF +setGT $v -Oz -o $out -- -t q -n . -e "FMT/GQ>=${PG_GQ_MIN}"
            $TABIX -p vcf $out
        fi
        NEW_PG_VCFS="$NEW_PG_VCFS $out"
    done
    PG_VCFS="$NEW_PG_VCFS"
fi

PG_MERGED=$OUT_DIR/pangenie_151.vcf.gz
if [ ! -s $PG_MERGED ]; then
    $BCF merge $PG_VCFS -Oz -o $PG_MERGED --threads 8
    $TABIX -p vcf $PG_MERGED
fi

# ---- Step C: merge cactus_80 + pangenie_151 → 231-founder VCF ----------------
COMBINED=$OUT_DIR/founders_231_chr.vcf.gz   # chroms already Chr1..Chr5
echo "[$(date)] merge cactus_80 + pangenie_151 → $COMBINED"
$BCF merge $RENAMED_CACTUS $PG_MERGED -Oz -o $COMBINED --threads 8
$TABIX -p vcf $COMBINED

N_FINAL=$($BCF query -l $COMBINED | wc -l)
N_REC=$($BCF view -H $COMBINED | wc -l)
echo "[$(date)] DONE  records=$N_REC  samples=$N_FINAL"
[ "$N_FINAL" = "231" ] || { echo "ERROR: expected 231 founders in merged VCF, got $N_FINAL"; exit 1; }
ls -lh $COMBINED

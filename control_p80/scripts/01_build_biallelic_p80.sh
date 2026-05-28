#!/bin/bash
#SBATCH --job-name=p80_a1_bial
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=logs/01_biallelic_%j.out
#SBATCH --error=logs/01_biallelic_%j.err

# =============================================================================
# Phase A (arch decomposition) -- canonical p80 biallelic VCF for Chr1.
#
# Shared prerequisite: panel/arch3/chr1's A1 outputs. A1 took the FULL 135-assembly
# pangenome VCF (pang_1001gplus_all.vcf.gz) + the 135-asm GFA, ran annotate_vcf,
# and produced the symbolic INFO/ID annotation + biallelic catalog. Built once;
# we reuse those files here so this script is just the cactus-side conversion.
#
# Why no transfer_id (cf. panel/arch3/chr1/jobA3): production cactus_78 came from a
# different source (cactus_pang69_1001g) than the 135-asm annotated catalog,
# so transfer_id was needed to propagate INFO/ID. Here we subset the 135-asm
# annotated VCF DIRECTLY -- INFO/ID is already present.
#
# Pipeline:
#   1. Build 80-Asm_ID inclusion list = (82-acc samples) − (data/exclude_list.txt)
#   2. bcftools view -S 80-list panel/arch3/chr1/chr1_135_annotated.sorted.vcf.gz
#      -> 80-sample annotated multi-allelic VCF
#   3. convert-to-biallelic.py against arch3 biallelic catalog
#      -> biallelic VCF (one row per atomic variant ID; carrier GT per sample)
#   4. sort + bgzip + tabix
#   5. fill-tags AC,AN,F_MISSING (convert-to-biallelic strips AC)
#   6. drop AC=0 (records carried only by excluded pair become orphans)
#   7. reheader Assembly_ID -> Accession_ID
#
# Output: data/pangenome_p80_chr1.vcf.gz (biallelic, 80 samples, Acc_ID-named)
# =============================================================================
set -euo pipefail

CTRL=/global/scratch/users/tbellg/kmate/control_p80
BASE=/global/scratch/users/tbellg/kmate

BCF=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/gwas/bin/tabix
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CONVERT=$BASE/external_tools/pangenie-tools/pipelines/run-from-callset/scripts/convert-to-biallelic.py

# Shared prerequisite: arch3 Chr1 A1 outputs (135-asm annotated + biallelic catalog)
A1_ANNOT=$BASE/panel/arch3/chr1/chr1_135_annotated.sorted.vcf.gz
A1_BIAL=$BASE/panel/arch3/chr1/chr1_135_annotated_biallelic.sorted.vcf.gz

# 82-acc source VCF (used only to get the canonical 82 Asm_ID list)
SRC_82=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.vcf.gz

PANEL_TSV=$BASE/data/sv_panel_to_accession_id.tsv
EXCLUDE_LIST=$BASE/data/exclude_list.txt

WORK=$CTRL/data
mkdir -p $WORK $CTRL/logs

for f in "$A1_ANNOT" "$A1_BIAL" "$SRC_82" "$PANEL_TSV" "$EXCLUDE_LIST" "$BCF" "$BGZIP" "$TABIX" "$PY" "$CONVERT"; do
    [ -e "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

OUT=$WORK/pangenome_p80_chr1.vcf.gz
if [ -s "$OUT" ] && [ -s "$OUT.tbi" ]; then
    echo "[$(date)] $OUT already exists -- nothing to do."
    $BCF query -l $OUT | wc -l | awk '{print "  samples: "$0}'
    $BCF view -H $OUT 2>/dev/null | wc -l | awk '{print "  records: "$0}'
    exit 0
fi

EXCLUDE_CSV=$(paste -sd, $EXCLUDE_LIST)
echo "[$(date)] exclude Assembly_IDs: $EXCLUDE_CSV"

# Step 1: build 80-Asm_ID inclusion list (82-acc samples − exclude_list)
KEEP_LIST=$WORK/samples_80_asm.txt
$BCF query -l $SRC_82 | awk -v ex="$EXCLUDE_CSV" 'BEGIN{n=split(ex,a,",");for(i=1;i<=n;i++)E[a[i]]=1} !($1 in E){print}' > $KEEP_LIST
N_KEEP=$(wc -l < $KEEP_LIST)
echo "[$(date)] Step 1: keep list -> $KEEP_LIST ($N_KEEP entries, expected 80)"
[ "$N_KEEP" -eq 80 ] || { echo "ERROR: keep list has $N_KEEP rows, expected 80" >&2; exit 1; }

# Build rename map (Asm_ID -> Acc_ID) for the 80 we keep
RENAME=$WORK/samples_80_rename.tsv
tail -n +2 $PANEL_TSV \
    | awk -F'\t' -v ex="$EXCLUDE_CSV" 'BEGIN{n=split(ex,a,",");for(i=1;i<=n;i++)E[a[i]]=1} !($1 in E){print $1"\t"$3}' \
    > $RENAME
[ "$(wc -l < $RENAME)" -eq 80 ] || { echo "ERROR: rename map has $(wc -l < $RENAME) rows" >&2; exit 1; }

# Step 2: subset 135-asm annotated VCF to the 80 Asm_IDs (preserves INFO/ID)
echo ""
echo "[$(date)] === Step 2: subset 135-asm annotated VCF to 80 Asm_IDs ==="
STEP2=$WORK/_p80_step2_subset.vcf
$BCF view -S $KEEP_LIST $A1_ANNOT > $STEP2
N_SAMP=$($BCF query -l $STEP2 | wc -l)
N_REC=$(grep -vc '^#' $STEP2)
echo "  samples: $N_SAMP (expected 80)"
echo "  records: $N_REC"
[ "$N_SAMP" -eq 80 ] || { echo "ERROR: got $N_SAMP samples after subset" >&2; exit 1; }

# Step 3: convert-to-biallelic (symbolic decomposition using A1 biallelic catalog)
echo ""
echo "[$(date)] === Step 3: convert-to-biallelic ==="
STEP3=$WORK/_p80_step3_biallelic.vcf
cat $STEP2 | $PY $CONVERT $A1_BIAL > $STEP3 2> $WORK/_p80_convert.log
echo "  records emitted: $(grep -vc '^#' $STEP3)"
rm -f $STEP2

# Step 4: sort + bgzip + tabix
echo ""
echo "[$(date)] === Step 4: sort + bgzip + tabix ==="
STEP4_BG=$WORK/_p80_step4_sorted.vcf.gz
awk '$1 ~ /^#/ {print $0; next} {print $0 | "sort -k1,1 -k2,2n"}' $STEP3 \
    | $BGZIP -c > $STEP4_BG
$TABIX -p vcf $STEP4_BG
rm -f $STEP3

# Step 5: fill-tags AC,AN,F_MISSING
echo ""
echo "[$(date)] === Step 5: fill-tags AC,AN,F_MISSING ==="
STEP5=$WORK/_p80_step5_filled.vcf.gz
$BCF +fill-tags $STEP4_BG --threads 4 -Oz -o $STEP5 -- -t AC,AN,F_MISSING
$TABIX -p vcf $STEP5
rm -f $STEP4_BG ${STEP4_BG}.tbi

# Step 6: drop AC=0
echo ""
echo "[$(date)] === Step 6: drop AC=0 ==="
STEP6=$WORK/_p80_step6_acpos.vcf.gz
N_PRE=$($BCF view -H $STEP5 | wc -l)
$BCF view -e 'INFO/AC=0' $STEP5 --threads 4 -Oz -o $STEP6
$TABIX -p vcf $STEP6
N_POST=$($BCF view -H $STEP6 | wc -l)
echo "  records pre-AC=0:  $N_PRE"
echo "  records post-AC=0: $N_POST  (dropped $((N_PRE - N_POST)))"
rm -f $STEP5 ${STEP5}.tbi

# Step 7: reheader Asm_ID -> Acc_ID
echo ""
echo "[$(date)] === Step 7: reheader Assembly_ID -> Accession_ID ==="
$BCF reheader -s $RENAME $STEP6 > $OUT
$TABIX -p vcf $OUT
rm -f $STEP6 ${STEP6}.tbi

# Sanity
echo ""
echo "[$(date)] DONE -> $OUT"
echo "  samples (Accession_IDs, first 5):"
$BCF query -l $OUT | head -5 | awk '{print "    "$0}'
echo -n "  n samples: "; $BCF query -l $OUT | wc -l
echo -n "  n records: "; $BCF view -H $OUT 2>/dev/null | wc -l
ls -lh $OUT $OUT.tbi | awk '{print "  "$5" "$NF}'

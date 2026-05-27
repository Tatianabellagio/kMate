#!/bin/bash
#SBATCH --job-name=p82_nobub
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/logs/p82_nobub_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/logs/p82_nobub_%j.err

# Build p82 panel using the SAME processing as v3 (no vcfbub, just bcftools norm).
# Mirrors v3's cactus-side pipeline: subset → no vcfbub → norm -m -any → cn_var.
# Output: scratch/p82_nobub/cn_var_p82_nobub.{cn_var,meta}.npz + ref_alt.tsv

set -euo pipefail
ROOT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
SRC=/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.vcf.gz
REF=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.numeric.iupacN.fa
PANEL_TSV=$ROOT/data/sv_panel_to_accession_id.tsv

BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix
BGZIP=/home/tbellagio/miniforge3/envs/pang/bin/bgzip
PYBIN=/home/tbellagio/miniforge3/envs/hapfm/bin/python

WORK=$ROOT/scratch/p82_nobub
mkdir -p $WORK $ROOT/scratch/logs

echo "[$(date)] Step 1: Chr1 subset (NO vcfbub)"
$BCF view -r Chr1 $SRC -O u --threads 4 > $WORK/_step1_chr1.bcf
echo "  raw Chr1 records: $($BCF view -H $WORK/_step1_chr1.bcf | wc -l)"

echo ""
echo "[$(date)] Step 2: bcftools norm -m -any (atomize multi-allelics) | sort"
$BCF norm -m -any -f $REF --threads 4 -O u $WORK/_step1_chr1.bcf \
    | $BCF sort -m 16G -T $WORK/sort_tmp -O z -o $WORK/_step2_norm.vcf.gz
$TABIX -p vcf $WORK/_step2_norm.vcf.gz
rm -rf $WORK/sort_tmp $WORK/_step1_chr1.bcf
echo "  post-norm records: $($BCF view -H $WORK/_step2_norm.vcf.gz | wc -l)"

echo ""
echo "[$(date)] Step 3: rename Assembly_ID -> Accession_ID"
RENAME=$WORK/samples_82_rename.tsv
tail -n +2 $PANEL_TSV | awk -F'\t' '{print $1"\t"$3}' > $RENAME
$BCF reheader -s $RENAME $WORK/_step2_norm.vcf.gz > $WORK/p82_nobub_chr1.vcf.gz
$TABIX -p vcf $WORK/p82_nobub_chr1.vcf.gz
rm -f $WORK/_step2_norm.vcf.gz $WORK/_step2_norm.vcf.gz.tbi
echo "  samples: $($BCF query -l $WORK/p82_nobub_chr1.vcf.gz | wc -l)"

echo ""
echo "[$(date)] Step 4: build cn_var"
$PYBIN -u $ROOT/poolfreq/src/build_cn_var.py \
    --vcf $WORK/p82_nobub_chr1.vcf.gz --out $WORK/cn_var_p82_nobub --chrom Chr1

echo ""
echo "[$(date)] Step 5: extract ref_alt.tsv (chrom, pos, ref, alt)"
$BCF query -f '%CHROM\t%POS\t%REF\t%ALT\n' $WORK/p82_nobub_chr1.vcf.gz > $WORK/cn_var_p82_nobub.ref_alt.tsv

echo ""
echo "[$(date)] DONE"
echo "  cn_var size: $(ls -l $WORK/cn_var_p82_nobub.cn_var.npz | awk '{print $5}')"
echo "  ref_alt rows: $(wc -l < $WORK/cn_var_p82_nobub.ref_alt.tsv)"
ls -lh $WORK/

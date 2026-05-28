#!/bin/bash
#SBATCH --job-name=chr1_pg
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/A2_pg_%j.out
#SBATCH --error=logs/A2_pg_%j.err
mkdir -p logs
set -euo pipefail

# Phase 2 Job A2: PG-side full pipeline on full Chr1 (depends on A1).
# Steps: subset → transfer_id → convert-to-biallelic → fill-tags → haploidize (het→.)
# No V4 record filter; per-cell het→missing only.

cd /global/scratch/users/tbellg/kmate/arch3/chr1

BCF=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/gwas/bin/tabix
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CONVERT=/global/scratch/users/tbellg/kmate/external_tools/pangenie-tools/pipelines/run-from-callset/scripts/convert-to-biallelic.py
TRANSFER=/global/scratch/users/tbellg/kmate/scratch/arch3_test/transfer_id_annotation.py

CHR1_ANNOT=chr1_135_annotated.sorted.vcf.gz
BIAL_CATALOG=chr1_135_annotated_biallelic.sorted.vcf.gz
PG_RAW=/global/scratch/users/tbellg/kmate/pangenie_genotyping/data/v3qc_tmp/pangenie_153_raw.vcf.gz

[ -s $CHR1_ANNOT ] || { echo "ERROR: missing $CHR1_ANNOT (A1 not done)"; exit 1; }

echo "[$(date)] === Step 1: subset PG (all 153 samples) to full Chr1 ==="
PG_CHR1=pg_153_chr1.vcf.gz
if [ ! -s $PG_CHR1 ]; then
  $BCF view -r Chr1 $PG_RAW -Oz -o $PG_CHR1
  $TABIX -p vcf $PG_CHR1
fi
echo "  PG Chr1 records: $($BCF view -H $PG_CHR1 | wc -l)"
echo "  samples: $($BCF query -l $PG_CHR1 | wc -l)"

echo
echo "[$(date)] === Step 2: transfer INFO/ID from chr1 catalog ==="
PG_ANNOT=pg_153_chr1_annotated.vcf
$PY $TRANSFER --cactus $CHR1_ANNOT --pg $PG_CHR1 --out $PG_ANNOT

echo
echo "[$(date)] === Step 3: convert-to-biallelic (this is the slow step, ~30-60 min) ==="
PG_BIAL=pg_153_chr1_biallelic.vcf
$BCF view $PG_ANNOT 2>/dev/null | $PY $CONVERT $BIAL_CATALOG > $PG_BIAL 2> pg_convert.log
echo "  convert exit: $?"
echo "  records emitted: $(grep -vc '^#' $PG_BIAL)"

echo
echo "[$(date)] === Step 4: fill-tags ==="
PG_FILLED=pg_153_chr1_biallelic_filled.vcf.gz
$BCF +fill-tags $PG_BIAL --threads 4 -Oz -o $PG_FILLED -- -t AC,AN,AC_Het,F_MISSING
$TABIX -p vcf $PG_FILLED

echo
# NO V4 record filter. See haploidize step below for rationale.
echo "[$(date)] === Step 5: SKIP V4 filter (per-cell het→missing handles everything) ==="

echo
echo "[$(date)] === Step 6: haploidize (0/0→0, 1/1→1, het→. [MISSING], ./.→.) ==="
PG_HAP=pg_153_chr1_haploid.vcf.gz
N_HET_FLIPS_FILE=pg_153_chr1_het_flips.tsv
$BCF view --threads 2 -Ov $PG_FILLED | \
awk -v flips_file=$N_HET_FLIPS_FILE 'BEGIN{OFS="\t"; n_het_to_miss=0; n_total=0}
    /^##/ { print; next }
    /^#CHROM/ { print; next }
    {
        $9 = "GT"
        for (i=10; i<=NF; i++) {
            split($i, parts, ":"); g = parts[1]
            if (g == "0/0" || g == "0|0") { $i = "0" }
            else if (g == "1/1" || g == "1|1") { $i = "1" }
            else if (g == "./." || g == ".|." || g == ".") { $i = "." }
            else if (g == "0/1" || g == "0|1" || g == "1/0" || g == "1|0") {
                $i = ".";  # het → missing (Arouisse 2020 Arabidopsis precedent)
                n_het_to_miss++;
            }
            else { $i = "." }
            n_total++;
        }
        print
    }
    END { print "het_to_missing_flips\ttotal_cells" > flips_file; print n_het_to_miss"\t"n_total >> flips_file }
' | $BGZIP -@ 4 -c > $PG_HAP
$TABIX -p vcf $PG_HAP
echo "  het→missing flip count:"
cat $N_HET_FLIPS_FILE

echo
echo "[$(date)] DONE A2 (PG Chr1)"
ls -lh $PG_HAP

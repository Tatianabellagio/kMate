#!/bin/bash
#SBATCH --job-name=haplo_vcf
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=logs/haplo_%j.out
#SBATCH --error=logs/haplo_%j.err

# =============================================================================
# haploidize_merged_vcf.sh
#
# Produce a fully-haploid sibling of the production merged VCF. Carrier-status
# convention: any-alt → 1, hom_ref → 0, missing → '.'. The original mixed-
# ploidy VCF is untouched on disk.
#
# Pipeline:
#   1. bcftools norm -m -any   → split multi-allelic records into biallelic
#                                so cactus integer alleles (33, 51, …) and
#                                PanGenie 1/2 calls all collapse cleanly.
#   2. awk on the GT field     → cactus haploid stays {0,1,.}; PanGenie
#                                {0/0,0/1,1/1,./.} → {0,1,1,.}.
#   3. bgzip + tabix.
#
# Why: PanGenie het rate is 0.79% of all calls and 8.4% of alt-carrying calls
# (Chr1:1-2Mb sample, 2026-05-08). 92% of carriers are already 1/1; for inbred
# A. thaliana the residual het is dominated by PanGenie call noise on noisy
# Cao 2011 GAII libraries (RESULTS_LOG 2026-05-03). Carrier-status haploid is
# the convention build_cn_var.py already uses, so the haploid VCF makes the
# on-disk catalog match what every downstream tool already does in memory.
# =============================================================================
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
SRC=$BASE/data/merged/founders_231_chr.vcf.gz
OUT=$BASE/data/merged/founders_231_chr.haploid.vcf.gz

BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix

[ -s "$SRC" ] || { echo "ERROR: missing $SRC" >&2; exit 1; }
[ ! -s "$OUT" ] || { echo "[$(date)] $OUT already exists — exiting (delete to rerun)"; exit 0; }

echo "[$(date)] decompose multi-allelics → haploidize → bgzip"
echo "  src: $SRC"
echo "  out: $OUT"

$BCF norm -m -any --threads 4 -Ou "$SRC" 2>/dev/null | \
$BCF view --threads 2 -Ov - | \
awk 'BEGIN{OFS="\t"}
    /^##/ { print; next }
    /^#CHROM/ { print; next }
    {
        # FORMAT field after merge is "GT:GQ:GL:KC" (cactus side fills
        # GQ/GL/KC with `.`). For the haploid sibling we keep only GT.
        $9 = "GT"
        for (i=10; i<=NF; i++) {
            split($i, parts, ":")
            g = parts[1]
            # carrier-status haploid map
            if (g == "0" || g == "0/0" || g == "0|0") {
                $i = "0"
            } else if (g == "." || g == "./." || g == ".|.") {
                $i = "."
            } else {
                # any allele containing a non-zero alt index → carrier=1
                # covers haploid 1, diploid 0/1, 1/0, 0|1, 1|0, 1/1, 1|1.
                # bcftools norm -m -any has already split multi-allelics so
                # 1/2 etc. cannot appear here.
                $i = "1"
            }
        }
        print
    }' | \
/global/home/users/tbellg/miniforge3/envs/pang/bin/bgzip -@ 4 -c > "$OUT"

$TABIX -p vcf "$OUT"

N_REC=$($BCF view -H "$OUT" 2>/dev/null | wc -l)
N_SAM=$($BCF query -l "$OUT" 2>/dev/null | wc -l)
echo "[$(date)] DONE  records=$N_REC  samples=$N_SAM"

echo ""
echo "=== sanity check on first 3 records (should be all 0/1/. — no /) ==="
$BCF view -H "$OUT" 2>/dev/null | head -3 | awk '{
    for (i=10; i<=NF && i<=14; i++) printf "%s ", $i
    print ""
}'

ls -lh "$OUT"

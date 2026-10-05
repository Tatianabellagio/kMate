#!/bin/bash
#SBATCH --job-name=decompose_check
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
# Check `kmate decompose` against the arch3 panel build on a Chr1 slice.
#
#   --gfa route:      Minigraph-Cactus VCF + GFA -> decompose          (one command)
#   reference route:  arch3's annotated multi-allelic Chr1 + catalog -> convert-to-biallelic -> sort
#
# Both must give identical records, IDs and genotypes, and build-var-pa must accept
# the output. 2026-10-05, Chr1:1-300000: 12,719 multi-allelic -> 19,725 biallelic
# records, 135 samples, identical; decompose --gfa 3 min / 7.7 GB.
#
# Usage: sbatch tests/run_decompose_check.sh   (or run inside an allocation)
#        WORK=<dir> REGION=Chr1:1-300000 can be overridden.
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PYTHONPATH=$ROOT/src
PANG=${PANG:-/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output}
A3=${A3:-$ROOT/panel/arch3/chr1}            # chr1_135_annotated{,_biallelic}.sorted.vcf.gz
REGION=${REGION:-Chr1:1-300000}
WORK=${WORK:-$(mktemp -d)}
cd "$WORK"; echo "work dir: $WORK"

bcftools view -r $REGION $PANG/pang_1001gplus_all.vcf.gz -Oz -o slice.vcf.gz
bcftools index -t -f slice.vcf.gz

python -m kmate.cli decompose --gfa $PANG/pang_1001gplus_all.gfa.gz \
    --genotyped-vcf slice.vcf.gz --out new.vcf.gz --threads 4

bcftools view -r $REGION $A3/chr1_135_annotated.sorted.vcf.gz \
  | python $ROOT/src/kmate/_vendor/pangenie/convert-to-biallelic.py \
      $A3/chr1_135_annotated_biallelic.sorted.vcf.gz > ref.unsorted.vcf
bcftools sort -Oz -o ref.vcf.gz ref.unsorted.vcf

q='%CHROM\t%POS\t%REF\t%ALT\t%INFO/ID[\t%GT]\n'
bcftools query -f "$q" new.vcf.gz | sort > new.q
bcftools query -f "$q" ref.vcf.gz | sort > ref.q
echo "records: new=$(wc -l < new.q) ref=$(wc -l < ref.q)"
cmp new.q ref.q && echo "PASS: records + genotypes identical"

python -m kmate.cli build-var-pa --vcf new.vcf.gz --chrom ${REGION%%:*} --out varpa > varpa.log
echo "PASS: build-var-pa accepted the output"

#!/bin/bash
#SBATCH --job-name=archA_pg
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/jobA_pg_%j.out
#SBATCH --error=logs/jobA_pg_%j.err
mkdir -p logs
set -euo pipefail

# Job A: full PG-side Arch 3 pipeline on test region (Chr1:5.8M-14M), all 153 samples.
# Steps: subset → transfer ID → convert-to-biallelic → fill-tags → V4 filter → haploidize.

cd /global/scratch/users/tbellg/hapfire_sv/scratch/arch3_test

BCF=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/gwas/bin/tabix
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CONVERT=/global/scratch/users/tbellg/hapfire_sv/external_tools/pangenie-tools/pipelines/run-from-callset/scripts/convert-to-biallelic.py
TRANSFER=/global/scratch/users/tbellg/hapfire_sv/scratch/arch3_test/transfer_id_annotation.py

CACTUS_ANNOT=full135_test_annotated.sorted.vcf.gz       # from job 63045
BIAL_CATALOG=full135_test_annotated_biallelic.sorted.vcf.gz
PG_RAW=/global/scratch/users/tbellg/hapfire_sv/pangenie_genotyping/data/v3qc_tmp/pangenie_153_raw.vcf.gz

[ -s $CACTUS_ANNOT ] || { echo "ERROR: missing $CACTUS_ANNOT"; exit 1; }

echo "[$(date)] === Step 1: subset PG (all 153 samples) to test region ==="
PG_TEST=pg_153_test_region.vcf.gz
if [ ! -s $PG_TEST ]; then
  $BCF view -r Chr1:5800000-14000000 $PG_RAW -Oz -o $PG_TEST
  $TABIX -p vcf $PG_TEST
fi
echo "  PG 153-sample, test region records: $($BCF view -H $PG_TEST | wc -l)"
echo "  samples: $($BCF query -l $PG_TEST | wc -l)"

echo
echo "[$(date)] === Step 2: transfer INFO/ID from cactus catalog ==="
PG_ANNOT=pg_153_test_annotated.vcf
$PY $TRANSFER --cactus $CACTUS_ANNOT --pg $PG_TEST --out $PG_ANNOT

echo
echo "[$(date)] === Step 3: convert-to-biallelic ==="
PG_BIAL=pg_153_test_biallelic.vcf
$BCF view $PG_ANNOT 2>/dev/null | $PY $CONVERT $BIAL_CATALOG > $PG_BIAL 2> pg_convert.log
echo "  convert exit: $?"
echo "  records emitted: $(grep -vc '^#' $PG_BIAL)"

echo
echo "[$(date)] === Step 4: fill-tags (AC, AN, AC_Het, F_MISSING) ==="
PG_FILLED=pg_153_test_biallelic_filled.vcf.gz
$BCF +fill-tags $PG_BIAL --threads 4 -Oz -o $PG_FILLED -- -t AC,AN,AC_Het,F_MISSING
$TABIX -p vcf $PG_FILLED

echo
# ---------------------------------------------------------------------------
# Step 5: het handling design decision (2026-05-20).
#
# We intentionally do NOT apply a V4 (AC_Het / AC) record-level filter, even
# though xwu's GrENE-Net pipeline does (with V4 < 0.01).
#
# Why no V4 record drop:
#   - xwu's V4 < 0.01 threshold is an unpublished in-house convention; not
#     documented in hapFIRE README, grenepipe paper (Czech 2022 Bioinformatics),
#     or GrENE-Net Science 2026.
#   - LOO concordance analysis showed hets are concentrated at SPECIFIC sites
#     (paralogue mis-mapping per Jaegle/Nordborg 2023 Genome Biology), not
#     spread randomly. Dropping the entire record loses VALID hom calls from
#     other founders that don't have hets at that site.
#   - Per-cell het masking (next step) handles the bad cells while preserving
#     all the good ones in the record.
#   - All records preserved enables maximum information for downstream pool-seq
#     AF estimation and ecotype identification.
#
# Het handling lives entirely in the haploidize step below.
# ---------------------------------------------------------------------------
echo "[$(date)] === Step 5: SKIP V4 record filter (per-cell het→missing handles everything) ==="
PG_V4=$PG_FILLED  # passthrough — no row filtering
echo "  no V4 filter applied; all records retained ($($BCF view -H $PG_V4 | wc -l) records)"

echo
# ---------------------------------------------------------------------------
# Step 6: diploid → haploid conversion + het handling.
#
# Genotype mapping (applied per-cell to each diploid PG GT):
#   0/0 (REF hom)            → 0      (REF call)
#   1/1 (ALT hom)             → 1      (ALT carrier)
#   ./.  (missing)            → .      (missing — preserved)
#   0/1, 1/0 (heterozygous)   → .      (HET → MISSING)  ← see rationale
#   anything else (defensive) → .      (treat as missing rather than guess)
#
# Why hets → MISSING (not REF, not ALT, not kept):
#   1. *A. thaliana* is highly selfing — expected genome-wide het ~1%. Most
#      0/1 calls in 1001G short-read PG output are sequencing/mapping
#      artefacts (~3-5% per accession per LOO; Jaegle 2023 shows ~44% of all
#      called SNP records have at least one het call, concentrated at
#      paralogue/CNV sites).
#   2. In haploid representation a "0/1" cell cannot exist — must collapse
#      to either 0, 1, or "." (missing). Each carries a different epistemic
#      claim:
#        het → 0 (REF):  asserts "this founder is REF" — biased toward 0
#                        when residual hets are real (e.g. recent outcrossing)
#        het → 1 (ALT):  asserts "this founder is ALT carrier" — biases
#                        AF estimate upward; not justified by evidence
#        het → .  (missing):  honest — "we don't know what this founder
#                             actually is at this position"
#   3. Arouisse 2020 (Plant J., Arabidopsis 1001G + RegMap + Beagle):
#      direct published precedent for het → NA on this same data class.
#      Their next step was Beagle imputation; we don't impute, so cells
#      stay as missing (cn_var_called marks them as not-called, AN
#      reduced by 1 cell per masked het, AC unchanged).
#   4. cn_full is sensitive to "." vs "0":
#        "." cell → build_kmer_cn.py writes "N" into the founder consensus
#                  FASTA (current behavior; matches v3qc_v3-era fix). K-mers
#                  spanning the N position get dropped → founder contributes
#                  NO k-mer evidence at that locus.
#        "0" cell → consensus has TAIR10 REF at the position. K-mers include
#                  the REF base; founder contributes evidence as if confidently
#                  REF (potentially overclaiming when truth is unknown).
#      So het→missing is even MORE conservative than the cn_var_called
#      story alone: founders we're not sure about also drop out of cn_full
#      at those positions, not just out of the AC/AN denominator.
#      (TODO: verify build_kmer_cn.py current default; see memory entry
#       project_v3qcv2_treat_missing_as_n_bug for history.)
#
# Filter combination (final):
#   GQ < 20 cell mask          (applied UPSTREAM in pangenie_153_raw)
#   per-cell het → missing      (THIS STEP)
#   NO per-side V4/AC=0 drop
#   post-merge AN=0 only       (Job B, after cactus+PG merge)
# ---------------------------------------------------------------------------
echo "[$(date)] === Step 6: haploidize (0/0→0, 1/1→1, het→. [MISSING], ./.→.) ==="
PG_HAP=pg_153_test_haploid.vcf.gz
N_HET_FLIPS_FILE=pg_153_het_flips.tsv
$BCF view --threads 2 -Ov $PG_V4 | \
awk -v flips_file=$N_HET_FLIPS_FILE 'BEGIN{OFS="\t"; n_het_to_miss=0; n_total=0}
    /^##/ { print; next }
    /^#CHROM/ { print; next }
    {
        $9 = "GT"
        for (i=10; i<=NF; i++) {
            split($i, parts, ":")
            g = parts[1]
            if (g == "0/0" || g == "0|0") { $i = "0" }
            else if (g == "1/1" || g == "1|1") { $i = "1" }
            else if (g == "./." || g == ".|." || g == ".") { $i = "." }
            else if (g == "0/1" || g == "0|1" || g == "1/0" || g == "1|0") {
                $i = ".";   # het → MISSING (Arouisse 2020 precedent)
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
# Spot-check via Python — avoids SIGPIPE issues from awk-exit-early in pipe under `set -e -o pipefail`.
echo "[$(date)] === Step 7+8: spot-check 3 positions on 10 samples + panel-level AC ==="
$PY <<PYEOF
import gzip
SAMPLES = ['100001', '9596', '9965', '6209', '5151', '10011', '6013', '9595', '6025', '265']
SPOT = [('5870018','T','A'), ('10421645','T','C'),
        ('13843898','C','T'), ('13843898','CT','TC'), ('13843898','CT','TG')]
results = {}
samples = None
sample_idx = {}
with gzip.open('$PG_HAP', 'rt') as f:
    for line in f:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            samples = line.rstrip().split('\t')[9:]
            sample_idx = {s:i for i,s in enumerate(samples)}
            continue
        p = line.rstrip().split('\t')
        if (p[1], p[3], p[4]) in SPOT:
            results[(p[1], p[3], p[4])] = p
print(f'PG haploid samples: {len(samples)}')
for k in SPOT:
    if k not in results:
        print(f'\n--- Chr1:{k[0]} {k[1]}>{k[2]} ---  RECORD NOT FOUND')
        continue
    r = results[k]
    gts = r[9:]
    ac = sum(1 for g in gts if g == '1')
    an = sum(1 for g in gts if g in ('0','1'))
    miss = sum(1 for g in gts if g == '.')
    print(f'\n--- Chr1:{k[0]} {k[1]}>{k[2]} ---')
    print(f'  panel: AC={ac}/{an}  AF={ac/max(an,1):.4f}  missing={miss}/{len(gts)}')
    for s in SAMPLES:
        if s in sample_idx:
            print(f'    {s}: GT={gts[sample_idx[s]]}')
PYEOF

echo
echo "[$(date)] DONE Job A"
ls -lh pg_153_test_haploid.vcf.gz pg_153_het_flips.tsv

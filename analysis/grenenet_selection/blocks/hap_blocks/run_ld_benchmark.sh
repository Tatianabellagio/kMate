#!/bin/bash
# Block-size error benchmark on the equimolar g0-231 sim (all 231 founders).
# kMate per-block EM (local-only, global-free) at each r2 partition vs global
# baseline; score per-record AF vs recomb_truth. Count k-mers ONCE (--kmer-db).
# Uses kMate's OWN panel + self-contained var_pa LD blocks (no external panel).
source ~/.bashrc 2>/dev/null || true
mamba activate kmate 2>/dev/null || true
set -euo pipefail
echo "python=$(which python)"
ROOT=/global/scratch/users/tbellg/kmate
HB=$ROOT/analysis/grenenet_selection/hap_blocks
SIM=$ROOT/benchmarks/p231/sims/cov10_n231_g0_s42_hotspots_p231_chr1
KMER=$ROOT/benchmarks/p231/data/kmer_pa_p231_filt2inv/kmer_pa
VARPA=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
VARMETA=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
OUT=$HB/bench_g0_231; mkdir -p "$OUT"
DB=$OUT/pool.jf
PART_TAG=${1:-alltype}   # which partition set: alltype (ld_blocks_r2_*) or snp (ld_blocks_snp_r2_*)
PREFIX=$HB/ld_blocks
[ "$PART_TAG" = "snp" ] && PREFIX=$HB/ld_blocks_snp

echo "[$(date)] build k-mer DB once"
[ -s "$DB" ] || python -m kmate.build_kmer_db --reads $SIM/reads/r1.fq $SIM/reads/r2.fq --out "$DB" --threads 8

echo "[$(date)] GLOBAL baseline (production, all 231 founders)"
python -m kmate.per_sample_per_chrom --kmer-pa-prefix "$KMER" --var-pa "$VARPA" --var-meta "$VARMETA" \
    --kmer-db "$DB" --sample global --out "$OUT/global.tsv" --threads 8 --chroms Chr1 \
    --block-mode global --kmer-weight uniform --normalize per_founder

for r2 in 0.10 0.20 0.30 0.40; do
    BT=${PREFIX}_r2_${r2}.tsv
    [ -s "$BT" ] || { echo "SKIP r2=$r2 (no $BT)"; continue; }
    echo "[$(date)] block-mode window local-only  r2=$r2  ($PART_TAG)"
    python -m kmate.per_sample_per_chrom --kmer-pa-prefix "$KMER" --var-pa "$VARPA" --var-meta "$VARMETA" \
        --kmer-db "$DB" --sample r2_${r2} --out "$OUT/${PART_TAG}_r2_${r2}.tsv" --threads 8 --chroms Chr1 \
        --block-mode window --blocks-tsv "$BT" --local-only --min-kmers-per-block 1 \
        --kmer-weight uniform --normalize per_founder
done
echo "[$(date)] DONE"

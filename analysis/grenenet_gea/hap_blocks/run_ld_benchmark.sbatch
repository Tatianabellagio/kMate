#!/bin/bash
#SBATCH --job-name=ld_bench
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/hap_blocks/logs/ld_bench_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/hap_blocks/logs/ld_bench_%j.out
# Block-size error benchmark on equimolar g0-231 sim. kMate per-block EM
# (local-only, global-free) at each r2 partition vs global baseline; score
# per-record AF vs recomb_truth. Self-contained var_pa LD blocks (audited).
mkdir -p /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/hap_blocks/logs
source ~/.bashrc 2>/dev/null || true
mamba activate kmate 2>/dev/null || true
set -uo pipefail
echo "python=$(which python)  host=$(hostname)"
ROOT=/global/scratch/users/tbellg/kmate
HB=$ROOT/analysis/grenenet_gea/hap_blocks
SIM=$ROOT/benchmarks/p231/sims/cov10_n231_g0_s42_hotspots_p231_chr1
KMER=$ROOT/benchmarks/p231/data/kmer_pa_p231_filt2inv/kmer_pa
VARPA=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
VARMETA=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
OUT=$HB/bench_g0_231; mkdir -p "$OUT"
DB=$OUT/pool.jf
run() { echo "[$(date)] $*"; "$@"; echo "  exit=$?"; }

echo "[$(date)] build k-mer DB once (hash 1G -- cov10 chr1 is small)"
[ -s "$DB" ] || python -m kmate.build_kmer_db --reads $SIM/reads/r1.fq $SIM/reads/r2.fq \
    --out "$DB" --threads 8 --hash-size 1G

echo "[$(date)] GLOBAL baseline (production, all 231 founders)"
python -m kmate.per_sample_per_chrom --kmer-pa-prefix "$KMER" --var-pa "$VARPA" --var-meta "$VARMETA" \
    --reads $SIM/reads/r1.fq $SIM/reads/r2.fq --kmer-db "$DB" --sample global --out "$OUT/global.tsv" --threads 8 --chroms Chr1 \
    --block-mode global --kmer-weight uniform --normalize per_founder

for r2 in 0.10 0.20 0.30 0.40; do
    BT=$HB/ld_blocks_r2_${r2}.tsv
    [ -s "$BT" ] || { echo "SKIP r2=$r2 (no $BT)"; continue; }
    echo "[$(date)] block-mode window local-only  r2=$r2"
    python -m kmate.per_sample_per_chrom --kmer-pa-prefix "$KMER" --var-pa "$VARPA" --var-meta "$VARMETA" \
        --reads $SIM/reads/r1.fq $SIM/reads/r2.fq --kmer-db "$DB" --sample r2_${r2} --out "$OUT/alltype_r2_${r2}.tsv" --threads 8 --chroms Chr1 \
        --block-mode window --blocks-tsv "$BT" --local-only --min-kmers-per-block 1 \
        --kmer-weight uniform --normalize per_founder
done
echo "[$(date)] ALL DONE"; ls -lh "$OUT"/*.tsv

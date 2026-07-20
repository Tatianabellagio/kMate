#!/bin/bash
#SBATCH --job-name=e2e_hb
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/hap_blocks/logs/e2e_hb_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/hap_blocks/logs/e2e_hb_%j.out
# End-to-end validation of the haploblock-collapse code path (commit b4d6ce0):
# run per_sample_per_chrom GLOBAL mode on the g0-231 equimolar sim (Chr1), reusing
# the prebuilt pool.jf DB, then confirm the new code reproduces the pre-change
# bench global-h EXACTLY (K_b=231 => collapse is a byte-identical no-op).
mkdir -p /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/hap_blocks/logs
source ~/.bashrc 2>/dev/null || true
mamba activate kmate 2>/dev/null || true
set -uo pipefail
echo "host=$(hostname)  python=$(which python)  $(date)"

ROOT=/global/scratch/users/tbellg/kmate
KMER=$ROOT/benchmarks/p231/data/kmer_pa_p231_filt2inv/kmer_pa
VARPA=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
VARMETA=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
DB=$ROOT/analysis/grenenet_gea/hap_blocks/bench_g0_231/pool.jf
SIM=$ROOT/benchmarks/p231/sims/cov10_n231_g0_s42_hotspots_p231_chr1
E2E=$ROOT/analysis/grenenet_gea/hap_blocks/e2e_newcode
mkdir -p "$E2E"

echo "[$(date)] running global-mode estimation (new haploblock code) ..."
python -m kmate.per_sample_per_chrom --kmer-pa-prefix "$KMER" \
  --var-pa "$VARPA" --var-meta "$VARMETA" \
  --reads "$SIM/reads/r1.fq" "$SIM/reads/r2.fq" --kmer-db "$DB" \
  --sample global_newcode --out "$E2E/global_newcode.tsv" \
  --threads 8 --chroms Chr1 --block-mode global \
  --kmer-weight uniform --normalize per_founder

echo "[$(date)] comparing new h to the pre-change bench global-h (expect ~0 diff) ..."
python - <<'PY'
import numpy as np
ROOT="/global/scratch/users/tbellg/kmate"
new=np.load(f"{ROOT}/analysis/grenenet_gea/hap_blocks/e2e_newcode/global_newcode.h_per_chrom.npz",allow_pickle=True)
old=np.load(f"{ROOT}/analysis/grenenet_gea/hap_blocks/bench_g0_231/global.h_per_chrom.npz",allow_pickle=True)
hn=new["Chr1"].astype(float); ho=old["Chr1"].astype(float)
assert (new["founders"].astype(str)==old["founders"].astype(str)).all(), "founder order mismatch"
d=np.abs(hn-ho)
print(f"F={len(hn)}  new sum={hn.sum():.6f}  old sum={ho.sum():.6f}")
print(f"max|new-old|={d.max():.3e}  mean|new-old|={d.mean():.3e}")
print(f"corr={np.corrcoef(hn,ho)[0,1]:.8f}")
print("VERDICT:", "NO-OP OK (haploblock collapse reproduces global at K_b=231)"
      if d.max()<1e-4 else "DIFFERS -- investigate")
PY
echo "[$(date)] DONE"

#!/bin/bash
#SBATCH --job-name=p231_a3c_filt2inv
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --output=logs/03c_filt2inv_%j.out
#SBATCH --error=logs/03c_filt2inv_%j.out

# benchmarks/p231 Phase A3c -- filt2inv (drop ac<2 singletons AND ac==F
# invariants), matching the REAL production K_pa filter (data/kmer_pa_231_arch3_
# filt2inv; PIPELINE_STATE.md Sec.0). The benchmark's existing "filt2" arm
# (03b_build_kmer_pa_filt2_p231.sh) only drops singletons -- NOT apples-to-apples
# with production, which is why absorbed-founder counts here don't match the
# 0-absorbed result from the controlled validation notebook (which used the
# real filt2inv panel). 2026-07-07.
mkdir -p logs
set -euo pipefail
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p231
SRC=$CTRL/data/kmer_pa_p231
OUT=$CTRL/data/kmer_pa_p231_filt2inv
mkdir -p $OUT

if [ -s "$OUT/kmer_pa_Chr1.kmer_pa.npz" ] && [ -s "$OUT/kmer_pa_Chr1.meta.npz" ]; then
    echo "[$(date)] kmer_pa_p231_filt2inv already present -- skip"; ls -lh $OUT/kmer_pa_Chr1.*; exit 0
fi
echo "[$(date)] filt2inv (2<=ac<=F-1) build from $SRC -> $OUT"

$PY << EOF
import numpy as np
from scipy.sparse import load_npz, save_npz
from pathlib import Path
chrom = "Chr1"; SRC = Path("${SRC}"); OUT = Path("${OUT}"); OUT.mkdir(exist_ok=True)
kmer_pa = load_npz(SRC / f"kmer_pa_{chrom}.kmer_pa.npz")
meta = np.load(SRC / f"kmer_pa_{chrom}.meta.npz", allow_pickle=True)
F, K = kmer_pa.shape
print(f"[{chrom}] kmer_pa ({F},{K:,}) nnz={kmer_pa.nnz:,}", flush=True)
ac = np.asarray(kmer_pa.sum(axis=0)).flatten().astype(np.int32)
keep = (ac >= 2) & (ac <= F - 1)
print(f"[{chrom}] keep 2<=ac<=F-1: {keep.sum():,}/{K:,} "
      f"(ac=0:{(ac==0).sum():,} ac=1:{(ac==1).sum():,} ac=F:{(ac==F).sum():,})", flush=True)
cn_f = kmer_pa.tocsc()[:, keep].tocsr()
new_meta = {}
for k in meta.keys():
    a = meta[k]
    new_meta[k] = a[keep] if (a.shape and a.shape[0] == K) else a
save_npz(OUT / f"kmer_pa_{chrom}.kmer_pa.npz", cn_f)
np.savez(OUT / f"kmer_pa_{chrom}.meta.npz", **new_meta)
print(f"[{chrom}] DONE filt2inv nnz={cn_f.nnz:,} shape={cn_f.shape}", flush=True)
EOF
echo "[$(date)] DONE"; ls -lh $OUT/kmer_pa_Chr1.kmer_pa.npz $OUT/kmer_pa_Chr1.meta.npz

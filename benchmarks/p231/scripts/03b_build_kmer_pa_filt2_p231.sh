#!/bin/bash
#SBATCH --job-name=p231_a3b_filt2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/03b_filt2_%j.out
#SBATCH --error=logs/03b_filt2_%j.err

# benchmarks/p231 Phase A3b -- filt2 (drop ac<2 singleton k-mers), mirrors
# benchmarks/p80/03b_build_kmer_pa_filt2_p80.sh. This is the front-runner kmer_pa base.
mkdir -p logs
set -uo pipefail
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p231
SRC=$CTRL/data/kmer_pa_p231
OUT=$CTRL/data/kmer_pa_p231_filt2
mkdir -p $OUT

if [ -s "$OUT/kmer_pa_Chr1.kmer_pa.npz" ] && [ -s "$OUT/kmer_pa_Chr1.meta.npz" ]; then
    echo "[$(date)] kmer_pa_p231_filt2 already present -- skip"; ls -lh $OUT/kmer_pa_Chr1.*; exit 0
fi
echo "[$(date)] filt2 (ac>=2) build from $SRC -> $OUT"

$PY << EOF
import numpy as np
from scipy.sparse import load_npz, save_npz
from pathlib import Path
chrom="Chr1"; SRC=Path("${SRC}"); OUT=Path("${OUT}"); OUT.mkdir(exist_ok=True)
kmer_pa=load_npz(SRC/f"cn_{chrom}.kmer_pa.npz"); meta=np.load(SRC/f"cn_{chrom}.meta.npz",allow_pickle=True)
F,K=kmer_pa.shape; print(f"[{chrom}] kmer_pa ({F},{K:,}) nnz={kmer_pa.nnz:,}",flush=True)
ac=np.asarray(kmer_pa.sum(axis=0)).flatten().astype(np.int32); keep=ac>=2
print(f"[{chrom}] keep ac>=2: {keep.sum():,}/{K:,} (ac=0:{(ac==0).sum():,} ac=1:{(ac==1).sum():,})",flush=True)
cn_f=kmer_pa.tocsc()[:,keep].tocsr()
new_meta={}
for k in meta.keys():
    a=meta[k]
    new_meta[k]=a[keep] if (a.shape and a.shape[0]==K) else a
save_npz(OUT/f"cn_{chrom}.kmer_pa.npz",cn_f); np.savez(OUT/f"cn_{chrom}.meta.npz",**new_meta)
print(f"[{chrom}] DONE filt2 nnz={cn_f.nnz:,} shape={cn_f.shape}",flush=True)
EOF
echo "[$(date)] DONE"; ls -lh $OUT/kmer_pa_Chr1.kmer_pa.npz $OUT/kmer_pa_Chr1.meta.npz

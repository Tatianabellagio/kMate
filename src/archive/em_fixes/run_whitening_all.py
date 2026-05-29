"""
Run whitening IRLS/GLS over all 13 test samples (5 g0 + 8 SEEDMIX) loading the
cn matrix + bubble sizes ONCE. Writes per-sample h to scratch/h_fixes/whitening/.
"""
from __future__ import annotations
import os, sys, time, gc
import numpy as np
from scipy.sparse import load_npz

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from whitening import whitened_irls

CN_PREFIX = "data/cn_full_231_v3qc_v3_filt2/cn_Chr1"
OUTDIR = "scratch/h_fixes/whitening"

G0_SIMS = [
    "g0_n231_rep0_rand", "g0_n200_rep0_rand",
    "g0_n50_rep0_cact", "g0_n50_rep1_bal", "g0_n50_rep2_pg",
]
SEEDMIX = [f"S{i}" for i in range(1, 9)]


def counts_path(sample):
    if sample.startswith("g0_"):
        return f"scratch/g0_sweep_h_test/filt2_{sample}.counts.npy"
    else:
        return f"scratch/seedmix_h_test/filt2_{sample}.counts.npy"


def main():
    t0 = time.time()
    os.makedirs(OUTDIR, exist_ok=True)
    print("Loading cn + meta (once) ...", flush=True)
    cn = load_npz(CN_PREFIX + ".cn.npz").tocsr()
    meta = np.load(CN_PREFIX + ".meta.npz", allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    bubble_id = np.asarray(meta["bubble_id"]).astype(np.int64)
    F, K = cn.shape
    mb_full = np.bincount(bubble_id, minlength=bubble_id.max() + 1)[bubble_id].astype(np.float64)
    print(f"  cn F={F} K={K:,} nnz={cn.nnz:,}  [{time.time()-t0:.0f}s]", flush=True)
    del meta, bubble_id
    gc.collect()

    samples = G0_SIMS + SEEDMIX
    for sample in samples:
        cp = counts_path(sample)
        if not os.path.exists(cp):
            print(f"!! MISSING counts for {sample}: {cp}", flush=True)
            continue
        ts = time.time()
        counts = np.load(cp)
        assert counts.shape[0] == K, f"{sample}: counts len {counts.shape[0]} != K {K}"
        nz = counts > 0
        cn_nz = cn[:, nz].tocsr()
        cn_nz.sort_indices()
        counts_nz = counts[nz].astype(np.float64)
        mb_nz = mb_full[nz]
        print(f"[{sample}] nz={nz.sum():,} sum={counts.sum():,}  building done [{time.time()-ts:.0f}s]",
              flush=True)
        h, info = whitened_irls(
            counts_nz, None, cn_nz, mb_nz, F,
            max_iter=200, inner_iter=5, tol=1e-7, verbose=True,
        )
        out = os.path.join(OUTDIR, f"{sample}.whiten.npz")
        np.savez(out, founders=founders, h=h.astype(np.float64), sample=sample,
                 iterations=info["iterations"], converged=info["converged"],
                 lam=info["lambda"], delta_history=np.array(info["delta_history"]))
        print(f"[{sample}] DONE iters={info['iterations']} conv={info['converged']} "
              f"lam={info['lambda']:.4g} -> {out}  [{time.time()-ts:.0f}s]\n", flush=True)
        del cn_nz, counts, counts_nz, mb_nz
        gc.collect()

    print(f"ALL DONE  [TOTAL {time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()

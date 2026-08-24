"""NEW GLOBAL on the real seed-mix (Chr1): r2=0.1 blocks -> exact-identity
haploblocks (K_b) -> per-block EM (production solve_em, full-panel Kf_w) ->
K_b-collapse + equal-split -> aggregate to a genome-wide (Chr1) per-founder h.

Compared against OLD GLOBAL = chromosome-wise EM with the correct Kf_w
(perchrom_new[0] in persample_SEEDMIX_S*.npz, already computed by
run_persample_kfw.py). Same cached counts, same panel, same solver, same
normalization -- the ONLY thing that changes is block+haploblock vs whole-chrom.

Real seed-mix has no hard ground truth; intended mix is ~equimolar (1/231), so
we report absorption (founders pushed to ~0) and whether new global rescues the
founders old global absorbs.
"""
import os, sys, time
import numpy as np
import scipy.sparse as sp

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/qc/seedmix_validation/fix_kfw_fullpanel"
CDIR = f"{ROOT}/analysis/grenenet_gea/qc/seedmix_validation/fix_seedmix/counts"
PANEL = f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
BLOCKTSV = f"{ROOT}/analysis/grenenet_gea/blocks/hap_blocks/ld_blocks_r2_0.10.tsv"
CHROM = "Chr1"
F = 231
sys.path.insert(0, f"{ROOT}/src")
from kmate.em_solver import solve_em

SAMPLES = ["SEEDMIX_S1", "SEEDMIX_S2", "SEEDMIX_S3", "SEEDMIX_S4", "SEEDMIX_S6", "SEEDMIX_S8"]


def load_blocks():
    rows = []
    with open(BLOCKTSV) as fh:
        next(fh)
        for ln in fh:
            c, s, e, n = ln.split()
            if c == CHROM:
                rows.append((int(s), int(e)))
    return rows


def main():
    t0 = time.time()
    meta = np.load(f"{PANEL}.meta.npz", allow_pickle=True)
    founders = meta["founders"].astype(str)
    bub = meta["bubble_id"].astype(np.int64)                       # per-kmer bubble
    bs = meta["bubble_start"].astype(np.int64); be = meta["bubble_end"].astype(np.int64)
    kmer_pos = ((bs[bub] + be[bub]) // 2)                          # per-kmer midpoint
    K = sp.load_npz(f"{PANEL}.kmer_pa.npz").astype(np.float32).tocsc()
    print(f"panel Chr1: {K.shape[0]} founders x {K.shape[1]:,} k-mers ({time.time()-t0:.0f}s)", flush=True)

    blocks = load_blocks()
    # assign k-mers to blocks by midpoint (sorted)
    order = np.argsort(kmer_pos); spos = kmer_pos[order]
    block_cols = []
    for (a, b) in blocks:
        lo = np.searchsorted(spos, a, "left"); hi = np.searchsorted(spos, b, "right")
        block_cols.append(order[lo:hi])
    nk = np.array([len(c) for c in block_cols])
    print(f"{len(blocks)} r2=0.1 blocks on Chr1; k-mers/block: "
          f"min={nk.min()} median={int(np.median(nk))} max={nk.max()}", flush=True)

    # exact-identity haploblock labels per block (panel-defined, sample-independent)
    labels, Kbs = [], []
    for cols in block_cols:
        if len(cols) == 0:
            labels.append(None); Kbs.append(0); continue
        sub = K[:, cols].tocsr()
        sig = [hash(sub.indices[sub.indptr[i]:sub.indptr[i + 1]].tobytes()) for i in range(F)]
        _, lab = np.unique(sig, return_inverse=True)
        labels.append(lab); Kbs.append(lab.max() + 1)
    print(f"K_b per block: {Kbs} ({time.time()-t0:.0f}s)", flush=True)

    results = {}
    for sid in SAMPLES:
        c = np.load(f"{CDIR}/{sid}_{CHROM}.npy").astype(np.float32)
        clean_blocks = []
        for j, cols in enumerate(block_cols):
            if labels[j] is None:
                continue
            Ksub = np.asarray(K[:, cols].todense(), dtype=np.float32)   # F x nk_block
            kfw_b = Ksub.sum(axis=1)                                    # full-block Kf_w (Kf_w fix)
            cb = c[cols]; nz = cb > 0
            if nz.sum() == 0:
                continue
            h_b, _ = solve_em(cb[nz], Ksub[:, nz], coverage=1.0, max_iter=400, tol=1e-7,
                              normalize="per_founder", kfw=kfw_b)
            lab = labels[j]; Kb = lab.max() + 1
            csum = np.bincount(lab, weights=h_b, minlength=Kb)         # identifiable class sums
            size = np.bincount(lab, minlength=Kb)
            hclean = (csum / size)[lab]                                # equal-split back to 231
            clean_blocks.append(hclean)
        clean = np.vstack(clean_blocks)
        used_nk = nk[[j for j in range(len(block_cols)) if labels[j] is not None]][:len(clean)]
        h_uniform = clean.mean(0)
        h_nkw = (clean * used_nk[:, None]).sum(0) / used_nk.sum()
        results[sid] = dict(h_uniform=h_uniform, h_nkw=h_nkw)
        print(f"[{sid}] new global: uniform n_abs(<1e-3)={int((h_uniform<1e-3).sum())} "
              f"nkw n_abs={int((h_nkw<1e-3).sum())} ({time.time()-t0:.0f}s)", flush=True)

    np.savez_compressed(f"{OUT}/newglobal_blocks_chr1.npz",
                        founders=founders,
                        samples=np.array(SAMPLES),
                        h_uniform=np.vstack([results[s]["h_uniform"] for s in SAMPLES]),
                        h_nkw=np.vstack([results[s]["h_nkw"] for s in SAMPLES]),
                        Kbs=np.array(Kbs), nk=nk)
    print(f"saved newglobal_blocks_chr1.npz ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

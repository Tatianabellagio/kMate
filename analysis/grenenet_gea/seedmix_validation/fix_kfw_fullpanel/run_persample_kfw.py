"""Real seed-mix, before/after commit 9669be7 (full-panel Kf_w), through the REAL
production solver (src/kmate/em_solver.py) -- not a prototype.

Reuses the counts already cached by fix_seedmix/run_persample.py
(results/.../fix_seedmix/counts/{SID}_{C}.npy), counted against the production
filt2inv panel -- no jellyfish re-run needed. Current production config for
global mode: normalize="per_founder", omega=None (--kmer-weight uniform).

  OLD (pre-fix emulation): kmer_pa pre-sliced to observed (c_k>0) columns, kfw=None
      -> solve_em's fallback computes Kf_w over just that observed slice.
  NEW (the fix):           kfw computed over the FULL per-chrom panel (all k-mers,
      incl. c_k=0) BEFORE the nz slice, passed explicitly.
"""
import os, sys, time
import numpy as np
import scipy.sparse as sp

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/results/grenenet_gea/seedmix_validation/fix_kfw_fullpanel"
CDIR = f"{ROOT}/results/grenenet_gea/seedmix_validation/fix_seedmix/counts"
PANEL = f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_{{C}}"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
sys.path.insert(0, f"{ROOT}/src")
from kmate.em_solver import solve_em

SID = sys.argv[1]  # e.g. SEEDMIX_S1


def main():
    t0 = time.time()
    Hold, Hnew = [], []
    fo0 = None
    for C in CHROMS:
        tc = time.time()
        meta = np.load(PANEL.format(C=C) + ".meta.npz", allow_pickle=True)
        fo = np.asarray(meta["founders"]).astype(str)
        if fo0 is None:
            fo0 = fo
        assert np.array_equal(fo, fo0)

        c = np.load(f"{CDIR}/{SID}_{C}.npy").astype(np.float32)
        Ksp = sp.load_npz(PANEL.format(C=C) + ".kmer_pa.npz").astype(np.float32).tocsr()
        K = np.asarray(Ksp.todense(), dtype=np.float32)
        del Ksp

        kfw_full = K.sum(axis=1).astype(np.float32)  # full-panel Kf_w (production default: omega=None)
        nz = c > 0
        Knz = K[:, nz]
        cnz = c[nz]
        del K

        h_old, i_old = solve_em(cnz, Knz, coverage=1.0, max_iter=400, tol=1e-7,
                                 normalize="per_founder", kfw=None)
        h_new, i_new = solve_em(cnz, Knz, coverage=1.0, max_iter=400, tol=1e-7,
                                 normalize="per_founder", kfw=kfw_full)
        Hold.append(h_old); Hnew.append(h_new)
        print(f"[{SID} {C}] old n_abs={int((h_old<1e-3).sum())} ({i_old['iterations']}it) "
              f"new n_abs={int((h_new<1e-3).sum())} ({i_new['iterations']}it) "
              f"[{time.time()-tc:.0f}s]", flush=True)
        del Knz
        import gc; gc.collect()

    np.savez_compressed(f"{OUT}/persample_{SID}.npz", founders=fo0,
                         perchrom_old=np.vstack(Hold), perchrom_new=np.vstack(Hnew))
    print(f"[{SID}] saved persample_{SID}.npz ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

"""One seed-mix sample x 5 chromosomes, production(multinomial) vs fix(per_founder),
global mode, omega=1/mb. Saves per-sample npz with per-chrom h for both modes.
Invoked per-sample (array) to avoid concurrent-query I/O contention. Reuses prebuilt DBs.
"""
import os, sys, json, time
import numpy as np
import scipy.sparse as sp

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/qc/seedmix_validation/fix_seedmix"
CDIR = f"{OUT}/counts"; DBDIR = f"{OUT}/dbs"
PANEL = f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_{{C}}"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/analysis/grenenet_gea/qc/seedmix_validation/fix_norm")
from kmate.kmer_count import query_kmer_db
from em_variants import solve_em, make_omega

SID = sys.argv[1]                    # e.g. SEEDMIX_S1
os.makedirs(CDIR, exist_ok=True)


def counts_for(sid, C, ki):
    cache = f"{CDIR}/{sid}_{C}.npy"
    if os.path.exists(cache):
        return np.load(cache)
    cd = query_kmer_db(f"{DBDIR}/{sid}.jf", list(ki), k=31)
    c = np.array([cd[km] for km in ki], dtype=np.int64)
    np.save(cache, c)
    return c


def main():
    t0 = time.time()
    Hp, Hf = [], []; fo0 = None
    for C in CHROMS:
        tc = time.time()
        meta = np.load(PANEL.format(C=C) + ".meta.npz", allow_pickle=True)
        ki = meta["kmer_index"]; fo = np.asarray(meta["founders"]).astype(str)
        bub = np.asarray(meta["bubble_id"]).astype(np.int64)
        if fo0 is None: fo0 = fo
        assert np.array_equal(fo, fo0)
        c = counts_for(SID, C, ki).astype(np.float32)
        K = sp.load_npz(PANEL.format(C=C) + ".kmer_pa.npz").astype(np.float32).tocsr()
        om = make_omega("1/mb", None, bub)
        hp, ip = solve_em(c, K, mode="multinomial", omega=om, max_iter=400, tol=1e-7)
        hf, iff = solve_em(c, K, mode="poisson", omega=om, max_iter=400, tol=1e-7)
        Hp.append(hp); Hf.append(hf)
        print(f"[{SID} {C}] prod n_abs={int((hp<1e-3).sum())} fix n_abs={int((hf<1e-3).sum())} "
              f"[{time.time()-tc:.0f}s]", flush=True)
        del K
        import gc; gc.collect()
    np.savez_compressed(f"{OUT}/persample_{SID}.npz", founders=fo0,
                        perchrom_prod=np.vstack(Hp), perchrom_fix=np.vstack(Hf))
    print(f"[{SID}] saved persample_{SID}.npz ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

"""Genome-wide REAL seed-mix h estimates for ALL 8 replicates x 5 chromosomes,
production (multinomial) vs fix (poisson), global mode, omega=1/m_b.

Saves EVERY per-(sample,chrom,mode) h vector so we can visualize the per-chrom
scatter and check whether the chromosome-average approaches the equimolar 1/231.
Counts each sample's reads once (build DB, query each chrom's filt2inv panel;
counts cached per sample/chrom).
"""
import os, sys, json, time, gc
import numpy as np
import scipy.sparse as sp
from multiprocessing import Pool

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/seedmix_validation/fix_seedmix"
CDIR = f"{OUT}/counts"; DBDIR = os.environ.get("JF_DIR", OUT) + "/dbs"
PANEL = f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_{{C}}"
MANIFEST = f"{ROOT}/data/seedmix_manifest_arch3.tsv"
SPLIT = json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))
CACTUS = set(map(str, SPLIT["cactus"]))
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]

sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/analysis/grenenet_gea/seedmix_validation/fix_norm")
from kmate.kmer_count import build_kmer_db, query_kmer_db
from em_variants import solve_em, make_omega

# globals for workers
K_G = None; OM_G = None; KI_G = None; C_G = None


def load_manifest():
    s = []
    with open(MANIFEST) as fh:
        for i, line in enumerate(fh):
            if i == 0: continue
            sid, r1, r2 = line.rstrip("\n").split("\t")
            s.append((sid, r1, r2))
    return s


def worker(sample):
    sid, r1, r2 = sample
    cpath = f"{CDIR}/{sid}_{C_G}.npy"
    if os.path.exists(cpath):
        c = np.load(cpath)
    else:
        cd = query_kmer_db(f"{DBDIR}/{sid}.jf", list(KI_G), k=31)
        c = np.array([cd[km] for km in KI_G], dtype=np.int64)
        np.save(cpath, c)
    c = c.astype(np.float32)
    hp, _ = solve_em(c, K_G, mode="multinomial", omega=OM_G, max_iter=400, tol=1e-7)
    hf, _ = solve_em(c, K_G, mode="poisson", omega=OM_G, max_iter=400, tol=1e-7)
    return sid, hp, hf, int((hp < 1e-3).sum()), int((hf < 1e-3).sum())


def main():
    global K_G, OM_G, KI_G, C_G
    t0 = time.time()
    os.makedirs(CDIR, exist_ok=True); os.makedirs(DBDIR, exist_ok=True)
    samples = load_manifest()

    # Phase 1: build a jellyfish DB per sample (once).
    for sid, r1, r2 in samples:
        db = f"{DBDIR}/{sid}.jf"
        if os.path.exists(db):
            print(f"  DB cache hit {sid}", flush=True); continue
        t = time.time(); build_kmer_db([r1, r2], db, k=31, threads=8, hash_size="3G")
        print(f"  built DB {sid} [{time.time()-t:.0f}s]", flush=True)

    # Phase 2: per chrom, load panel once, run all samples in parallel.
    H = {sid: {} for sid, _, _ in samples}     # H[sid][chrom] = (h_prod, h_fix)
    fo0 = None
    for C in CHROMS:
        tc = time.time()
        meta = np.load(PANEL.format(C=C) + ".meta.npz", allow_pickle=True)
        KI_G = meta["kmer_index"]; fo = np.asarray(meta["founders"]).astype(str)
        bub = np.asarray(meta["bubble_id"]).astype(np.int64)
        if fo0 is None: fo0 = fo
        assert np.array_equal(fo, fo0)
        K_G = sp.load_npz(PANEL.format(C=C) + ".kmer_pa.npz").astype(np.float32).tocsr()
        ac = np.asarray(K_G.sum(axis=0)).ravel().astype(np.float32)
        OM_G = make_omega("1/mb", ac, bub); C_G = C
        print(f"[{C}] panel loaded ({K_G.shape[1]:,} kmers), running {len(samples)} samples ...", flush=True)
        with Pool(min(8, len(samples))) as p:
            res = p.map(worker, samples)
        for sid, hp, hf, nap, naf in res:
            H[sid][C] = (hp, hf)
            print(f"    {sid} {C}: n_abs prod={nap} fix={naf}", flush=True)
        del K_G; gc.collect()
        print(f"[{C}] done [{time.time()-tc:.0f}s]", flush=True)

    # Save everything + aggregates
    u = 1.0 / len(fo0); is_cac = np.array([f in CACTUS for f in fo0])
    save = {"founders": fo0, "chroms": np.array(CHROMS), "expected": np.full(len(fo0), u)}
    summ = {}
    for sid, _, _ in samples:
        Hp = np.vstack([H[sid][C][0] for C in CHROMS])   # 5 x F
        Hf = np.vstack([H[sid][C][1] for C in CHROMS])
        save[f"perchrom_prod__{sid}"] = Hp
        save[f"perchrom_fix__{sid}"] = Hf
        for tag, Hm in [("prod", Hp), ("fix", Hf)]:
            hbar = Hm.mean(0); hbar = hbar / hbar.sum()
            save[f"agg_{tag}__{sid}"] = hbar
            summ[f"{sid}_{tag}"] = dict(
                n_abs=int((hbar < 1e-3).sum()), eff_n=float(1/np.sum(hbar**2)),
                rmse_vs_uniform=float(np.sqrt(np.mean((hbar-u)**2))),
                cac_ratio=float(hbar[is_cac].sum()/is_cac.mean()))
    np.savez_compressed(f"{OUT}/seedmix_allsamples_h.npz", **save)
    json.dump(summ, open(f"{OUT}/seedmix_allsamples_summary.json", "w"), indent=2)
    print("\n=== aggregate (chrom-mean) n_absorbed per sample ===", flush=True)
    for sid, _, _ in samples:
        print(f"  {sid}: prod {summ[f'{sid}_prod']['n_abs']:3d} -> fix {summ[f'{sid}_fix']['n_abs']:3d} "
              f"| rmse_vs_1/231 {summ[f'{sid}_prod']['rmse_vs_uniform']:.2e} -> "
              f"{summ[f'{sid}_fix']['rmse_vs_uniform']:.2e}", flush=True)
    print(f"\nsaved -> seedmix_allsamples_h.npz + summary.json  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

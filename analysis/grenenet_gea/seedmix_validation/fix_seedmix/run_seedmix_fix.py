"""REAL seed-mix S1 test of the per-founder-normalization fix (Chr1).

Faithful before/after on identical input:
  - count k-mers from the REAL SEEDMIX_S1 reads against the production
    filt2inv Chr1 panel using the driver's OWN counting code (byte-identical c).
  - run the SAME EM twice, changing ONLY the M-step normalization:
      * mode="multinomial"  == PRODUCTION (h_new = em/total_c)          -> h_prod
      * mode="poisson"      == FIX        (h_new = em/Kf_w, renormalized)-> h_fix
    both with the production omega=1/m_b weighting, global mode.
Truth for the seed mix is ~uniform 1/231 (all founders equimolar by design);
hapFIRE is the external ruler (does NOT collapse founders). We compare both
kMate variants to expected 1/231 and to hapFIRE.
"""
import os, sys, json, time
import numpy as np
import scipy.sparse as sp

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/seedmix_validation/fix_seedmix"
PANEL = f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
R1 = "/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz"
R2 = "/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz"
SPLIT = json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))
CACTUS = set(map(str, SPLIT["cactus"]))
WATCH = ["9977", "9985", "10013", "9941", "9507", "9761"]
JF_DB = os.environ.get("JF_DIR", OUT) + "/seedmix_S1.jf"
COUNTS_CACHE = f"{OUT}/seedmix_S1_chr1_filt2inv_counts.npy"

sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/analysis/grenenet_gea/seedmix_validation/fix_norm")
from kmate.kmer_count import build_kmer_db, query_kmer_db
from em_variants import solve_em, make_omega


def get_counts(kmer_index):
    if os.path.exists(COUNTS_CACHE):
        print(f"  counts cache hit {COUNTS_CACHE}", flush=True)
        return np.load(COUNTS_CACHE)
    t = time.time()
    if not os.path.exists(JF_DB):
        print("  building jellyfish DB from seed-mix reads (once) ...", flush=True)
        build_kmer_db([R1, R2], JF_DB, k=31, threads=8, hash_size="3G")
        print(f"  DB built [{time.time()-t:.0f}s]", flush=True)
    print("  querying panel k-mers ...", flush=True)
    cd = query_kmer_db(JF_DB, list(kmer_index), k=31)
    counts = np.array([cd[km] for km in kmer_index], dtype=np.int64)
    np.save(COUNTS_CACHE, counts)
    print(f"  counts: nz {(counts>0).sum():,}/{len(counts):,} mean {counts.mean():.1f} "
          f"[{time.time()-t:.0f}s]", flush=True)
    return counts


def score(h, fo, hap=None):
    u = 1.0 / len(h)
    is_cac = np.array([f in CACTUS for f in fo])
    d = h - u
    m = dict(rmse_vs_uniform=float(np.sqrt(np.mean(d**2))),
             n_abs_1e3=int((h < 1e-3).sum()),
             n_lt_u10=int((h < u/10).sum()),
             eff_n=float(1.0/np.sum(h**2)),
             hmin=float(h.min()), hmax=float(h.max()),
             cac_mass=float(h[is_cac].sum()), cac_exp=float(is_cac.mean()))
    m["cac_ratio"] = m["cac_mass"]/m["cac_exp"]
    if hap is not None:
        ok = np.isfinite(hap) & np.isfinite(h)
        m["pearson_vs_hapfire"] = float(np.corrcoef(h[ok], hap[ok])[0, 1])
        m["spearman_vs_hapfire"] = float(
            np.corrcoef(np.argsort(np.argsort(h[ok])), np.argsort(np.argsort(hap[ok])))[0, 1])
    return m


def main():
    t0 = time.time()
    meta = np.load(f"{PANEL}.meta.npz", allow_pickle=True)
    kmer_index = meta["kmer_index"]
    fo = np.asarray(meta["founders"]).astype(str)
    bub = np.asarray(meta["bubble_id"]).astype(np.int64)
    F = len(fo); u = 1.0/F
    print(f"panel Chr1: F={F} K={len(kmer_index):,}", flush=True)

    counts = get_counts(kmer_index).astype(np.float32)

    print("loading sparse kmer_pa ...", flush=True)
    K = sp.load_npz(f"{PANEL}.kmer_pa.npz").astype(np.float32).tocsr()
    ac = np.asarray(K.sum(axis=0)).ravel().astype(np.float32)
    omega = make_omega("1/mb", ac, bub)  # production ω_k = 1/m_b

    # hapFIRE ruler + expected, aligned to founder order
    import csv
    hapmap, expmap = {}, {}
    with open(f"{ROOT}/analysis/grenenet_gea/seedmix_validation/seedmix_kmate_vs_hapfire_ordered.csv") as fh:
        for row in csv.DictReader(fh):
            hapmap[row["founder"]] = float(row["hapfire"])
            expmap[row["founder"]] = float(row["expected"])
    hap = np.array([hapmap.get(f, np.nan) for f in fo])

    results = {}
    hsave = {}
    for label, mode in [("prod_multinomial", "multinomial"), ("fix_poisson", "poisson")]:
        t = time.time()
        h, info = solve_em(counts, K, mode=mode, omega=omega, max_iter=400, tol=1e-7)
        s = score(h, fo, hap)
        s.update(iters=info["iters"], conv=info["conv"], secs=round(time.time()-t, 1))
        results[label] = s
        hsave[label] = h
        print(f"\n[{label}] {info['iters']}it {s['secs']}s conv={info['conv']}", flush=True)
        print(f"    rmse_vs_uniform={s['rmse_vs_uniform']:.3e}  n_abs(<1e-3)={s['n_abs_1e3']}  "
              f"eff_n={s['eff_n']:.1f}  cac_ratio={s['cac_ratio']:.2f}  "
              f"pearson_vs_hapFIRE={s.get('pearson_vs_hapfire', float('nan')):.3f}", flush=True)
        print(f"    WATCH founders (h x231, expect ~1.0):", flush=True)
        for w in WATCH:
            j = list(fo).index(w)
            print(f"      {w} {'LR' if w in CACTUS else 'SR'}: h*231={h[j]*231:8.3f}  "
                  f"(h={h[j]:.3e}, hapFIRE*231={hap[j]*231:.3f})", flush=True)

    np.savez_compressed(f"{OUT}/seedmix_fix_h.npz", founders=fo, hap=hap,
                        expected=np.array([expmap.get(f, u) for f in fo]),
                        **{f"h__{k}": v for k, v in hsave.items()})
    json.dump(results, open(f"{OUT}/seedmix_fix_summary.json", "w"), indent=2)
    print(f"\nsaved -> seedmix_fix_h.npz + seedmix_fix_summary.json  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

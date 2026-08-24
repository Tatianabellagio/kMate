"""Genome-wide REAL seed-mix S1 test of the per-founder-normalization fix.

Reproduces production's genome-wide aggregation: run per-chrom global EM on all 5
chromosomes, then average the 5 per-chrom h vectors (the 'kmate' column in
seedmix_kmate_vs_hapfire is exactly this mean; kmate_sd is its across-chrom SD).
Do it for BOTH modes, changing ONLY the M-step normalization:
    multinomial = PRODUCTION ;  poisson = FIX (divide by per-founder Kf_w).
Both use production omega=1/m_b, global mode. Count k-mers ONCE from the real
seed-mix reads (build DB, query each chrom's filt2inv panel).
Compare aggregated h to expected 1/231 and to hapFIRE (external ruler).
"""
import os, sys, json, time
import numpy as np
import scipy.sparse as sp

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/qc/seedmix_validation/fix_seedmix"
PANEL = f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_{{C}}"
R1 = "/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz"
R2 = "/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz"
SPLIT = json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))
CACTUS = set(map(str, SPLIT["cactus"]))
WATCH = ["9977", "9985", "10013", "9941", "9507", "9761"]
JF_DB = os.environ.get("JF_DIR", OUT) + "/seedmix_S1_gw.jf"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]

sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/analysis/grenenet_gea/qc/seedmix_validation/fix_norm")
from kmate.kmer_count import build_kmer_db, query_kmer_db
from em_variants import solve_em, make_omega


def counts_for(chrom, kmer_index):
    cache = f"{OUT}/seedmix_S1_{chrom}_filt2inv_counts.npy"
    if os.path.exists(cache):
        return np.load(cache)
    cd = query_kmer_db(JF_DB, list(kmer_index), k=31)
    c = np.array([cd[km] for km in kmer_index], dtype=np.int64)
    np.save(cache, c)
    return c


def main():
    t0 = time.time()
    if not os.path.exists(JF_DB):
        print("building jellyfish DB from seed-mix reads (once) ...", flush=True)
        build_kmer_db([R1, R2], JF_DB, k=31, threads=8, hash_size="3G")
        print(f"  DB built [{time.time()-t0:.0f}s]", flush=True)

    fo0 = None
    h_per = {"prod_multinomial": [], "fix_poisson": []}
    for C in CHROMS:
        tc = time.time()
        meta = np.load(PANEL.format(C=C) + ".meta.npz", allow_pickle=True)
        ki = meta["kmer_index"]; fo = np.asarray(meta["founders"]).astype(str)
        bub = np.asarray(meta["bubble_id"]).astype(np.int64)
        if fo0 is None: fo0 = fo
        assert np.array_equal(fo, fo0), f"founder order differs on {C}"
        c = counts_for(C, ki).astype(np.float32)
        K = sp.load_npz(PANEL.format(C=C) + ".kmer_pa.npz").astype(np.float32).tocsr()
        ac = np.asarray(K.sum(axis=0)).ravel().astype(np.float32)
        om = make_omega("1/mb", ac, bub)
        for label, mode in [("prod_multinomial", "multinomial"), ("fix_poisson", "poisson")]:
            h, info = solve_em(c, K, mode=mode, omega=om, max_iter=400, tol=1e-7)
            h_per[label].append(h)
            print(f"  [{C} {label}] {info['iters']}it n_abs={int((h<1e-3).sum())} "
                  f"eff_n={1/np.sum(h**2):.1f} [{time.time()-tc:.0f}s]", flush=True)
        del K
        import gc; gc.collect()

    # hapFIRE ruler + expected
    import csv
    hapmap, expmap = {}, {}
    with open(f"{ROOT}/analysis/grenenet_gea/qc/seedmix_validation/seedmix_kmate_vs_hapfire_ordered.csv") as fh:
        for row in csv.DictReader(fh):
            hapmap[row["founder"]] = float(row["hapfire"]); expmap[row["founder"]] = float(row["expected"])
    hap = np.array([hapmap.get(f, np.nan) for f in fo0])
    u = 1.0 / len(fo0); is_cac = np.array([f in CACTUS for f in fo0])

    out = {}
    hsave = {}
    for label in h_per:
        H = np.vstack(h_per[label])            # 5 x F
        h = H.mean(0); h = h / h.sum()         # genome-wide aggregate (production mean)
        sd = H.std(0)
        hsave[f"h__{label}"] = h; hsave[f"sd__{label}"] = sd
        ok = np.isfinite(hap) & np.isfinite(h)
        rk_h = np.argsort(np.argsort(h[ok])); rk_p = np.argsort(np.argsort(hap[ok]))
        m = dict(n_abs=int((h < 1e-3).sum()), eff_n=float(1/np.sum(h**2)),
                 rmse_vs_uniform=float(np.sqrt(np.mean((h-u)**2))),
                 cac_ratio=float(h[is_cac].sum()/is_cac.mean()),
                 pearson_vs_hapfire=float(np.corrcoef(h[ok], hap[ok])[0, 1]),
                 spearman_vs_hapfire=float(np.corrcoef(rk_h, rk_p)[0, 1]))
        out[label] = m
        print(f"\n[{label}] GENOME-WIDE (mean of 5 chr)", flush=True)
        print(f"    n_abs(<1e-3)={m['n_abs']}  eff_n={m['eff_n']:.1f}  "
              f"pearson_vs_hapFIRE={m['pearson_vs_hapfire']:.3f}  "
              f"spearman={m['spearman_vs_hapfire']:.3f}  cac_ratio={m['cac_ratio']:.2f}", flush=True)
        for w in WATCH:
            j = list(fo0).index(w)
            print(f"      {w} {'LR' if w in CACTUS else 'SR'}: h*231={h[j]*231:7.3f} "
                  f"(hapFIRE*231={hap[j]*231:.3f})", flush=True)

    np.savez_compressed(f"{OUT}/seedmix_gw_h.npz", founders=fo0, hap=hap,
                        expected=np.array([expmap.get(f, u) for f in fo0]), **hsave)
    json.dump(out, open(f"{OUT}/seedmix_gw_summary.json", "w"), indent=2)
    print(f"\nsaved -> seedmix_gw_h.npz + seedmix_gw_summary.json  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

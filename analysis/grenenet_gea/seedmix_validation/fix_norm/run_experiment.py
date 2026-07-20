"""Reproduce the completeness bias + test per-founder ('poisson') normalization fix.

Root cause under test: production em_solver.solve_em normalizes each M-step by a
GLOBAL scalar (total_c), so h_new_f ~= h_f * Kf_w_f / sum_f'(h_f' * Kf_w_f') at the
noiseless fixed point (Kf_w_f = founder f's omega-weighted total k-mer content).
h_true is an exact fixed point only if Kf_w is EQUAL across founders. It isn't
(PG vs cactus, filt2inv retention varies per founder) -> deterministic drift toward
high-Kf_w founders, collapsing low-Kf_w founders even with ZERO read noise. The
'poisson' M-step divides by each founder's OWN Kf_w before renormalizing, making
h_true an exact fixed point for ANY Kf_w heterogeneity (this is the effective-length
correction kallisto/RSEM/salmon apply and ALGORITHM.md's own RNA-seq analogy implies
but production is missing).

Tests, on Chr1, RAW (ac==1 kept) vs PROD (filt2inv):
  - noiseless uniform truth (identifiability floor)
  - REALISTIC NOISY uniform truth (Poisson counts at seedmix-like depth T)
  - noiseless skewed truth (does the fix preserve real non-uniform signal?)
  - noisy skewed truth
Scores: per-founder h RMSE vs truth, n absorbed (<1e-3), cactus mass ratio, and
the smoking-gun check: Spearman(h_final - h_true, Kf_w) per config -- should be
~0 for 'poisson' mode and strongly positive for 'multinomial' if the hypothesis
is right.
"""
import numpy as np, json, time, sys, gc
from multiprocessing import Pool
from scipy.sparse import load_npz
from scipy.stats import spearmanr

sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/seedmix_validation/fix_norm")
from em_variants import solve_em, make_omega

ROOT = "/global/scratch/users/tbellg/kmate"
RAW = f"{ROOT}/benchmarks/p231/data/kmer_pa_p231/kmer_pa_Chr1"
PROD = f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
SPLIT = json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))
CACTUS = set(map(str, SPLIT["cactus"]))
HARD = ["9977", "9985", "10013", "9941", "9507", "9761"]
OUT = f"{ROOT}/analysis/grenenet_gea/seedmix_validation/fix_norm"

# ---- populated in parent before fork (COW-shared with workers) ----
K_G = None
FO_G = None


def load_panel(prefix):
    K = load_npz(prefix + ".kmer_pa.npz").astype(np.float32).tocsr()
    meta = np.load(prefix + ".meta.npz", allow_pickle=True)
    fo = np.asarray(meta["founders"]).astype(str)
    bub = meta["bubble_id"] if "bubble_id" in meta.files else None
    return K, fo, bub


def score(h, fo, h_true, Kf_w):
    is_cac = np.array([f in CACTUS for f in fo])
    d = h - h_true
    rmse = float(np.sqrt(np.mean(d ** 2)))
    nabs = int((h < 1e-3).sum())
    cac_mass = float(h[is_cac].sum())
    cac_true = float(h_true[is_cac].sum())
    hard = {f: float(h[list(fo).index(f)] * 231) for f in HARD if f in fo}
    out = dict(rmse=rmse, nabs=nabs, cac_mass=cac_mass, cac_true=cac_true, hard=hard,
               cac_ratio=cac_mass / max(cac_true, 1e-12))
    if Kf_w is not None and np.std(Kf_w) > 0:
        out["spearman_resid_Kfw"] = float(spearmanr(d, Kf_w).statistic)
    if h_true.std() > 1e-9:
        out["slope"] = float(np.polyfit(h_true, h, 1)[0])
        out["pearson"] = float(np.corrcoef(h, h_true)[0, 1])
    return out


def run_one(job):
    t = time.time()
    om = make_omega(job["omk"], job["ac"], job["bub"]) if job["omk"] else None
    h, info = solve_em(job["counts"], K_G, mode=job["mode"], omega=om,
                        max_iter=job.get("max_iter", 400), tol=job.get("tol", 1e-9))
    s = score(h, FO_G, job["h_true"], info.get("Kf_w"))
    s.update(name=job["name"], iters=info["iters"], conv=info["conv"],
              secs=round(time.time() - t, 1))
    print(f"  [{job['name']:32s}] {s['iters']:3d}it {s['secs']:6.0f}s rmse={s['rmse']:.3e} "
          f"nabs={s['nabs']:3d} cac_ratio={s['cac_ratio']:.2f} "
          f"spear(resid,Kfw)={s.get('spearman_resid_Kfw', float('nan')):+.2f} "
          f"slope={s.get('slope', float('nan')):.3f}", flush=True)
    return s, h


def make_counts(h_true, K, T, seed=None):
    mu = np.asarray(h_true.astype(np.float32) @ K).ravel()
    if seed is None:
        return (T * mu).astype(np.float32)
    rng = np.random.default_rng(seed)
    return rng.poisson(T * mu).astype(np.float32)


def run_panel(name, K, fo, bub, ac):
    global K_G, FO_G
    K_G, FO_G = K, fo
    F, Kn = K.shape
    print(f"\n===== {name}  ({F}x{Kn:,}, nnz {K.nnz:,}) =====", flush=True)

    h_uni = np.full(F, 1.0 / F)
    rng = np.random.default_rng(0)
    fo_l = list(fo)
    collapsers = [fo_l.index(w) for w in HARD if w in fo_l]
    others = rng.choice([i for i in range(F) if i not in collapsers],
                         20 - len(collapsers), replace=False)
    hi_idx = np.array(sorted(set(collapsers) | set(others.tolist())))
    h_skew = np.full(F, 1.0); h_skew[hi_idx] = 6.0; h_skew /= h_skew.sum()

    T = 77  # realistic seedmix-like pool depth
    counts = {
        "uni_noiseless": make_counts(h_uni, K, T),
        "uni_noisy": make_counts(h_uni, K, T, seed=0),
        "skew_noiseless": make_counts(h_skew, K, T),
        "skew_noisy": make_counts(h_skew, K, T, seed=2),
    }
    h_true_map = {"uni_noiseless": h_uni, "uni_noisy": h_uni,
                  "skew_noiseless": h_skew, "skew_noisy": h_skew}

    configs = [("multinomial", None), ("poisson", None),
               ("multinomial", "1/ac"), ("poisson", "1/ac")]
    if bub is not None:
        configs += [("multinomial", "1/mb"), ("poisson", "1/mb")]

    jobs = []
    for cname, c in counts.items():
        for mode, omk in configs:
            jobs.append(dict(name=f"{name}_{cname}_{mode}_{omk}", counts=c,
                              h_true=h_true_map[cname], mode=mode, omk=omk,
                              ac=ac, bub=bub))

    nproc = min(8, len(jobs))
    with Pool(nproc) as p:
        res = p.map(run_one, jobs)

    summ = {}
    harr = {}
    for (s, h) in res:
        summ[s["name"]] = s
        harr[s["name"]] = h
    return summ, harr


if __name__ == "__main__":
    t0 = time.time()
    all_summ = {}
    all_h = {}
    for name, prefix in [("PROD", PROD), ("RAW", RAW)]:
        K, fo, bub = load_panel(prefix)
        ac = np.asarray(K.sum(axis=0)).ravel().astype(np.float32)
        summ, harr = run_panel(name, K, fo, bub, ac)
        all_summ.update(summ)
        for k, v in harr.items():
            all_h[f"h__{k}"] = v
        del K, fo, bub, ac
        gc.collect()

    json.dump(all_summ, open(f"{OUT}/results.json", "w"), indent=2)
    np.savez_compressed(f"{OUT}/results_h.npz", **all_h)
    print(f"\nwrote results.json + results_h.npz  ({time.time()-t0:.0f}s total)")

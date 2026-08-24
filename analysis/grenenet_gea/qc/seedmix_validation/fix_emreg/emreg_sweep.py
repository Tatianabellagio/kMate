#!/usr/bin/env python3
"""EM-regularization sweep for the kMate global-mode founder-collapse fix.

Read-free, Chr1-only, ideal-count experiment. Does a Dirichlet(alpha) prior or a
prior_h anchor stop the ~tens-of-founders collapse-to-0 WITHOUT (a) reintroducing
long-read (cactus) over-credit or (b) erasing real non-uniform / selection signal?

Key identifiability fact: with h_true=uniform, the NOISELESS ideal counts make
uniform an *exact EM fixed point*; from a uniform init the EM never moves (delta=0)
so collapse cannot appear.  The real collapse is READ NOISE pushing the estimate
off a flat likelihood ridge and sliding to the simplex boundary.  So we reproduce
collapse with Poisson counts at realistic pool depth T (expected count_k = T*mu_k),
and ALSO probe the pure ridge with a noiseless-but-perturbed-init control.

We do NOT edit src/kmate/*.  solve_em is re-implemented here for scipy-sparse K
(math identical to src/kmate/em_solver.solve_em: Dirichlet pseudo-count
alpha-1, prior_h anchor lambda*N*prior_h, omega per-kmer weights).

Runs are independent -> parallel via fork (K shared copy-on-write; workers read only).
"""
import os, sys, json, time
import numpy as np
import scipy.sparse as sp
from multiprocessing import Pool

ROOT = "/global/scratch/users/tbellg/kmate"
OUT  = f"{ROOT}/analysis/grenenet_gea/qc/seedmix_validation/fix_emreg"
PANEL = f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
SPLIT = json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))
WATCH = ["9977", "9985", "10013", "9941", "9507", "9761"]

MAXIT, TOL = 300, 1e-7

# ---- globals populated in parent before fork (shared COW) ----
Kf = None       # F x K  csr float32
Ktf = None      # K x F  csc float32 (== Kf.T, for h@K product)
FO = None        # founders (str)
IS_CAC = None    # bool F  (cactus / long-read)
OMEGA = None     # 1/m_b  K-vec float32
AC = None        # per-kmer carrier count


def solve_em_sparse(counts, h_init, alpha=1.0, prior_h=None, prior_weight=0.0,
                    omega=None, max_iter=MAXIT, tol=TOL):
    """Faithful sparse port of src/kmate/em_solver.solve_em."""
    F = Kf.shape[0]
    h = h_init.astype(np.float32).copy()
    wc = counts.astype(np.float32) if omega is None else (omega * counts).astype(np.float32)
    total_c = wc.sum()
    prior_pseudo = np.float32(max(0.0, alpha - 1.0))
    use_anchor = prior_h is not None and prior_weight > 0 and total_c > 0
    anchor = (np.float32(prior_weight) * np.float32(total_c) * prior_h.astype(np.float32)
              ) if use_anchor else None
    a_sum = float(anchor.sum()) if anchor is not None else 0.0
    delta = 1.0
    for it in range(max_iter):
        denom = np.maximum(Ktf.dot(h), np.float32(1e-7))      # = h @ K
        cw = wc / denom
        em = h * Kf.dot(cw)                                    # = K @ cw
        if anchor is not None:
            em = em + anchor
        if prior_pseudo > 0:
            h_new = (em + prior_pseudo) / (total_c + F * prior_pseudo + a_sum)
        else:
            h_new = em / max(total_c + a_sum, np.float32(1e-12))
        h_new = h_new / h_new.sum()
        delta = float(np.linalg.norm(h_new - h))
        h = h_new
        if delta < tol:
            break
    return h.astype(np.float64), it + 1, delta


def metrics(h, h_true):
    F = len(h); u = 1.0 / F
    d = h - h_true
    m = dict(
        rmse=float(np.sqrt(np.mean(d**2))),
        n_abs=int((h < 1e-3).sum()),
        n_lt_u10=int((h < u / 10).sum()),
        n_lt_u2=int((h < u / 2).sum()),
        eff_n=float(1.0 / np.sum(h**2)),
        hmin=float(h.min()), hmax=float(h.max()),
        mass_cac=float(h[IS_CAC].sum()),
        mass_pg=float(h[~IS_CAC].sum()),
        exp_cac=float(IS_CAC.mean()),
    )
    m["cac_ratio"] = m["mass_cac"] / m["exp_cac"]
    m["pg_ratio"] = m["mass_pg"] / (1 - m["exp_cac"])
    # signal fidelity vs truth (meaningful for skewed truth)
    if h_true.std() > 1e-9:
        m["pearson"] = float(np.corrcoef(h, h_true)[0, 1])
        from scipy.stats import spearmanr
        m["spearman"] = float(spearmanr(h, h_true).statistic)
        # regression slope h_est ~ h_true (1.0 = faithful, <1 = shrunk/flattened)
        A = np.vstack([h_true, np.ones_like(h_true)]).T
        m["slope"] = float(np.linalg.lstsq(A, h, rcond=None)[0][0])
    watch = {}
    fo = list(FO)
    for w in WATCH:
        if w in fo:
            j = fo.index(w)
            watch[w] = dict(h=float(h[j]), h_over_u=float(h[j] / u),
                            lr=bool(IS_CAC[j]))
    m["watch"] = watch
    return m


def run_one(job):
    t = time.time()
    h, iters, delta = solve_em_sparse(
        job["counts"], job["h_init"], alpha=job.get("alpha", 1.0),
        prior_h=job.get("prior_h"), prior_weight=job.get("prior_weight", 0.0),
        omega=job.get("omega"))
    m = metrics(h, job["h_true"])
    m.update(name=job["name"], iters=iters, delta=delta,
             secs=round(time.time() - t, 1), h=h)
    print(f"  [{job['name']}] {iters}it {m['secs']}s  eff_n={m['eff_n']:.1f} "
          f"n_abs={m['n_abs']} rmse={m['rmse']:.2e} cac={m['cac_ratio']:.2f} "
          f"slope={m.get('slope', float('nan')):.3f}", flush=True)
    return m


def make_counts(h_true, T, seed=0, noiseless=False):
    mu = np.asarray(Ktf.dot(h_true.astype(np.float32))).ravel()   # = K.T @ h_true
    if noiseless:
        return (T * mu).astype(np.float32)
    rng = np.random.default_rng(seed)
    return rng.poisson(T * mu).astype(np.float32)


def main():
    global Kf, Ktf, FO, IS_CAC, OMEGA, AC
    t0 = time.time()
    K = sp.load_npz(f"{PANEL}.kmer_pa.npz").astype(np.float32)     # F x K csr
    Kf = K
    Ktf = K.T.tocsc()
    meta = np.load(f"{PANEL}.meta.npz", allow_pickle=True)
    FO = meta["founders"].astype(str)
    bid = np.asarray(meta["bubble_id"]).astype(np.int64)
    m_b = np.bincount(bid)[bid].astype(np.float32)
    OMEGA = (1.0 / m_b).astype(np.float32)
    CAC = set(map(str, SPLIT["cactus"]))
    IS_CAC = np.array([f in CAC for f in FO])
    AC = np.asarray(K.sum(axis=0)).ravel()
    F = K.shape[0]; u = 1.0 / F
    print(f"loaded K={K.shape} nnz={K.nnz:,} cactus={IS_CAC.sum()} PG={(~IS_CAC).sum()} "
          f"[{time.time()-t0:.0f}s]", flush=True)

    h_uni = np.full(F, u)
    # skewed "selection" truth: 20 founders elevated. Pick 10 that collapse + 10 random
    # so the test is fair (not just protecting the collapsers).
    rng = np.random.default_rng(1)
    fo = list(FO)
    collapsers = [fo.index(w) for w in WATCH if w in fo]
    others = rng.choice([i for i in range(F) if i not in collapsers],
                        20 - len(collapsers), replace=False)
    hi_idx = np.array(sorted(set(collapsers) | set(others.tolist())))
    # strong skew: hi founders ~ 6x low
    h_skew = np.full(F, 1.0); h_skew[hi_idx] = 6.0; h_skew /= h_skew.sum()
    # mild skew (subtle selection, ~1.5x) to test detection of small effects
    h_mild = np.full(F, 1.0); h_mild[hi_idx] = 1.5; h_mild /= h_mild.sum()

    T = 77          # realistic pool depth (expected mean count ~= real seedmix)
    c_uni = make_counts(h_uni, T, seed=0)
    c_uni_noiseless = make_counts(h_uni, T, noiseless=True)
    c_skew = make_counts(h_skew, T, seed=2)
    c_mild = make_counts(h_mild, T, seed=3)
    print(f"count scales: uni mean={c_uni.mean():.1f} nz={100*(c_uni>0).mean():.0f}% | "
          f"skew mean={c_skew.mean():.1f}", flush=True)

    jobs = []
    def add(**kw): jobs.append(kw)

    # ---- G1: alpha sweep, uniform truth, omega=1/m_b (PRODUCTION weighting) ----
    for a in [1.0, 1.001, 1.01, 1.05, 1.1, 1.5]:
        add(name=f"uni_om_a{a}", counts=c_uni, h_init=h_uni, h_true=h_uni,
            alpha=a, omega=OMEGA)
    # ---- G2: alpha sweep, uniform truth, omega=None (MLE) ----
    for a in [1.0, 1.05, 1.5]:
        add(name=f"uni_none_a{a}", counts=c_uni, h_init=h_uni, h_true=h_uni,
            alpha=a, omega=None)
    # ---- G3: STRONG skew truth, omega=1/m_b (does alpha flatten real signal?) ----
    for a in [1.0, 1.01, 1.05, 1.1, 1.5]:
        add(name=f"skew_om_a{a}", counts=c_skew, h_init=h_uni, h_true=h_skew,
            alpha=a, omega=OMEGA)
    # ---- G3b: MILD skew (subtle selection) ----
    for a in [1.0, 1.05, 1.1, 1.5]:
        add(name=f"mild_om_a{a}", counts=c_mild, h_init=h_uni, h_true=h_mild,
            alpha=a, omega=OMEGA)
    # ---- G4: prior_h=uniform anchor sweep (uniform truth) ----
    for w in [0.05, 0.2, 0.5]:
        add(name=f"uni_om_pw{w}", counts=c_uni, h_init=h_uni, h_true=h_uni,
            alpha=1.0, prior_h=h_uni, prior_weight=w, omega=OMEGA)
    # ---- G5: EVOLVED — anchor to p0 (=h_uni, the seed-mix start) while truth is
    #         skewed by selection. Does p0-anchor regularize WITHOUT biasing the
    #         estimated selection Δp=h-p0?  Compare vs uniform-dirichlet (G3). ----
    for w in [0.2, 0.5]:
        add(name=f"skew_p0anchor_pw{w}", counts=c_skew, h_init=h_uni,
            h_true=h_skew, prior_h=h_uni, prior_weight=w, omega=OMEGA)
    # ---- G6: noiseless ridge control (perturbed init) ----
    r_init = (h_uni * rng.uniform(0.5, 1.5, F)); r_init /= r_init.sum()
    r_init = r_init.astype(np.float32)
    for a in [1.0, 1.05]:
        add(name=f"noiseless_perturb_a{a}", counts=c_uni_noiseless, h_init=r_init,
            h_true=h_uni, alpha=a, omega=OMEGA)
    add(name="noiseless_uniinit_a1.0", counts=c_uni_noiseless, h_init=h_uni,
        h_true=h_uni, alpha=1.0, omega=OMEGA)

    print(f"running {len(jobs)} EM jobs ...", flush=True)
    nproc = min(8, len(jobs))
    with Pool(nproc) as p:
        res = p.map(run_one, jobs)

    # save
    save = {"founders": FO, "is_cactus": IS_CAC, "hi_idx": hi_idx,
            "h_skew_true": h_skew, "h_mild_true": h_mild, "T": T, "ac": AC}
    summ = {}
    for m in res:
        save[f"h__{m['name']}"] = m.pop("h")
        summ[m["name"]] = m
    np.savez_compressed(f"{OUT}/emreg_sweep.npz",
                        summaries=json.dumps(summ), **{k: v for k, v in save.items()})
    with open(f"{OUT}/emreg_summary.json", "w") as fh:
        json.dump(summ, fh, indent=1)
    print(f"\nsaved -> {OUT}/emreg_sweep.npz  ({time.time()-t0:.0f}s total)")


if __name__ == "__main__":
    main()

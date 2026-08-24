#!/usr/bin/env python
"""Beta-binomial (overdispersion-aware) latent-factor climate-GEA.

WHY: the plain plot-level binomial GLM  [alt,ref] ~ const + z(bio1) + LF1..LFK
treats AF*flowers*2 (up to ~1350 genomes/pool) as that many INDEPENDENT Bernoulli
draws, so the sampling variance of the observed pool AF is ~100x too small ->
GIF ~= 105 raw, and even K=16 latent factors (which fix STRUCTURE, not the
LIKELIHOOD) leave it hugely inflated (min p ~ 1e-172). The correct model for
pool-seq allele counts is a BETA-BINOMIAL with extra-binomial variance rho.

MODEL (per variant, across the 355 pools):
    counts_j ~ BetaBinomial(N_j, mu_j, rho),   logit(mu_j)=const+beta*z(bio1)+LF1..LFK
where N_j = round(total_flowers_j * 2). Test beta (climate slope).

TRACTABILITY (per-variant full BB-MLE over ~2.7M variants is infeasible):
  We use the FIXED-rho beta-binomial, which for the GLM is EXACTLY an effective-
  sample-size binomial. For a beta-binomial the variance of the observed
  proportion p_hat=s/N is
        Var(p_hat) = mu(1-mu)/N * [1 + (N-1)*rho]         (design effect DE)
  A binomial GLM run on EFFECTIVE trials
        N_eff_j = N_j / (1 + (N_j-1)*rho)
  reproduces exactly this Var(p_hat) (mu(1-mu)/N_eff = mu(1-mu)*DE/N), hence the
  correct Fisher information and Wald test for beta -- at plain-binomial-GLM cost,
  reusing the existing _chunk machinery verbatim. As N->inf, N_eff -> 1/rho (the
  pool-seq effective-sample-size cap). This is the standard Williams (1982)
  quasi-likelihood / design-effect correction for overdispersed binomials.

rho ESTIMATION (--estimate-rho mode): Williams' POOLED MOMENT estimator from the
FULL (K=16) design residuals on a random subsample of variants:
    fit binomial GLM per variant, Pearson X2_v = sum_j (s-N mu)^2/(N mu(1-mu));
    E[X2_v] ~= (n_v - p) + rho * sum_j (N_j - 1)
    rho_hat = sum_v (X2_v - (n_v-p)) / sum_v sum_j (N_j - 1)
computed GLOBAL and PER-CHROMOSOME; written to rho_{cls}.json. Estimating rho from
the K=16 residuals strips structure-driven inflation, leaving pure overdispersion;
the SAME rho is then applied to both the K=16 and K=0 genome passes so the
K0-vs-K16 contrast isolates STRUCTURE with overdispersion held fixed (2x2 decomp).

Usage (basic env, statsmodels):
  $PY run_betabinom_latent.py --class nonsnp --estimate-rho --nsub 15000 --threads 16
  $PY run_betabinom_latent.py --class nonsnp --k 16 --rho-json rho_nonsnp.json --threads 16
  $PY run_betabinom_latent.py --class nonsnp --k 0  --rho-json rho_nonsnp.json --threads 16
"""
from __future__ import annotations
import argparse, os, sys, json, warnings
import numpy as np
import pandas as pd
from multiprocessing import Pool
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

warnings.filterwarnings("ignore")

_X = None       # [pools x (2+K)] design: const + z(clim) + K latent factors
_SUCC = None    # [rec x pools] (effective) alt counts
_FAIL = None    # [rec x pools] (effective) ref counts
_N = None       # [pools] raw trials (for rho moment estimator)
_CHROM = None   # [rec] chrom codes (for per-chrom rho accumulation)


# --------------------------------------------------------------------------- #
#  latent factors  (identical construction to run_binomial_latent.py)
# --------------------------------------------------------------------------- #
def latent_factors(af: np.ndarray, k: int) -> np.ndarray:
    if k == 0:
        return np.empty((af.shape[0], 0), dtype=np.float64)
    Xc = af.astype(np.float64)
    Xc -= Xc.mean(axis=0, keepdims=True)
    G = Xc @ Xc.T
    w, V = np.linalg.eigh(G)
    U = V[:, ::-1][:, :k]
    return StandardScaler().fit_transform(U)


# --------------------------------------------------------------------------- #
#  genome pass: binomial GLM on EFFECTIVE counts  (== fixed-rho beta-binomial)
# --------------------------------------------------------------------------- #
def _chunk(rng):
    a, b = rng
    slope = np.full(b - a, np.nan); pv = np.full(b - a, np.nan)
    for i in range(a, b):
        succ = _SUCC[i]; fail = _FAIL[i]; n = succ + fail
        ok = n > 0
        if ok.sum() < _X.shape[1] + 1:
            continue
        s = succ[ok]; f = fail[ok]
        if s.sum() == 0 or f.sum() == 0:
            continue
        X = _X[ok]
        if np.ptp(X[:, 1]) == 0:
            continue
        try:
            res = sm.GLM(np.column_stack([s, f]), X,
                         family=sm.families.Binomial()).fit()
            slope[i - a] = res.params[1]; pv[i - a] = res.pvalues[1]
        except Exception:
            pass
    return a, slope, pv


# --------------------------------------------------------------------------- #
#  rho pass: Williams pooled moment estimator (numerator, denom per variant)
# --------------------------------------------------------------------------- #
def _rho_chunk(rng):
    """Return list of (chrom_code, rho_hat_variant) Williams moment estimates.

    Per variant rho_hat = (X2 - (n_ok - p)) / sum_j (N_j - 1). The pooled
    sum-of-numerators / sum-of-denominators is NOT robust: a few rare variants
    with a near-zero fitted-variance pool give astronomically large Pearson X2
    and dominate. We return per-variant estimates and aggregate by MEDIAN, which
    is robust. To keep the moment estimator well conditioned we (a) only fit
    common variants (MAF gate applied before dispatch) and (b) drop pools whose
    fitted mu is within EPS of 0/1 (var ~ 0 -> unstable Pearson term)."""
    a, b = rng
    p = _X.shape[1]
    EPS = 1e-3
    out = []
    for i in range(a, b):
        s = _SUCC[i]; f = _FAIL[i]; N = s + f          # raw counts here
        ok = N > 0
        if ok.sum() < p + 5:
            continue
        so = s[ok]; fo = f[ok]; No = N[ok]
        if so.sum() == 0 or fo.sum() == 0:
            continue
        X = _X[ok]
        if np.ptp(X[:, 1]) == 0:
            continue
        try:
            res = sm.GLM(np.column_stack([so, fo]), X,
                         family=sm.families.Binomial()).fit()
            mu = res.fittedvalues                       # fitted proportions
            good = (mu > EPS) & (mu < 1 - EPS)          # drop separated pools
            if good.sum() < p + 5:
                continue
            var = No[good] * mu[good] * (1.0 - mu[good])
            X2 = float((((so[good] - No[good] * mu[good]) ** 2) / var).sum())
            v_num = X2 - (int(good.sum()) - p)
            v_den = float((No[good] - 1.0).sum())
            if v_den <= 0:
                continue
            out.append((int(_CHROM[i]), v_num / v_den))
        except Exception:
            continue
    return out


# --------------------------------------------------------------------------- #
def load(cls, gen, climate, cmdir, k):
    recs = pd.read_csv(f"{cmdir}/{cls}_gen{gen}.records.csv")
    pools = pd.read_csv(f"{cmdir}/gen{gen}.pools.csv")
    clim = pools[climate].to_numpy(float)
    N = np.round(pools["total_flowers"].to_numpy(float) * 2).astype(np.float64)
    af = np.load(f"{cmdir}/{cls}_gen{gen}_af.npy")
    assert af.shape[0] == len(clim) and af.shape[1] == len(recs)
    U = latent_factors(af, k)
    z = StandardScaler().fit_transform(clim.reshape(-1, 1))
    X = np.column_stack([np.ones(len(clim)), z.ravel(), U])
    aft = np.ascontiguousarray(af.T); del af
    return recs, N, X, aft


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True,
                    choices=["snp", "sv", "smallindel", "nonsnp"])
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--cmdir", default=f"{lib.GEA}/phase1_replication/results/class_matrices")
    ap.add_argument("--out", default=f"{lib.GEA}/gea_newpanel/results/betabinom")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--estimate-rho", action="store_true",
                    help="rho-estimation mode: fit K=16 subsample, write rho_{cls}.json")
    ap.add_argument("--nsub", type=int, default=15000, help="variants sampled for rho")
    ap.add_argument("--rho-maf", type=float, default=0.05,
                    help="MAF floor for the (well-conditioned) rho subsample")
    ap.add_argument("--rho-json", default=None, help="per-chrom rho JSON for genome pass")
    ap.add_argument("--rho", type=float, default=None, help="override: single global rho")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    global _X, _SUCC, _FAIL, _N, _CHROM

    # -------------------------------------------------- rho estimation mode
    if args.estimate_rho:
        recs, N, X, aft = load(args.cls, args.gen, args.climate, args.cmdir, 16)
        chroms = sorted(recs["chrom"].unique())
        cmap = {c: i for i, c in enumerate(chroms)}
        chrom_code = recs["chrom"].map(cmap).to_numpy()
        rng = np.random.default_rng(args.seed)
        n_rec = aft.shape[0]                      # aft is [rec x pools]
        # MAF gate: estimator is only well conditioned on common variants
        maf = recs["maf"].to_numpy(float)
        pool_of = np.flatnonzero(maf >= args.rho_maf)
        idx = rng.choice(pool_of, size=min(args.nsub, len(pool_of)), replace=False)
        idx.sort()
        finite = np.isfinite(aft[idx])
        succ = np.where(finite, aft[idx] * N[None, :], 0.0)   # RAW counts
        fail = np.where(finite, (1.0 - aft[idx]) * N[None, :], 0.0)
        _X = X; _SUCC = succ; _FAIL = fail; _N = N
        _CHROM = chrom_code[idx]
        print(f"[rho] {args.cls}: fitting {len(idx):,} subsampled variants "
              f"(MAF>={args.rho_maf}, K=16, {X.shape[1]} cols) for Williams moment rho",
              flush=True)
        step = 500
        ranges = [(a, min(a + step, len(idx))) for a in range(0, len(idx), step)]
        vals = []
        with Pool(args.threads) as pool:
            for res in pool.imap_unordered(_rho_chunk, ranges):
                vals.extend(res)
        vals = np.array(vals) if vals else np.zeros((0, 2))
        rho_g = float(max(np.median(vals[:, 1]), 0.0)) if len(vals) else 0.0
        per = {}
        for c in chroms:
            m = vals[vals[:, 0] == cmap[c], 1] if len(vals) else np.array([])
            per[c] = float(max(np.median(m), 0.0)) if len(m) else rho_g
        out = {"class": args.cls, "global": rho_g, "per_chrom": per,
               "nsub": int(len(idx)), "n_used": int(len(vals)),
               "design": "K=16", "n_params": int(X.shape[1]),
               "estimator": "median per-variant Williams moment",
               "rho_maf": args.rho_maf}
        jp = f"{args.out}/rho_{args.cls}.json"
        json.dump(out, open(jp, "w"), indent=2)
        print(f"[rho] {args.cls} GLOBAL rho={rho_g:.5f} | per-chrom "
              + " ".join(f"{c}={per[c]:.4f}" for c in chroms) + f"\n  -> {jp}", flush=True)
        # report an effective-sample-size cap for intuition
        Nmed = float(np.median(N))
        print(f"[rho] median N={Nmed:.0f} -> DE={1+(Nmed-1)*rho_g:.1f}, "
              f"N_eff_cap=1/rho={1/rho_g if rho_g>0 else np.inf:.0f}", flush=True)
        return

    # -------------------------------------------------- genome pass
    recs, N, X, aft = load(args.cls, args.gen, args.climate, args.cmdir, args.k)
    n_rec = aft.shape[0]                          # aft is [rec x pools]

    # per-record rho -> effective trials
    if args.rho is not None:
        rho_rec = np.full(n_rec, float(args.rho))
        rho_desc = f"global {args.rho:.5f}"
    elif args.rho_json is not None:
        jp = args.rho_json
        if not os.path.isabs(jp):
            jp = f"{args.out}/{jp}"
        rj = json.load(open(jp))
        per = rj["per_chrom"]
        rho_rec = recs["chrom"].map(per).to_numpy(float)
        rho_rec = np.where(np.isfinite(rho_rec), rho_rec, rj["global"])
        rho_desc = f"per-chrom from {os.path.basename(jp)} (global {rj['global']:.5f})"
    else:
        raise SystemExit("genome pass needs --rho or --rho-json")

    finite = np.isfinite(aft)
    # N_eff_ij = N_j / (1 + (N_j-1)*rho_i)   [rec x pools]
    Neff = N[None, :] / (1.0 + (N[None, :] - 1.0) * rho_rec[:, None])
    succ = np.where(finite, aft * Neff, 0.0)
    fail = np.where(finite, (1.0 - aft) * Neff, 0.0)
    del Neff, aft, finite

    _X = X; _SUCC = succ; _FAIL = fail

    print(f"betabinom {args.cls} gen{args.gen}: {X.shape[0]} pools x {n_rec:,} "
          f"records vs {args.climate}, K={args.k} ({X.shape[1]} cols), rho={rho_desc}",
          flush=True)

    step = 5000
    ranges = [(a, min(a + step, n_rec)) for a in range(0, n_rec, step)]
    slope = np.full(n_rec, np.nan); pv = np.full(n_rec, np.nan)
    with Pool(args.threads) as pool:
        for a, sl, p in pool.imap_unordered(_chunk, ranges):
            slope[a:a + len(sl)] = sl; pv[a:a + len(p)] = p

    recs = recs.rename(columns={"maf": "MAF"}).assign(
        slope=slope, pval=pv, rho=rho_rec)
    out = f"{args.out}/betabinom_lf{args.k}_{args.cls}_gen{args.gen}_{args.climate}.csv"
    recs[["chrom", "pos", "ref_len", "alt_len", "MAF", "block",
          "slope", "pval", "rho"]].to_csv(out, index=False)
    fin = np.isfinite(pv)
    from scipy.stats import chi2
    lam = np.median(chi2.isf(pv[fin], 1)) / chi2.ppf(0.5, 1)
    print(f"  tested {int(fin.sum()):,}/{n_rec:,} | GIF={lam:.2f} | "
          f"p<1e-5: {int((pv[fin]<1e-5).sum()):,} | p<0.05: {int((pv[fin]<0.05).sum()):,} "
          f"| min p={np.nanmin(pv):.2e}\n  -> {out}", flush=True)


if __name__ == "__main__":
    main()

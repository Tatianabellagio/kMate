#!/usr/bin/env python3
"""Calibrate the AF identifiability floor c for the total AF SE.

From af_certainty.py: AF error is bias-dominated (a coverage-independent ~0.01
identifiability floor), so the pure-variance Fisher delta-method SE under-covers
(95% interval covers only 13-32%). The fix chosen by TB: report a CALIBRATED total
SE  SE_total = sqrt(SE_Fisher^2 + c^2)  where c is a benchmarked panel-level floor.

This finds c per panel: the value making the empirical 95% interval cover the TRUE
AF in 95% of records, on a closed-loop g0 pool (truth = pool weights; no recomb so
global h_true is exact). Reports c at cov 10 and 30 (should be ~stable since c
captures the coverage-independent bias) for `defined` and `well-called` records.

Usage:  af_calibrate_floor.py --panel {p80,231}   (run via sbatch)

2026-07-06: GLOBAL mode now normalizes per_founder (the Kf_w fix) and drops
omega=1/m_b (PIPELINE_STATE.md Sec.0) -- omega is loaded below for reference
but NOT passed to solve_em/fisher_cov_h, matching production's --kmer-weight
uniform default. Re-run after any EM/weighting change: this calibrates the
AF_ID_FLOOR_DEFAULT constant baked into per_sample_per_chrom.py.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse

ROOT = Path("/global/scratch/users/tbellg/kmate")
sys.path.insert(0, str(ROOT / "src"))
from kmate.em_solver import solve_em                              # noqa
from kmate.h_uncertainty import fisher_cov_h, af_se_from_cov      # noqa

PANELS = {
    "p80": dict(
        kpre=ROOT / "benchmarks/p80/data/kmer_pa_p80_filt2/kmer_pa_Chr1",
        vpre=ROOT / "benchmarks/p80/data/var_pa_p80",
        pool=ROOT / "benchmarks/p80/sims/cov10_n50_g0_s42_hotspots_p80_chr1/pool_weights.tsv"),
    "231": dict(
        kpre=ROOT / "data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1",
        vpre=ROOT / "panel/arch3/chr1/var_pa_231_arch3_chr1",
        pool=ROOT / "benchmarks/p231/sims/cov10_n50_g0_s42_hotspots_p231_chr1/pool_weights.tsv"),
}
COVS = [10.0, 30.0]
REPS = 2
SUBSAMPLE = 200_000
SEED = 20260630


def load(panel):
    p = PANELS[panel]
    kp = sparse.load_npz(f"{p['kpre']}.kmer_pa.npz")
    K = np.asarray(kp.todense(), np.float32) if sparse.issparse(kp) else np.asarray(kp, np.float32)
    m = np.load(f"{p['kpre']}.meta.npz", allow_pickle=True)
    fk = np.asarray(m["founders"]).astype(str); bid = np.asarray(m["bubble_id"])
    omega = (1.0 / np.bincount(bid)[bid]).astype(np.float32)
    V = np.asarray(sparse.load_npz(f"{p['vpre']}.var_pa.npz").todense(), np.float64)
    U = np.asarray(sparse.load_npz(f"{p['vpre']}.var_called.npz").todense(), np.float64)
    fv = np.asarray(np.load(f"{p['vpre']}.meta.npz", allow_pickle=True)["founders"]).astype(str)
    pos = {f: i for i, f in enumerate(fv)}; ridx = np.array([pos[f] for f in fk])
    V, U = V[ridx], U[ridx]
    w = pd.read_csv(p["pool"], sep="\t"); wm = {str(f): float(x) for f, x in zip(w["founder"], w["weight"])}
    h = np.array([wm.get(str(f), 0.0) for f in fk], np.float64); h /= h.sum()
    return K, omega, V, U, h


def solve_c_for_coverage(err, se, target=0.95, lo=0.0, hi=0.2):
    """Binary-search c so mean(|err| <= 1.96*sqrt(se^2+c^2)) == target."""
    def cov(c):
        return float((err <= 1.96 * np.sqrt(se * se + c * c)).mean())
    if cov(hi) < target:
        return hi, cov(hi)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if cov(mid) < target:
            lo = mid
        else:
            hi = mid
    return hi, cov(hi)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--panel", choices=list(PANELS), required=True)
    args = ap.parse_args(); t0 = time.time()
    K, omega, V, U, h_true = load(args.panel)
    F, Kn = K.shape; R = V.shape[1]
    mu = (h_true @ K).astype(np.float64)
    af_true = (h_true @ V) / np.maximum(h_true @ U, 1e-12)
    defined = (h_true @ U) > 0
    callrate = (U[h_true > 0] > 0).sum(0) / (h_true > 0).sum()
    well = defined & (callrate >= 0.99)
    chunk = 1_000_000 if F > 120 else None
    print(f"panel={args.panel} F={F} K={Kn:,} R={R:,} defined={defined.sum():,} "
          f"well={well.sum():,} chunk={chunk}  load {time.time()-t0:.0f}s", flush=True)

    rng = np.random.default_rng(SEED)
    # fixed record subsample (defined) for the SE delta-method cost
    didx = np.flatnonzero(defined)
    sub = didx if didx.size <= SUBSAMPLE else rng.choice(didx, SUBSAMPLE, replace=False)
    well_sub = well[sub]
    Vs, Us, aft_s = V[:, sub], U[:, sub], af_true[sub]

    rows = []
    for lam in COVS:
        cs_def, cs_well, fish95, mae = [], [], [], []
        for r in range(REPS):
            c = rng.poisson(lam * mu)
            nz = c > 0
            # GLOBAL mode: uniform kmer-weight (omega=None), per_founder normalize
            # (solve_em default) -- matches the production recipe.
            h, _ = solve_em(c[nz].astype(np.float32), K[:, nz], lam,
                            max_iter=1500, tol=1e-9, omega=None)
            Sig, sp = fisher_cov_h(h, K, c.astype(np.float64), omega=None, chunk=chunk)
            af_hat = (h @ Vs) / np.maximum(h @ Us, 1e-12)
            _, se = af_se_from_cov(h, Sig, Vs, Us, support=sp)
            err = np.abs(af_hat - aft_s)
            fish95.append(float((err <= 1.96 * se).mean()))
            mae.append(float(err.mean()))
            c_def, _ = solve_c_for_coverage(err, se)
            c_well, _ = solve_c_for_coverage(err[well_sub], se[well_sub])
            cs_def.append(c_def); cs_well.append(c_well)
        rows.append(dict(cov=lam, af_mae=np.mean(mae), fisher95=np.mean(fish95),
                         c_defined=np.mean(cs_def), c_wellcalled=np.mean(cs_well)))
        print(f"  cov {lam:>4g}: AF MAE={np.mean(mae):.5f}  Fisher95={np.mean(fish95):.3f}  "
              f"-> c(defined)={np.mean(cs_def):.4f}  c(wellcalled)={np.mean(cs_well):.4f}", flush=True)

    df = pd.DataFrame(rows)
    out = ROOT / f"benchmarks/h_uncertainty/results/af_floor_{args.panel}.tsv"
    df.to_csv(out, sep="\t", index=False)
    c_final = float(df["c_wellcalled"].max())   # conservative: cover the higher-coverage (tighter-SE) case
    print(f"\n=== CALIBRATED FLOOR (panel {args.panel}) ===", flush=True)
    print(f"  recommended c = {c_final:.4f}  (max c_wellcalled across cov; conservative)", flush=True)
    print(f"  -> SE_total = sqrt(SE_Fisher^2 + {c_final:.4f}^2);  wrote {out}  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()

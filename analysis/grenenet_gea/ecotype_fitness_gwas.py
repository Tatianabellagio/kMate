"""Tier 2 — 231-ecotype kinship-corrected fitness-GWAS (the SV-selection test at the
identifiable unit).

Phenotype = per-founder realized fitness (ecotype_fitness.csv, 8 axes: {relative,census} x
{global, cold, mid, hot} zone means). Genotype = the founder carrier matrix var_pa (231 x
variants, all classes). Kinship from genome-wide SNPs absorbs clade-sorting (the confound that
made every prior per-variant test look significant).

Model: EMMAX-style LMM with LOCO kinship. Per chromosome c the GRM is built from the OTHER four
chromosomes' SNPs (leave-one-chromosome-out) so a variant never deflates its own signal;
K_c = Q diag(d) Q', estimate delta per phenotype (REML, 1-D), vectorized GLS Wald z per variant.
We ALSO emit a NAIVE (no-kinship, ordinary OLS) z for the SAME variants/phenotypes — the user
wants the contrast: naive is inflated by clade structure (ecotype sorting), LOCO tames it, and
the gap quantifies how much apparent SV/SNP signal is just sorting.

Two streaming passes over the 5 per-chrom var_pa npz:
  pass 1 -> accumulate per-chrom SNP GRM contribution (SNPs, MAF>=0.05, call>=0.9, standardized)
  pass 2 -> per chrom: LOCO-K eigendecomp + per-variant z (LMM) and z_naive (OLS)

Variant class from panel meta ref_len/alt_len: SNP (1,1) | indel (non-SNP, |dlen|<=50) |
SV (|dlen|>50). Uncalled founders imputed to the column allele frequency.

Output -> results/grenenet_gea/ecotype_fitness/gwas/
  gwas_z.npz : chrom,pos,ref_len,alt_len,mac,vclass, Z[n_var x P] (LOCO), Znaive[n_var x P],
               pheno_names, delta[C x P], h2[C x P]
Env: kmate. Heavy -> run via run_ecotype_gwas.sbatch (NOT the login node).
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize_scalar
from scipy import stats
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

OUT = f"{lib.GEA}/ecotype_fitness/gwas"
CHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5"]
PHENOS = [f"{a}_{f}" for f in ("rel", "cen")
          for a in ("w_global", "w_cold", "w_mid", "w_hot")]
MIN_MAC = 12          # founder MAF>=~5% floor (COMMON_MAC) for TESTED variants
CALL_MIN = 0.9        # call-rate floor
BLK = 200_000         # variant block size for streaming
N = 231


def qn(x):
    """Rank inverse-normal transform of one phenotype (NaN-safe). Tames single-founder
    outliers (e.g. 9748, near-absent-at-founding -> top cold-winner, which LOCO-inflated the
    raw w_cold axis to lambda~2) while preserving the winner RANKING, and puts the relative
    and census flavours on a common scale."""
    x = np.asarray(x, float); out = np.full_like(x, np.nan); m = np.isfinite(x)
    r = stats.rankdata(x[m]); out[m] = stats.norm.ppf((r - 0.5) / m.sum())
    return out


def _load_chrom(cl):
    """Return (pos, ref_len, alt_len, vp csc founders x var, vc csc, n_alt, n_cal)."""
    base = f"{lib.PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
    meta = np.load(f"{base}.meta.npz", allow_pickle=True)
    vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
    vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
    n_alt = np.asarray(vp.sum(0)).ravel().astype(float)
    n_cal = np.asarray(vc.sum(0)).ravel().astype(float)
    return (meta["pos"].astype(np.int64), meta["ref_len"].astype(int),
            meta["alt_len"].astype(int), vp, vc, n_alt, n_cal)


def _dense_geno(vp, vc, n_alt, n_cal, cols):
    """Dense genotype block 231 x len(cols); uncalled cells -> column allele freq. + freq p."""
    g = vp[:, cols].toarray().astype(np.float64)
    called = vc[:, cols].toarray().astype(bool)
    p = np.where(n_cal[cols] > 0, n_alt[cols] / np.maximum(n_cal[cols], 1), 0.0)
    miss = ~called
    if miss.any():
        g[miss] = np.take(p, np.where(miss)[1])
    return g, p


def _snp_mask(rl, al, n_alt, n_cal):
    return (rl == 1) & (al == 1) & (n_alt >= MIN_MAC) & (n_alt <= N - MIN_MAC) \
        & (n_cal >= CALL_MIN * N)


def kinship_contrib(cl):
    """Unnormalized SNP GRM contribution (sum Z Z') and SNP count for one chromosome."""
    pos, rl, al, vp, vc, n_alt, n_cal = _load_chrom(cl)
    idx = np.where(_snp_mask(rl, al, n_alt, n_cal))[0]
    Kc = np.zeros((N, N))
    for s in range(0, len(idx), BLK):
        cols = idx[s:s + BLK]
        g, p = _dense_geno(vp, vc, n_alt, n_cal, cols)
        sd = np.sqrt(np.maximum(p * (1 - p), 1e-6))
        Z = (g - p) / sd
        Kc += Z @ Z.T
    return Kc, len(idx)


def _reml_delta(yt, Xt, d):
    """REML delta=sigma_e^2/sigma_g^2 in the rotated null LMM (1-D, log-scale). Returns w,s2."""
    n = len(yt); q = Xt.shape[1]

    def negll(logdelta):
        delta = np.exp(logdelta)
        w = 1.0 / (d + delta)
        XtW = Xt * w[:, None]
        A = Xt.T @ XtW
        b = np.linalg.solve(A, XtW.T @ yt)
        r = yt - Xt @ b
        s2 = float((w * r * r).sum()) / (n - q)
        return 0.5 * ((n - q) * np.log(s2) + np.log(d + delta).sum()
                      + np.linalg.slogdet(A)[1])

    res = minimize_scalar(negll, bounds=(-8, 8), method="bounded")
    delta = float(np.exp(res.x))
    w = 1.0 / (d + delta)
    XtW = Xt * w[:, None]
    b = np.linalg.solve(Xt.T @ XtW, XtW.T @ yt)
    r = yt - Xt @ b
    s2 = float((w * r * r).sum()) / (n - q)
    return delta, w, s2


def _wald_z(Xr, Yt, w, S1, Sy, Syy):
    """Vectorized single-regressor Wald z (intercept partialled out) under fixed weights w.

    Xr: 231 x b rotated (or raw) genotypes; Yt: 231 rotated (or raw) phenotype. Scalars
    S1,Sy,Syy are the weighted intercept sums for this phenotype. Returns z (b,)."""
    Sx = w @ Xr
    Sxx = w @ (Xr * Xr)
    Sxy = (w * Yt) @ Xr
    denom = Sxx - Sx * Sx / S1
    num = Sxy - Sx * Sy / S1
    beta = num / np.maximum(denom, 1e-12)
    syy = Syy - Sy * Sy / S1
    rss = syy - beta * num
    sigma2 = rss / (N - 2)
    se = np.sqrt(np.maximum(sigma2, 1e-30) / np.maximum(denom, 1e-12))
    z = beta / se
    z[denom <= 1e-10] = 0.0
    return z.astype(np.float32)


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    fit = pd.read_csv(f"{lib.GEA}/ecotype_fitness/ecotype_fitness.csv")
    Y = np.column_stack([qn(fit[p].to_numpy(float)) for p in PHENOS])  # 231 x P, rank-INT
    P = Y.shape[1]
    print(f"phenotypes ({P}, rank-inverse-normal): {PHENOS}", flush=True)

    # pass 1: per-chrom GRM contributions
    print("pass 1: kinship contributions", flush=True)
    Kc = {}; Mc = {}
    for cl in CHROMS:
        Kc[cl], Mc[cl] = kinship_contrib(cl)
        print(f"  {cl}: {Mc[cl]} SNPs ({time.time()-t0:.0f}s)", flush=True)
    Kraw = sum(Kc.values()); Mtot = sum(Mc.values())

    # naive (no-kinship) null constants: w=1, raw phenotype
    ones = np.ones(N)
    S1n = float(ones @ ones)
    Syn = ones @ Y                                  # P
    Syyn = ones @ (Y * Y)                           # P

    print("pass 2: per-chrom LOCO GWAS (LMM) + naive OLS", flush=True)
    Zout, Znout, CH, POS, RL, AL, MAC, VCL = [], [], [], [], [], [], [], []
    DELTA = np.zeros((len(CHROMS), P)); H2 = np.zeros((len(CHROMS), P))
    for ci, cl in enumerate(CHROMS):
        # LOCO kinship for this chrom
        Kl = (Kraw - Kc[cl]) / (Mtot - Mc[cl])
        d, Q = np.linalg.eigh(Kl); d = np.maximum(d, 1e-8)
        one_t = Q.T @ ones
        Yt = Q.T @ Y                                # 231 x P (rotated)
        Xt = one_t[:, None]
        W = np.zeros((N, P)); S1 = np.zeros(P); Sy = np.zeros(P); Syy = np.zeros(P)
        for j in range(P):
            delta, w, _ = _reml_delta(Yt[:, j], Xt, d)
            DELTA[ci, j] = delta; H2[ci, j] = 1.0 / (1.0 + delta)
            W[:, j] = w; S1[j] = w.sum(); Sy[j] = (w * Yt[:, j]).sum()
            Syy[j] = (w * Yt[:, j] * Yt[:, j]).sum()

        pos, rl, al, vp, vc, n_alt, n_cal = _load_chrom(cl)
        keep = (n_alt >= MIN_MAC) & (n_alt <= N - MIN_MAC) & (n_cal >= CALL_MIN * N)
        idx = np.where(keep)[0]
        dlen = np.abs(al - rl)
        vclass = np.where((rl == 1) & (al == 1), 0, np.where(dlen > 50, 2, 1))
        for s in range(0, len(idx), BLK):
            cols = idx[s:s + BLK]
            g, _ = _dense_geno(vp, vc, n_alt, n_cal, cols)
            Xr = Q.T @ g                            # rotated genotypes (LMM)
            b = Xr.shape[1]
            Zb = np.empty((b, P), np.float32); Zn = np.empty((b, P), np.float32)
            for j in range(P):
                Zb[:, j] = _wald_z(Xr, Yt[:, j], W[:, j], S1[j], Sy[j], Syy[j])
                Zn[:, j] = _wald_z(g, Y[:, j], ones, S1n, Syn[j], Syyn[j])
            Zout.append(Zb); Znout.append(Zn)
            CH.append(np.full(b, cl)); POS.append(pos[cols])
            RL.append(rl[cols]); AL.append(al[cols])
            MAC.append(np.minimum(n_alt[cols], N - n_alt[cols]).astype(int))
            VCL.append(vclass[cols])
        print(f"  {cl}: {len(idx)} variants, mean h2={H2[ci].mean():.2f} "
              f"({time.time()-t0:.0f}s)", flush=True)

    Z = np.vstack(Zout); Zn = np.vstack(Znout)
    np.savez(f"{OUT}/gwas_z.npz",
             chrom=np.concatenate(CH), pos=np.concatenate(POS),
             ref_len=np.concatenate(RL), alt_len=np.concatenate(AL),
             mac=np.concatenate(MAC), vclass=np.concatenate(VCL),
             Z=Z, Znaive=Zn, pheno_names=np.array(PHENOS),
             delta=DELTA, h2=H2, chroms=np.array(CHROMS), n_snp_K=Mtot)
    lam = {PHENOS[j]: [round(float(np.median(Z[:, j] ** 2) / stats.chi2.ppf(0.5, 1)), 3),
                       round(float(np.median(Zn[:, j] ** 2) / stats.chi2.ppf(0.5, 1)), 3)]
           for j in range(P)}
    json.dump({"lambda_gc_[loco,naive]": lam,
               "h2_mean": dict(zip(PHENOS, H2.mean(0).round(3).tolist())),
               "n_var": int(Z.shape[0]), "n_snp_K": int(Mtot)},
              open(f"{OUT}/gwas_meta.json", "w"), indent=2)
    print(f"\nn_var={Z.shape[0]:,}  n_snp_K={Mtot:,}")
    print("lambda_GC [LOCO, naive] per pheno:")
    for k, v in lam.items():
        print(f"  {k:16s} LOCO {v[0]:6.2f}   naive {v[1]:8.2f}")
    print(f"wrote {OUT}/gwas_z.npz  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()

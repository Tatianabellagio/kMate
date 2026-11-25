#!/usr/bin/env python
"""Burden LMM: founder fitness in GrENE-Net ~ variant burden + covariates + kinship random effect.

  y_f   = founder log-fitness: log((h0 + dh + eps)/(h0 + eps)), per garden and mean over gardens
          (founder_site_dh.npz, as _sv_founder_load_fitness.py)
  burden= RATE: variants of the class carried / variants of the class the founder is called at
          (founder_burden.csv; >= 2 carriers, variant called in >= 90% of founders)
  model = y ~ 1 + z(burden_test) + z(burden_ref) + callrate_sv + callrate_snp + cactus + u,  u ~ N(0, s2g K)
          burden_ref = same kind and origin, INTERGENIC -- absorbs total divergence from the
          reference, so burden_test measures the load specific to that context.
  fit   = REML via eigendecomposition of K (EMMA), Wald test on burden_test.
Null for the headline model: 1000 founder-label permutations of the burden_test residual after
regressing it on burden_ref (keeps the ref burden, breaks the test-context-specific part).
Env: kmate. Writes results/sv_adaptive/sv_burden_lmm{,_by_garden,_joint}.csv, sv_burden_lmm_perm.npz.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import optimize, stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
EF = f"{lib.GEA}/r3_persite_gwas/results/ecotype_fitness/founder_site_dh.npz"
EPS = 1e-5
CTX = ["CDS", "UTR", "intron", "TE", "promoter"]
UNIT = os.environ.get("SV_UNIT", "bubble")          # bubble | allele | event (SV burden unit)
TAG = {"bubble": "", "allele": "_allele", "event": "_event"}[UNIT]
if UNIT == "event" and os.environ.get("EVENT_MAP", "truvari") == "strict":
    TAG = "_eventstrict"
TAG += "" if float(os.environ.get("SV_CALL_MIN", 0.9)) == 0.9 else f"_cr{int(round(float(os.environ['SV_CALL_MIN']) * 100))}"


class LMM:
    def __init__(self, K):
        S, U = np.linalg.eigh(K); self.S = np.clip(S, 0, None); self.U = U

    def fit(self, y, X):
        Uy = self.U.T @ y; UX = self.U.T @ X; n, p = X.shape

        def nll(ld):
            d = np.exp(ld); w = 1 / (self.S + d)
            XtWX = UX.T @ (UX * w[:, None]); b = np.linalg.solve(XtWX, UX.T @ (Uy * w))
            r = Uy - UX @ b; s2 = (r * r * w).sum() / (n - p)
            return 0.5 * ((n - p) * np.log(s2) + np.log(self.S + d).sum() + np.linalg.slogdet(XtWX)[1])
        grid = np.linspace(-8, 8, 60); v = [nll(g) for g in grid]; g0 = grid[int(np.argmin(v))]
        ld = optimize.minimize_scalar(nll, bounds=(g0 - 0.3, g0 + 0.3), method="bounded").x
        d = np.exp(ld); w = 1 / (self.S + d)
        XtWX = UX.T @ (UX * w[:, None]); b = np.linalg.solve(XtWX, UX.T @ (Uy * w))
        r = Uy - UX @ b; s2 = (r * r * w).sum() / (n - p)
        se = np.sqrt(np.diag(np.linalg.inv(XtWX)) * s2)
        return b, se, 1 / (1 + d)                # h2-like share of kinship variance


def z(v):
    v = np.asarray(v, float); return (v - v.mean()) / v.std()


def main():
    D = pd.read_csv(f"{OUT}/founder_burden{TAG}.csv", dtype={"founder": str})
    K = np.load(f"{OUT}/founder_kinship.npy")
    Zf = np.load(EF, allow_pickle=True)
    assert np.array_equal(Zf["founders"].astype(str), D.founder.to_numpy())
    Y = np.log((Zf["H0"][None, :] + Zf["DH_rel"] + EPS) / (Zf["H0"][None, :] + EPS))   # sites x founders
    ymean = Y.mean(0); bio1_site = Zf["bio1"]
    L = LMM(K)
    cov = np.column_stack([np.ones(len(D)), z(D.callrate_sv), z(D.callrate_snp), D.is_cactus.astype(float)])
    rate = lambda c: (D[c] / D[f"called_{c}"]).to_numpy(float)

    rows = []
    for fq in ("all", "rare"):
        for kind in ("ins", "del", "indel", "snp"):
            for o in ("all", "shared", "cold", "warm"):
                ref = rate(f"{kind}_{o}_intergenic_{fq}")
                for x in CTX:
                    t = rate(f"{kind}_{o}_{x}_{fq}")
                    if t.std() == 0 or ref.std() == 0:
                        continue
                    X = np.column_stack([cov, z(t), z(ref)])
                    b, se, h2 = L.fit(ymean, X)
                    b0, se0, _ = L.fit(ymean, np.column_stack([cov, z(t)]))
                    rows.append(dict(freq=fq, kind=kind, origin=o, context=x, mean_burden=t.mean(),
                                     r_ref=np.corrcoef(t, ref)[0, 1],
                                     beta_alone=b0[-1], p_alone=2 * stats.norm.sf(abs(b0[-1] / se0[-1])),
                                     beta=b[-2], se=se[-2], p=2 * stats.norm.sf(abs(b[-2] / se[-2])),
                                     beta_ref=b[-1], p_ref=2 * stats.norm.sf(abs(b[-1] / se[-1])), h2=h2))
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/sv_burden_lmm{TAG}.csv", index=False)
    pd.set_option("display.width", 200)
    for fq in ("all", "rare"):
        print(f"\n=== freq={fq}: mean log-fitness ~ burden(context) + burden(intergenic) + callrate + cactus + K ===")
        print("beta per SD of context burden, conditional on intergenic burden of same kind/origin (p)")
        T = R[R.freq == fq].copy()
        T["cell"] = T.beta.map("{:+.3f}".format) + " (" + T.p.map("{:.2g}".format) + ")"
        print(T.pivot_table(index=["kind", "origin"], columns="context", values="cell", aggfunc="first")[CTX].to_string())

    # headline model, per garden and permutation null
    print("\n=== headline: shared CDS insertion burden (all freq) ===")
    t, ref = rate("ins_shared_CDS_all"), rate("ins_shared_intergenic_all")
    X = np.column_stack([cov, z(t), z(ref)])
    tr = t - np.polyval(np.polyfit(ref, t, 1), ref)
    b_obs = L.fit(ymean, X)[0][-2]
    rng = np.random.default_rng(0); null = []
    for _ in range(1000):
        tp = np.polyval(np.polyfit(ref, t, 1), ref) + rng.permutation(tr)
        null.append(L.fit(ymean, np.column_stack([cov, z(tp), z(ref)]))[0][-2])
    null = np.array(null)
    print(f"  beta {b_obs:+.4f}; permutation p = {(np.abs(null) >= abs(b_obs)).mean():.3f}")
    np.savez(f"{OUT}/sv_burden_lmm_perm{TAG}.npz", null=null, obs=b_obs)
    G = []
    for s in range(Y.shape[0]):
        b, se, h2 = L.fit(Y[s], X)
        G.append(dict(site=int(Zf["sites"][s]), bio1=bio1_site[s], zone=Zf["zone"][s], beta=b[-2], se=se[-2],
                      p=2 * stats.norm.sf(abs(b[-2] / se[-2])), h2=h2))
    G = pd.DataFrame(G); G.to_csv(f"{OUT}/sv_burden_lmm_by_garden{TAG}.csv", index=False)
    rho = stats.spearmanr(G.bio1, G.beta)
    print(f"  per garden: beta < 0 in {(G.beta < 0).sum()}/{len(G)}; nominal p<0.05 in {(G.p < 0.05).sum()}; "
          f"beta vs garden bio1 rho = {rho.statistic:+.2f} (p = {rho.pvalue:.2g})")
    # is the insertion effect independent of the SNP / indel burden of the same context?
    print("\n=== joint models, shared, all freq (rates) ===")
    J = []
    for lab, cols in (("ins + SNP", ["ins", "snp"]), ("ins + SNP + indel", ["ins", "snp", "indel"])):
        names = [f"{k}_shared_{x}_all" for k in cols for x in ("CDS", "intergenic")]
        b, se, _ = L.fit(ymean, np.column_stack([cov] + [z(rate(n)) for n in names]))
        J += [dict(model=lab, term=n, beta=bb, se=ss, p=2 * stats.norm.sf(abs(bb / ss)))
              for n, bb, ss in zip(names, b[-len(names):], se[-len(names):])]
        print(f" {lab}:"); [print(f"   {n:<28} {bb:+.3f}  p={2 * stats.norm.sf(abs(bb / ss)):.2g}")
                            for n, bb, ss in zip(names, b[-len(names):], se[-len(names):])]
    pd.DataFrame(J).to_csv(f"{OUT}/sv_burden_lmm_joint{TAG}.csv", index=False)
    for k in ("ins", "snp"):
        print(f"  corr({k} CDS rate, callrate_sv) = {np.corrcoef(rate(f'{k}_shared_CDS_all'), D.callrate_sv)[0, 1]:+.2f}")
    print(f"[wrote] {OUT}/sv_burden_lmm.csv, sv_burden_lmm_by_garden.csv")


if __name__ == "__main__":
    main()

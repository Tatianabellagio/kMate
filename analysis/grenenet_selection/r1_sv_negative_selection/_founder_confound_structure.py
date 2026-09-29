#!/usr/bin/env python
"""Structure + divergence controls for the founder insertion-load x climate correlation.

Follows `_founder_climate_confound.py`, which found that the founder climate response
gamma_f correlates with insertion FRACTION (bio1 rho=-0.537) and survives controlling for
founder home climate (partial -0.363), panel membership, and founder p0. Two weak links
in that result are tested here:

  (W1) The K_snp mixed model pinned its variance component at the grid boundary
       (delta=1e-5), i.e. REML wanted ALL variance in the kinship term. With a
       rank-deficient K that makes the fixed-effect SE unreliable. Diagnose the
       eigen-spectrum and the REML profile before trusting that p-value.

  (W2) "231 founders are ~10-25 clades" -- a Spearman over 231 structured points can be
       driven entirely by between-clade contrast (n_effective ~ 20), which is a far
       weaker claim than n=231 suggests. Test the correlation AT the clade level and
       WITHIN clades separately.

  (W3) Is insertion fraction just a proxy for genome-wide divergence from the TAIR10/Col-0
       reference? If divergent founders lose in hot gardens, "insertion load" is
       incidental. Internal control already present: deletion fraction shows NO
       correlation with gamma, which divergence alone cannot explain. Tested explicitly
       here by adding total carried-variant count as a covariate.

Env: kmate. Reads founder_climate_confound.npz, founder_load_test.npz, class_grms.npz.
Writes results/sv_adaptive/founder_confound_structure.csv + printed report.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

GEA = lib.GEA
OUT = f"{GEA}/r1_sv_negative_selection/results/sv_adaptive"
rng = np.random.default_rng(0)


def partial_spearman(x, y, Z):
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    RZ = np.column_stack([stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
    A = np.column_stack([np.ones(len(rx)), RZ])
    bx, *_ = np.linalg.lstsq(A, rx, rcond=None)
    by, *_ = np.linalg.lstsq(A, ry, rcond=None)
    ex = rx - A @ bx; ey = ry - A @ by
    r = float(np.corrcoef(ex, ey)[0, 1])
    dof = len(rx) - 2 - Z.shape[1]
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), dof))


def main():
    C = np.load(f"{OUT}/founder_climate_confound.npz", allow_pickle=True)
    F = np.load(f"{OUT}/founder_load_test.npz", allow_pickle=True)
    G = np.load(f"{GEA}/r3_persite_gwas/results/varexp/class_grms.npz", allow_pickle=True)

    keep = C["keep"].astype(bool)
    ins_frac = C["ins_frac"][keep]
    del_frac = C["del_frac"][keep]
    home1 = C["home_bio1"][keep]; home18 = C["home_bio18"][keep]
    g1 = C["gamma_bio1"][keep]; g18 = C["gamma_bio18"][keep]
    tot = F["tot_carry"][keep]
    K = G["K_snp"].astype(np.float64)
    K = K / np.mean(np.diag(K))
    Ks = K[np.ix_(keep, keep)]
    n = int(keep.sum())
    rows = []

    # ------------------------------------------------------------------ W1
    print("=" * 78)
    print("W1  Is the K_snp mixed model trustworthy? (it pinned delta at the boundary)")
    d = np.linalg.eigvalsh(Ks)[::-1]
    d = np.maximum(d, 0)
    print(f"  eigenvalues of K_snp[keep]: n={n}  top5={np.round(d[:5], 2)}")
    print(f"    #eig < 1e-8 : {int((d < 1e-8).sum())}   #eig < 1e-3 : {int((d < 1e-3).sum())}")
    print(f"    condition (d1/d_min_pos) = {d[0] / d[d > 0].min():.3g}")
    print(f"    variance in top 10 eigenvectors: {d[:10].sum() / d.sum():.3f}")

    # REML profile for the intercept-only model, to see if the boundary is real
    def reml_ll(y, X, delta):
        dd, U = np.linalg.eigh(Ks)
        dd = np.maximum(dd, 0)
        yt = U.T @ y; Xt = U.T @ X
        h = dd + delta
        Xh = Xt / h[:, None]
        A = Xt.T @ Xh
        beta = np.linalg.pinv(A) @ (Xh.T @ yt)
        r = yt - Xt @ beta
        rss = float((r * r / h).sum())
        sg2 = rss / (n - X.shape[1])
        _, logdetA = np.linalg.slogdet(A)
        return -0.5 * ((n - X.shape[1]) * np.log(sg2) + np.log(h).sum() + logdetA + (n - X.shape[1]))

    y = (g1 - g1.mean()) / g1.std()
    X0 = np.ones((n, 1))
    print("  REML profile over delta (intercept-only, gamma_bio1):")
    for dl in (1e-6, 1e-4, 1e-2, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0):
        print(f"    delta={dl:<8g} logL={reml_ll(y, X0, dl):+.3f}")
    print("  -> if logL keeps rising as delta falls, the VC is degenerate and the")
    print("     fixed-effect Wald p from that model should NOT be quoted.")

    # ------------------------------------------------------------------ W2
    print("\n" + "=" * 78)
    print("W2  Clade-level vs within-clade: is this a between-clade contrast?")
    dist = np.add.outer(np.diag(Ks), np.diag(Ks)) - 2 * Ks
    np.fill_diagonal(dist, 0.0)
    dist = np.maximum(dist, 0.0)
    Zl = linkage(squareform(dist, checks=False), method="average")

    for k in (5, 10, 15, 20, 25, 40):
        lab = fcluster(Zl, k, criterion="maxclust")
        uq = np.unique(lab)
        cm_ins = np.array([ins_frac[lab == c].mean() for c in uq])
        cm_g1 = np.array([g1[lab == c].mean() for c in uq])
        cm_g18 = np.array([g18[lab == c].mean() for c in uq])
        cm_h1 = np.array([home1[lab == c].mean() for c in uq])
        sz = np.array([(lab == c).sum() for c in uq])
        r1, p1 = stats.spearmanr(cm_ins, cm_g1)
        r18, p18 = stats.spearmanr(cm_ins, cm_g18)
        # within-clade: subtract clade means from both, then correlate
        wi = ins_frac - np.array([cm_ins[list(uq).index(c)] for c in lab])
        wg = g1 - np.array([cm_g1[list(uq).index(c)] for c in lab])
        rw, pw = stats.spearmanr(wi, wg)
        # clade-level partial, controlling clade-mean home bio1
        rp, pp = partial_spearman(cm_ins, cm_g1, cm_h1[:, None]) if k >= 10 else (np.nan, np.nan)
        print(f"  k={k:<3d} clades (max size {sz.max():>3d}, singletons {int((sz == 1).sum()):>2d}) | "
              f"BETWEEN bio1 rho={r1:+.3f} p={p1:.3f} | bio18 rho={r18:+.3f} p={p18:.3f} | "
              f"between|home rho={rp:+.3f} p={pp:.3f} | WITHIN bio1 rho={rw:+.3f} p={pw:.2e}")
        rows.append(dict(test="clade", k=k, between_bio1=r1, p_between=p1,
                         between_bio18=r18, within_bio1=rw, p_within=pw,
                         between_partial_home=rp, p_between_partial=pp))

    # ------------------------------------------------------------------ W3
    print("\n" + "=" * 78)
    print("W3  Is insertion fraction just genome-wide divergence?")
    for lbl, v in (("tot_carry (divergence)", tot), ("del_frac", del_frac)):
        r, p = stats.spearmanr(ins_frac, v)
        print(f"  Spearman(ins_frac, {lbl:22s}) = {r:+.3f}  p={p:.2e}")
    for nm, g in (("gamma_bio1", g1), ("gamma_bio18", g18)):
        r_i, p_i = stats.spearmanr(g, ins_frac)
        r_d, p_d = stats.spearmanr(g, del_frac)
        r_t, p_t = stats.spearmanr(g, tot)
        rp, pp = partial_spearman(ins_frac, g, np.column_stack([tot, del_frac]))
        rpa, ppa = partial_spearman(ins_frac, g, np.column_stack([tot, del_frac, home1, home18]))
        print(f"  {nm}: ins {r_i:+.3f} | del {r_d:+.3f} (p={p_d:.2f}) | tot {r_t:+.3f} (p={p_t:.2e})")
        print(f"      partial(ins | tot+del)            = {rp:+.3f}  p={pp:.2e}")
        print(f"      partial(ins | tot+del+home)       = {rpa:+.3f}  p={ppa:.2e}")
        rows.append(dict(test="divergence", axis=nm, rho_ins=r_i, rho_del=r_d, rho_tot=r_t,
                         partial_tot_del=rp, p_tot_del=pp,
                         partial_tot_del_home=rpa, p_tot_del_home=ppa))

    # ------------------------------------------------------- clade-permutation null
    print("\n" + "=" * 78)
    print("W2b Clade-restricted permutation null (preserves clade structure)")
    for k in (10, 20):
        lab = fcluster(Zl, k, criterion="maxclust")
        obs = stats.spearmanr(ins_frac, g1).correlation
        # permute gamma WITHIN clades only -> kills any within-clade signal, keeps
        # between-clade contrast intact. If obs stays extreme, signal is between-clade.
        null_w = []
        for _ in range(2000):
            perm = g1.copy()
            for c in np.unique(lab):
                idx = np.where(lab == c)[0]
                perm[idx] = rng.permutation(perm[idx])
            null_w.append(stats.spearmanr(ins_frac, perm).correlation)
        null_w = np.array(null_w)
        # permute whole clade labels (swap clade-mean assignments) -> tests between-clade
        p_w = float((np.abs(null_w) >= abs(obs)).mean())
        print(f"  k={k:<3d}: obs rho={obs:+.3f} | within-clade-permuted null "
              f"mean={null_w.mean():+.3f} sd={null_w.std():.3f} p={p_w:.4f}")
        rows.append(dict(test="perm_within_clade", k=k, obs=obs, null_mean=null_w.mean(),
                         null_sd=null_w.std(), p=p_w))

    pd.DataFrame(rows).to_csv(f"{OUT}/founder_confound_structure.csv", index=False)
    print(f"\n[wrote] {OUT}/founder_confound_structure.csv")


if __name__ == "__main__":
    main()

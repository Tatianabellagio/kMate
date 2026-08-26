#!/usr/bin/env python
"""Which measure of founder SV content tracks climate -- and is any of it just Col-0 divergence?

Follow-up to `_founder_climate_confound.py`, which used a single SV summary (`ins_frac`).
Two questions it could not answer:

  Q1  Why insertion FRACTION specifically? SV content can be counted as a NUMBER of SVs or
      as the number of BASE PAIRS affected, for insertions / deletions / both, raw or
      normalised by total variant load. Those are different variables and can disagree.

  Q2  How much of this is divergence from Col-0? The pangenome is TAIR10/Col-0-referenced,
      so an "insertion" is by definition sequence a founder has and Col-0 lacks. Insertion
      count is therefore partly a divergence statistic. If divergence from Col-0 is itself
      climate-structured, the SV reading could be redundant with it.

      NOTE on `ins_frac` specifically: tot_carry is ~77% SNPs (median 284k of 367k), so
      dividing by it ALREADY normalises insertion count by a divergence proxy. Part of the
      answer to Q2 is therefore built into the original variable -- which is exactly why the
      raw counts need to be tested beside it rather than assumed equivalent.

Tests every measure against (a) founder home bio1 and (b) the founder climate response
gamma, then partials out SNP divergence (n_snp = ALT alleles carried at SNPs) and home
climate.

Env: kmate. Reads founder_sv_content.npz + founder_climate_confound.npz.
Writes results/sv_adaptive/founder_sv_content_climate.csv + printed report.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"

MEASURES = [
    ("n_ins",       "SV insertions, count"),
    ("n_del",       "SV deletions, count"),
    ("n_sv",        "SVs, count"),
    ("bp_ins",      "SV insertions, bp"),
    ("bp_del",      "SV deletions, bp"),
    ("bp_sv",       "SVs, bp"),
    ("n_ins_frac",  "SV insertions, count / total load"),
    ("n_del_frac",  "SV deletions, count / total load"),
    ("n_sv_frac",   "SVs, count / total load"),
    ("bp_ins_frac", "SV insertions, bp / total load"),
    ("bp_del_frac", "SV deletions, bp / total load"),
    ("bp_sv_frac",  "SVs, bp / total load"),
    ("n_snp",       "SNP ALT count = divergence from Col-0"),
    ("tot_carry",   "all common variants carried"),
]


def partial_spearman(x, y, Z):
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    if Z is None:
        r = float(stats.spearmanr(x, y).correlation); dof = len(rx) - 2
    else:
        Z = np.atleast_2d(Z.T).T if Z.ndim == 1 else Z
        RZ = np.column_stack([stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
        A = np.column_stack([np.ones(len(rx)), RZ])
        bx, *_ = np.linalg.lstsq(A, rx, rcond=None); by, *_ = np.linalg.lstsq(A, ry, rcond=None)
        ex = rx - A @ bx; ey = ry - A @ by
        r = float(np.corrcoef(ex, ey)[0, 1]); dof = len(rx) - 2 - Z.shape[1]
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), dof))


def main():
    S = np.load(f"{OUT}/founder_sv_content.npz", allow_pickle=True)
    C = np.load(f"{OUT}/founder_climate_confound.npz", allow_pickle=True)
    assert list(S["founders"].astype("U6")) == list(C["founders"].astype("U6"))
    keep = C["keep"].astype(bool)
    g1 = C["gamma_bio1"][keep]; g18 = C["gamma_bio18"][keep]
    home1 = C["home_bio1"][keep]; home18 = C["home_bio18"][keep]
    M = {k: S[k][keep].astype(float) for k, _ in MEASURES}
    nsnp = M["n_snp"]
    HOME = np.column_stack([home1, home18])

    # sanity: is ins_frac here the same variable as the original ins_frac?
    orig = C["ins_frac"][keep]
    r_chk = stats.spearmanr(orig, M["n_ins_frac"]).correlation
    print(f"[check] rebuilt n_ins_frac vs original ins_frac: Spearman={r_chk:+.4f}\n")

    rows = []
    print("=" * 108)
    print("Q1  Which SV measure tracks climate?   (rho with founder HOME bio1, and with gamma_bio1)")
    print(f"{'measure':<34}{'home bio1':>13}{'gamma_bio1':>13}{'g|home':>11}{'g|nsnp':>11}{'g|home+nsnp':>14}")
    print("-" * 108)
    for k, lab in MEASURES:
        v = M[k]
        rh, ph = stats.spearmanr(v, home1)
        rg, pg = stats.spearmanr(v, g1)
        rgh, pgh = partial_spearman(v, g1, HOME)
        rgn, pgn = partial_spearman(v, g1, nsnp[:, None])
        rgb, pgb = partial_spearman(v, g1, np.column_stack([home1, home18, nsnp]))
        star = lambda p: "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else " "
        print(f"{lab:<34}{rh:>+9.3f}{star(ph):<4}{rg:>+9.3f}{star(pg):<4}"
              f"{rgh:>+8.3f}{star(pgh):<3}{rgn:>+8.3f}{star(pgn):<3}{rgb:>+10.3f}{star(pgb):<4}")
        rows.append(dict(measure=k, label=lab, rho_home_bio1=rh, p_home_bio1=ph,
                         rho_gamma=rg, p_gamma=pg, partial_home=rgh, p_partial_home=pgh,
                         partial_nsnp=rgn, p_partial_nsnp=pgn,
                         partial_home_nsnp=rgb, p_partial_home_nsnp=pgb))

    print("\n" + "=" * 108)
    print("Q2  Is the climate signal just divergence from Col-0?")
    rdh, pdh = stats.spearmanr(nsnp, home1)
    rdg, pdg = stats.spearmanr(nsnp, g1)
    print(f"  divergence (n_snp) vs home bio1 : rho={rdh:+.3f} p={pdh:.2e}")
    print(f"  divergence (n_snp) vs gamma_bio1: rho={rdg:+.3f} p={pdg:.2e}")
    rdgh, pdgh = partial_spearman(nsnp, g1, HOME)
    print(f"  divergence vs gamma | home      : rho={rdgh:+.3f} p={pdgh:.2e}")
    print("  (if divergence itself carried the climate signal, these would be large)")

    print("\n  How entangled is each SV measure with divergence?")
    for k, lab in MEASURES:
        if k in ("n_snp", "tot_carry"):
            continue
        r, p = stats.spearmanr(M[k], nsnp)
        print(f"    {lab:<34} vs n_snp: rho={r:+.3f} p={p:.1e}")

    print("\n" + "=" * 108)
    print("Insertions vs deletions, bp vs count -- the direct contrast")
    for k in ("n_ins_frac", "n_del_frac", "bp_ins_frac", "bp_del_frac"):
        r, p = stats.spearmanr(M[k], g1)
        rp, pp = partial_spearman(M[k], g1, np.column_stack([home1, home18, nsnp]))
        lab = dict(MEASURES)[k]
        print(f"  {lab:<36} gamma rho={r:+.3f} (p={p:.1e})  |  fully controlled {rp:+.3f} (p={pp:.1e})")

    pd.DataFrame(rows).to_csv(f"{OUT}/founder_sv_content_climate.csv", index=False)
    print(f"\n[wrote] {OUT}/founder_sv_content_climate.csv")


if __name__ == "__main__":
    main()

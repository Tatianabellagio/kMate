#!/usr/bin/env python
"""Could founder insertion content be an artifact of how the panel was BUILT?

The cold-origin -> more-inserted-sequence relationship (rho = -0.500 on bp_ins_frac) is a
statement about accessions, but every accession's SV content was *measured* by a pipeline.
Three ways that pipeline could manufacture the pattern:

  A1  PANEL HALF. 80 founders got their SVs from their own long-read assembly (cactus); 151
      were genotyped from short reads against the graph (PanGenie). These are completely
      different measurement processes. §1 of the notebook tested this on the insertion COUNT
      fraction only -- the headline measure is now insertion BP, which was never checked.

  A2  CALL QUALITY. A founder whose genotypes are better resolved can carry more of
      everything. Per-founder call rate over the common panel (and over SV records alone) is
      the direct covariate and is computable from var_called.

  A3  GRAPH CONTENT CEILING. A PanGenie founder can only carry insertions that some cactus
      assembly contributed to the graph. So PG insertion content is bounded by graph content
      in a way cactus insertion content is not. If the climate relationship holds WITHIN the
      PG half alone, it cannot be an assembly-quality effect -- those founders have no
      assembly. If it holds within the cactus half alone, it is not a genotyping effect.
      Holding in both is the strong result.

NOT TESTABLE HERE: the specific long-read technology (HiFi vs CLR vs ONT), assembly
contiguity (N50), and assembly length. None of that metadata is in this repo -- checked
`data/request_assemblies_for_Moi.csv` (Assembly_ID / Accession_ID / version flags only) and
the panel READMEs. Getting it means pulling the assembly metadata table from the source
pangenome release. That is the main open technical check and is stated as such.

Env: kmate. Reads founder_sv_content.npz + founder_climate_confound.npz.
Writes results/sv_adaptive/founder_sv_technical.csv + printed report.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
MEAS = [("bp_ins_frac", "ins, bp / load"), ("bp_ins", "ins, bp"),
        ("n_ins_frac", "ins, count / load"), ("n_ins", "ins, count"),
        ("bp_del_frac", "del, bp / load"), ("n_del_frac", "del, count / load")]


def partial_spearman(x, y, Z):
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    RZ = np.column_stack([stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
    A = np.column_stack([np.ones(len(rx)), RZ])
    bx, *_ = np.linalg.lstsq(A, rx, rcond=None); by, *_ = np.linalg.lstsq(A, ry, rcond=None)
    ex = rx - A @ bx; ey = ry - A @ by
    r = float(np.corrcoef(ex, ey)[0, 1]); dof = len(rx) - 2 - Z.shape[1]
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), dof))


def main():
    C = np.load(f"{OUT}/founder_climate_confound.npz", allow_pickle=True)
    S = np.load(f"{OUT}/founder_sv_content.npz", allow_pickle=True)
    keep = C["keep"].astype(bool)
    home1 = C["home_bio1"][keep]; home18 = C["home_bio18"][keep]
    g1 = C["gamma_bio1"][keep]
    cac = C["is_cactus"][keep].astype(bool)
    M = {k: S[k][keep].astype(float) for k in S.files if k != "founders"}
    cr_all = M["call_rate_all"]; cr_sv = M["call_rate_sv"]
    rows = []

    # ------------------------------------------------------------------ A1
    print("=" * 96)
    print("A1  Does the panel half predict SV content?  (§1 tested the COUNT fraction only)")
    for k, lab in MEAS:
        v = M[k]
        r, p = stats.spearmanr(v, cac.astype(float))
        mw = stats.mannwhitneyu(v[cac], v[~cac])
        print(f"  {lab:<20} cactus median {np.median(v[cac]):>12.5g} | PG {np.median(v[~cac]):>12.5g} "
              f"| rho={r:+.3f} p={p:.3f} MWU p={mw.pvalue:.3f}")
        rows.append(dict(test="A1_panel", measure=k, rho=r, p=p, mwu_p=mw.pvalue))

    # ------------------------------------------------------------------ A2
    print("\n" + "=" * 96)
    print("A2  Call quality -- do better-called founders carry more inserted sequence?")
    print(f"  call rate (all common): median {np.median(cr_all):.4f}  range {cr_all.min():.4f}-{cr_all.max():.4f}")
    print(f"  call rate (SV records): median {np.median(cr_sv):.4f}  range {cr_sv.min():.4f}-{cr_sv.max():.4f}")
    print(f"  cactus vs PG call rate (all): {np.median(cr_all[cac]):.4f} vs {np.median(cr_all[~cac]):.4f} "
          f"MWU p={stats.mannwhitneyu(cr_all[cac], cr_all[~cac]).pvalue:.2e}")
    for nm, cr in (("call_rate_all", cr_all), ("call_rate_sv", cr_sv)):
        r_ch, p_ch = stats.spearmanr(cr, home1)
        print(f"\n  --- {nm} ---   vs ORIGIN bio1: rho={r_ch:+.3f} p={p_ch:.2e}")
        for k, lab in MEAS:
            v = M[k]
            r, p = stats.spearmanr(v, cr)
            rp, pp = partial_spearman(v, home1, cr[:, None])
            print(f"    {lab:<20} vs {nm}: rho={r:+.3f} (p={p:.1e})   "
                  f"| ORIGIN-climate rho with {nm} held fixed: {rp:+.3f} (p={pp:.1e})")
            rows.append(dict(test=f"A2_{nm}", measure=k, rho=r, p=p,
                             partial_origin=rp, p_partial=pp))

    # ------------------------------------------------------------------ A3
    print("\n" + "=" * 96)
    print("A3  Does the ORIGIN-climate relationship hold WITHIN each panel half separately?")
    print("    (PG founders have no assembly at all -> if it holds there, not an assembly effect)")
    print(f"{'measure':<20}{'all 231':>18}{'cactus n=%d' % cac.sum():>20}{'PanGenie n=%d' % (~cac).sum():>22}")
    for k, lab in MEAS:
        v = M[k]
        ra, pa = stats.spearmanr(v, home1)
        rc, pc = stats.spearmanr(v[cac], home1[cac])
        rp_, pp_ = stats.spearmanr(v[~cac], home1[~cac])
        print(f"{lab:<20}{ra:>+11.3f} ({pa:>6.1e}){rc:>+11.3f} ({pc:>6.1e}){rp_:>+13.3f} ({pp_:>6.1e})")
        rows.append(dict(test="A3_within_panel", measure=k, rho_all=ra, p_all=pa,
                         rho_cactus=rc, p_cactus=pc, rho_pg=rp_, p_pg=pp_))

    print("\n  Same, for the SELECTION response gamma_bio1:")
    for k, lab in MEAS[:2]:
        v = M[k]
        ra = stats.spearmanr(v, g1)
        rc = stats.spearmanr(v[cac], g1[cac])
        rp_ = stats.spearmanr(v[~cac], g1[~cac])
        print(f"    {lab:<20} all {ra.correlation:+.3f} | cactus {rc.correlation:+.3f} "
              f"(p={rc.pvalue:.1e}) | PG {rp_.correlation:+.3f} (p={rp_.pvalue:.1e})")

    # full control including call rate
    print("\n" + "=" * 96)
    print("Headline measure with EVERY technical covariate held fixed at once")
    v = M["bp_ins_frac"]
    Z = np.column_stack([home1, home18, M["n_snp"], cr_all, cac.astype(float)])
    r, p = partial_spearman(v, g1, Z)
    print(f"  bp_ins_frac vs gamma | origin climate + Col-0 divergence + call rate + panel half:"
          f"  rho={r:+.3f} (p={p:.1e})")
    rows.append(dict(test="full_control", measure="bp_ins_frac", rho=r, p=p))

    pd.DataFrame(rows).to_csv(f"{OUT}/founder_sv_technical.csv", index=False)
    print(f"\n[wrote] {OUT}/founder_sv_technical.csv")


if __name__ == "__main__":
    main()

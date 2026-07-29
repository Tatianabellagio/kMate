#!/usr/bin/env python
"""Do the REAL-DATA block-WZA hits survive an honest (heavy-tail-aware) threshold?

Context
-------
`null_tail_shape.py` established that the production WZA reads block p-values off a
NORMAL reference (`1 - norm.cdf(Z_std)`), but within-block LD makes the null block-Z
HEAVY-TAILED (excess kurtosis ~0.86 for snp/nonsnp/smallindel; ~0.32 for SV). So the
nominal Bonferroni cutoff (~4.7-4.8 z) is ~170x too permissive, and the honest cutoff
from a signal-free site-permutation null + GPD tail fit is ~6.5 z (SV ~5.2). That was
all measured on the PERMUTATION null; it was never applied to the real hit lists.

This script closes that loop. It is pure post-processing of the production `_final`
WZA outputs -- NO refits -- and is exact by construction on two fronts:
  * MODEL-matched: it uses each model's OWN emitted block statistic. Production writes
    Z_pVal = norm.sf(Z_std), so Z_std = norm.isf(Z_pVal) recovers the exact production
    standardized block-Z (kendall / lfmm / quasibinomial each).
  * STANDARDISATION-matched: that Z_std already carries production's isotonic-SD +
    deg-5-clamped-mean correction. (Sanity: the real Z_std has sd~1.008 and q999~3.93,
    matching the isotonic-mean permutation null in null_tail_shape.csv to 2 dp -- so the
    mean-fit difference between the two is negligible in the tail, which is what licenses
    applying the permutation z_honest to the real Z.)

Per model x class x axis it counts, over the ~42-57k blocks:
  * nominal  : Z_std > norm.isf(0.05 / n)                      (the current production rule)
  * honest_perm : Z_std > z_honest[class]  from null_tail_shape.csv  (signal-free perm null)
  * honest_emp  : Z_std > z_emp, where z_emp is a per-axis peaks-over-threshold GPD
                  extrapolation of THIS axis's own block-Z tail to the 1-0.05/n quantile
                  (RepAdapt/Booker-style empirical null; model-matched by construction;
                  conservative if the axis carries real signal, exact under the null).
The two honest references are independent (one external+signal-free, one self+empirical);
agreement between them makes the verdict robust.

Writes honest_hits_check.csv (per model x class summary) and honest_hits_survivors.csv
(every block that clears EITHER honest threshold, so we can eyeball whether the handful
of survivors are real or noise-floor).
"""
from __future__ import annotations
import os, glob, re, warnings
import numpy as np
import pandas as pd
from scipy.stats import norm, genpareto

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.dirname(HERE)
WZA = f"{GEA}/phase1_replication/results/multiaxis/wza"

# production file naming uses model tag 'binomial' for the quasibinomial run
MODELS = ["kendall", "lfmm", "binomial"]
CLASSES = ["snp", "sv", "smallindel", "nonsnp"]
FNAME = re.compile(r"wza_(\w+?)_(snp|sv|smallindel|nonsnp)_gen9_(bio\d+)_final\.csv$")
TAIL_U_Q = 0.99          # peaks-over-threshold anchor for the per-axis empirical GPD
PMIN = 1e-300


def load_perm_honest():
    """z_honest / z_nominal per class from the signal-free permutation null."""
    d = pd.read_csv(f"{HERE}/null_tail_shape.csv").set_index("cls")
    return d["z_honest"].to_dict(), d["z_nominal"].to_dict(), d["excess_kurtosis"].to_dict()


def gpd_threshold(z, n, alpha):
    """Peaks-over-threshold GPD extrapolation of z's tail to the (1-alpha) quantile."""
    u = np.quantile(z, TAIL_U_Q)
    exc = z[z > u] - u
    if len(exc) < 25:
        return np.nan
    shape, _, scale = genpareto.fit(exc, floc=0)
    # P(Z > u+y) = (1-TAIL_U_Q)*SF(y); solve for the alpha-quantile
    return u + genpareto.ppf(1 - alpha / (1 - TAIL_U_Q), shape, loc=0, scale=scale)


def main():
    z_hon, z_nom_cls, kurt = load_perm_honest()
    rows, survivors = [], []

    for f in sorted(glob.glob(f"{WZA}/wza_*_gen9_*_final.csv")):
        m = FNAME.search(os.path.basename(f))
        if not m:
            continue
        model, cls, axis = m.group(1), m.group(2), m.group(3)
        if model not in MODELS or cls not in CLASSES:
            continue
        w = pd.read_csv(f)
        if "Z_pVal" not in w.columns or not len(w):
            continue
        w = w[w["Z_pVal"].notna()].copy()
        n = len(w)
        z = norm.isf(w["Z_pVal"].clip(PMIN, 1 - 1e-16).to_numpy())
        alpha = 0.05 / n

        z_nom = norm.isf(alpha)                 # per-file nominal Bonferroni z
        z_hp = z_hon[cls]                        # signal-free permutation honest z
        z_he = gpd_threshold(z, n, alpha)        # per-axis empirical honest z

        m_nom = z > z_nom
        m_hp = z > z_hp
        m_he = z > z_he if np.isfinite(z_he) else np.zeros(n, bool)
        rows.append(dict(model=model, cls=cls, axis=axis, n_blocks=n,
                         z_nominal=z_nom, z_honest_perm=z_hp, z_honest_emp=z_he,
                         max_z=float(np.nanmax(z)),
                         hits_nominal=int(m_nom.sum()),
                         hits_honest_perm=int(m_hp.sum()),
                         hits_honest_emp=int(m_he.sum())))
        keep = m_nom & (m_hp | m_he)             # nominal hits, flag which honest bars they clear
        for i in np.where(keep)[0]:
            r = w.iloc[i]
            survivors.append(dict(model=model, cls=cls, axis=axis, n_blocks=n,
                                  chrom=r.get("chrom"), pos=r.get("pos"),
                                  SNPs=r.get("SNPs"), z=float(z[i]),
                                  z_nominal=z_nom, z_honest_perm=z_hp, z_honest_emp=z_he,
                                  clears_perm=bool(z[i] > z_hp),
                                  clears_emp=bool(np.isfinite(z_he) and z[i] > z_he)))

    D = pd.DataFrame(rows)
    D.to_csv(f"{HERE}/honest_hits_check.csv", index=False)
    S = pd.DataFrame(survivors)
    S.to_csv(f"{HERE}/honest_hits_survivors.csv", index=False)

    # per model x class: sum hits over the 19 axes
    g = (D.groupby(["model", "cls"])
           .agg(axes=("axis", "nunique"),
                nominal=("hits_nominal", "sum"),
                honest_perm=("hits_honest_perm", "sum"),
                honest_emp=("hits_honest_emp", "sum"),
                max_z=("max_z", "max")).reset_index())
    g["kurtosis"] = g["cls"].map(kurt).round(2)
    g["z_hon_perm"] = g["cls"].map(z_hon).round(2)
    order = {c: i for i, c in enumerate(CLASSES)}
    g = g.sort_values(["model", "cls"], key=lambda s: s.map(order) if s.name == "cls" else s)

    print("REAL-DATA block-WZA hits, summed over 19 climate axes, per model x class")
    print("(nominal = current production rule; honest_perm = signal-free permutation null;")
    print(" honest_emp = per-axis empirical GPD self-null)\n")
    print(g[["model", "cls", "axes", "nominal", "honest_perm", "honest_emp",
             "max_z", "z_hon_perm", "kurtosis"]].to_string(index=False))
    print(f"\nwrote {HERE}/honest_hits_check.csv  ({len(D)} model x class x axis rows)")
    print(f"wrote {HERE}/honest_hits_survivors.csv  ({len(S)} surviving nominal hits)")
    if len(S):
        print("\nnominal hits that clear EITHER honest bar:")
        print(S.sort_values("z", ascending=False)[
            ["model", "cls", "axis", "chrom", "pos", "z",
             "z_honest_perm", "z_honest_emp", "clears_perm", "clears_emp"]
        ].round(3).to_string(index=False))


if __name__ == "__main__":
    main()

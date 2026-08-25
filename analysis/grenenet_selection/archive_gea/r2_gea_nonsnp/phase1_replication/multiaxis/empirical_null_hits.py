#!/usr/bin/env python
"""Empirical-null block hit-calling for the multiaxis WZA -- replaces the NORMAL
reference with the heavy-tailed empirical null, then BH-FDR as in production.

Why
---
Production WZA writes a block p-value off a NORMAL reference
(`wza_script.py`: Z_pVal = norm.sf(Z_std)). But our windows are BigLD clq0.9
LD-islands: within-block LD makes the null block-Z HEAVY-TAILED (excess kurtosis
~0.86 for snp/nonsnp/smallindel; ~0.32 for SV -- `null_tail_shape.csv`). Booker
et al. (2024, Mol. Ecol. Resour.) document this exact deviation and attribute it
to within-window linkage; they only ever claim the parametric p is *quasi*-uniform
and validate with ranking metrics (AUC-PR) at ~1000 windows -- NOT extreme-tail
calibration at our ~42-57k blocks x genome-wide alpha. So the Normal cutoff
(z~4.8) is ~170x too permissive here. The fix is on the null side, not the
windowing side (reducing within-block LD would also kill WZA's power).

What
----
Pure post-processing of the 240 production `_final` WZA files -- NO refits.
For each (model x class x axis) file:
  * recover the EXACT production standardized block-Z:  z_std = norm.isf(Z_pVal)
    (model- and standardisation-matched by construction; carries the isotonic-SD +
     deg-5-clamped-mean correction already).
  * EMPIRICAL block p-value (`p_emp`) from the genome-wide block-Z self-null:
      - peaks-over-threshold GPD fit above the 99th pct anchor u for the smooth
        extreme tail (RepAdapt/Booker empirical-p style),
      - plain right-tail ECDF at/below u.
    This is calibrated to the ACTUAL heavy tail; conservative if the axis carries
    real signal (the signal sits in its own tail), exact under the null.
  * BH-FDR on p_emp  ->  q_emp     (the production decision rule, honest null)
  * for reference: BH-FDR on the nominal Normal Z_pVal -> q_nom, plus the
    signal-free permutation gate  z_std > z_honest[class]  (`null_tail_shape.csv`).

A block is an EMPIRICAL HIT if q_emp < --fdr. We also flag Bonferroni-empirical
(p_emp < 0.05/n) and the perm gate, so the two independent honest nulls
(self-empirical + external-permutation) can be compared -- agreement makes the
verdict robust.

Outputs (results/multiaxis/empirical_null/):
  blocks_empirical_p.csv.gz   every block x combo: z_std, p_nom, p_emp, q_nom,
                              q_emp, hit_* flags  (the reusable honest table)
  hits_empirical.csv          only blocks that clear ANY honest bar, + span/genes
  summary_empirical.csv       per model x class: nominal vs empirical hit counts

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY phase1_replication/multiaxis/empirical_null_hits.py --fdr 0.05
"""
from __future__ import annotations
import argparse, glob, os, re, sys, warnings
import numpy as np
import pandas as pd
from scipy.stats import norm, genpareto
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import lib

warnings.filterwarnings("ignore")
GEA = lib.GEA
WZA = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza"
OUTDIR = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/empirical_null"
NULL_TAIL = f"{GEA}/wza/investigation/null_tail_shape.csv"

MODELS = ["kendall", "lfmm", "binomial"]
CLASSES = ["snp", "sv", "smallindel", "nonsnp"]
FNAME = re.compile(r"wza_(\w+?)_(snp|sv|smallindel|nonsnp)_gen9_(bio\d+|pc1)_final\.csv$")
TAIL_U_Q = 0.99          # peaks-over-threshold anchor
PMIN = 1e-300


def bh_fdr(p):
    """Benjamini-Hochberg q-values (same implementation as the production caller)."""
    p = np.asarray(p, float); n = np.isfinite(p).sum()
    q = np.full(len(p), np.nan)
    ok = np.where(np.isfinite(p))[0]
    if not len(ok):
        return q
    o = ok[np.argsort(p[ok])]
    ranked = p[o] * n / (np.arange(1, len(o) + 1))
    q[o] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    return q


def empirical_block_p(z):
    """Right-tail empirical p-value for each block-Z against the genome-wide self-null.

    ECDF right tail below the 99th-pct anchor u; GPD peaks-over-threshold
    extrapolation above u so the extreme blocks get a smooth p rather than a
    discrete 1/n floor. Continuous at u (both ~ 1-TAIL_U_Q).
    """
    z = np.asarray(z, float)
    n = len(z)
    zs = np.sort(z)
    # ECDF right-tail: P(Z >= z_i) = (# >= z_i)/n
    p = (n - np.searchsorted(zs, z, side="left")) / n
    u = np.quantile(z, TAIL_U_Q)
    exc = z[z > u] - u
    if len(exc) >= 25:
        c, _, scale = genpareto.fit(exc, floc=0)
        tail = z > u
        # P(Z > z) = (1-TAIL_U_Q) * SF_GPD(z-u)
        p_tail = (1 - TAIL_U_Q) * genpareto.sf(z[tail] - u, c, loc=0, scale=scale)
        p[tail] = np.clip(p_tail, PMIN, 1.0)
    return np.clip(p, PMIN, 1.0)


def load_perm_honest():
    if not os.path.exists(NULL_TAIL):
        return {}
    d = pd.read_csv(NULL_TAIL).set_index("cls")
    return d["z_honest"].to_dict()


def block_spans(r2: float = 0.9):
    """chrom/start/end/n_variants for every clq{r2} block; id 'Chr{n}_{idx}'."""
    tag = f"clq{r2}"; rows = []
    for ci in range(1, 6):
        c = f"Chr{ci}"
        f = f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv"
        if not os.path.exists(f):
            continue
        g = pd.read_csv(f, sep="\t").sort_values("start_pos").reset_index(drop=True)
        for idx, row in g.iterrows():
            rows.append((f"{c}_{idx}", c, int(row.start_pos), int(row.end_pos),
                         int(row.n_variants)))
    return pd.DataFrame(rows, columns=["block", "chrom", "start", "end",
                                       "n_variants"]).set_index("block")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--annotate", action="store_true",
                    help="map empirical hits to overlapping TAIR10 genes (local GFF)")
    ap.add_argument("--dump-blocks", action="store_true",
                    help="also write the full ~10M-row per-block p/q table (slow, large)")
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    z_hon = load_perm_honest()

    rows = []
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
        z_std = norm.isf(w["Z_pVal"].clip(PMIN, 1 - 1e-16).to_numpy())
        p_nom = w["Z_pVal"].to_numpy()
        p_emp = empirical_block_p(z_std)
        alpha = args.fdr / n
        zh = z_hon.get(cls, np.nan)
        sub = pd.DataFrame(dict(
            model=model, cls=cls, axis=axis, block=w["index"].astype(str),
            chrom=w.get("chrom"), pos=w.get("pos"), SNPs=w.get("SNPs"),
            z_std=z_std, p_nom=p_nom, p_emp=p_emp,
            q_nom=bh_fdr(p_nom), q_emp=bh_fdr(p_emp), n_blocks=n))
        sub["hit_nom_bh"] = sub["q_nom"] < args.fdr
        sub["hit_emp_bh"] = sub["q_emp"] < args.fdr
        sub["hit_emp_bonf"] = sub["p_emp"] < alpha
        sub["hit_perm"] = sub["z_std"] > zh if np.isfinite(zh) else False
        rows.append(sub)

    D = pd.concat(rows, ignore_index=True)
    if args.dump_blocks:
        D.to_csv(f"{OUTDIR}/blocks_empirical_p.csv.gz", index=False, compression="gzip")

    # per model x class summary, hits summed over axes
    g = (D.groupby(["model", "cls"]).agg(
            axes=("axis", "nunique"),
            nom_bh=("hit_nom_bh", "sum"),
            emp_bh=("hit_emp_bh", "sum"),
            emp_bonf=("hit_emp_bonf", "sum"),
            perm=("hit_perm", "sum"),
            max_z=("z_std", "max")).reset_index())
    order = {c: i for i, c in enumerate(CLASSES)}
    g = g.sort_values(["model", "cls"],
                      key=lambda s: s.map(order) if s.name == "cls" else s)
    g.to_csv(f"{OUTDIR}/summary_empirical.csv", index=False)

    # any block clearing any honest bar -> hits table (for gene follow-up)
    H = D[D["hit_emp_bh"] | D["hit_emp_bonf"] | D["hit_perm"]].copy()
    if args.annotate and len(H):
        spans = block_spans(); genes = lib.load_genes()
        gnames, spans_txt = [], []
        for blk in H["block"]:
            if blk in spans.index:
                sp = spans.loc[blk]
                gc = genes[(genes.chrom == sp.chrom) & (genes.end >= sp.start)
                           & (genes.start <= sp.end)]
                gnames.append(";".join(gc.gene)); spans_txt.append(f"{sp.chrom}:{sp.start}-{sp.end}")
            else:
                gnames.append(""); spans_txt.append("")
        H["region"] = spans_txt; H["genes"] = gnames
    H.sort_values("z_std", ascending=False).to_csv(f"{OUTDIR}/hits_empirical.csv", index=False)

    print("multiaxis block-WZA hits, summed over axes, per model x class")
    print("nom_bh = BH on the NORMAL Z_pVal (current production rule)")
    print("emp_bh = BH on the empirical heavy-tail p (honest production rule)")
    print("emp_bonf = empirical p < 0.05/n | perm = signal-free permutation gate\n")
    print(g[["model", "cls", "axes", "nom_bh", "emp_bh", "emp_bonf", "perm",
             "max_z"]].to_string(index=False))
    print(f"\nblocks total: {len(D):,} | any-honest-bar hits: {len(H)}")
    if args.dump_blocks:
        print(f"wrote {OUTDIR}/blocks_empirical_p.csv.gz")
    print(f"wrote {OUTDIR}/summary_empirical.csv")
    print(f"wrote {OUTDIR}/hits_empirical.csv")
    if len(H):
        cols = ["model", "cls", "axis", "block", "z_std", "q_emp", "hit_emp_bh",
                "hit_emp_bonf", "hit_perm"] + (["region", "genes"] if args.annotate else [])
        print("\nblocks clearing any honest bar:")
        print(H.sort_values("z_std", ascending=False)[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()

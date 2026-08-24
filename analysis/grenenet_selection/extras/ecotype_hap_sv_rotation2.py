"""Decisive: SV-enrichment of winning haplotypes under a null that controls BOTH frequency AND
space at once, plus a spatial-dispersion diagnostic to explain the disagreement with the
JOINT-axis rotation audit.

My first rotation null lost the k (founder-count) control (rot null r2 0.0235 < exact-k null
0.027). Fix: residualize tag on k (subtract the qbin(k,8) mean -> tag_resid removes the
frequency confound), THEN apply the genome-rotation null to tag_resid (removes the spatial
confound). If winning haplotypes still have positive residual tag beyond a random spatial phase
-> real, robust to both. If it collapses to ~0 -> the survival was the k-confound leaking and
the spatial audit is right.

Also: DISPERSION diagnostic — is my ecotype-fitness selected set spatially CONTIGUOUS (like a
JOINT sweep, where the spatial confound bites) or DISPERSED across the genome (ecotype sorting,
where it does not)? Adjacency = fraction of selected haplotypes whose genomic neighbor is also
selected, vs the expected fraction q.

Reuses hap_sv_table.csv. Output -> hap_sv_rotation2.csv. Env: kmate.
"""
from __future__ import annotations
import os, sys, json
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

G = f"{lib.GEA}/r3_persite_gwas/results/ecotype_fitness/gwas"
# SFX matches ecotype_hap_sv_enrichment.py: "" = all-class clq90, "_clq90nosv" = SNP+indel-only
# clustering (self-tagging control). Reads hap_sv_table{SFX}.csv, writes hap_sv_rotation2{SFX}.*
MEMB_TAG = os.environ.get("HAPMEMB_TAG", "clq90")
SFX = "" if MEMB_TAG == "clq90" else f"_{MEMB_TAG}"
NROT = 5000
CELLS = [("zone_union", "fit_zone", 0.005), ("zone_union", "fit_zone", 0.01),
         ("global", "fit_glob", 0.01)]


def qbin(x, n):
    r = pd.Series(x).rank(method="first")
    return pd.qcut(r, min(n, r.nunique()), labels=False, duplicates="drop").to_numpy()


def rot_null_resid(resid, chrom, sel_idx, nrot=NROT, seed=0):
    rng = np.random.default_rng(seed)
    obs = resid[sel_idx].mean()
    order = np.arange(len(resid))
    idx_by_c = {c: order[chrom == c] for c in np.unique(chrom)}
    null = np.empty(nrot)
    for r in range(nrot):
        rr = resid.copy()
        for c, ix in idx_by_c.items():
            if len(ix) > 1:
                rr[ix] = np.roll(resid[ix], rng.integers(1, len(ix)))
        null[r] = rr[sel_idx].mean()
    p = (1 + (null >= obs).sum()) / (nrot + 1)
    return float(obs), float(null.mean()), float(null.std()), float(p)


def adjacency(sel_mask, chrom):
    """fraction of selected haplotypes whose immediate genomic neighbor is also selected."""
    hits = 0; tot = 0
    order = np.arange(len(sel_mask))
    for c in np.unique(chrom):
        ix = order[chrom == c]
        s = sel_mask[ix]
        if s.sum() == 0:
            continue
        nb = np.zeros(len(s), bool)
        nb[:-1] |= s[1:]; nb[1:] |= s[:-1]
        hits += (s & nb).sum(); tot += s.sum()
    return hits / tot if tot else np.nan


def main():
    H = pd.read_csv(f"{G}/hap_sv_table{SFX}.csv").sort_values(
        ["chrom", "start", "k"]).reset_index(drop=True)
    chrom = H.chrom.to_numpy()
    kbin = qbin(H.k.to_numpy(), 8)
    rows = []
    for axis, fcol, q in CELLS:
        fc = H[fcol].to_numpy()
        sel_mask = fc >= np.quantile(fc, 1 - q)
        sel_idx = np.where(sel_mask)[0]
        adj = adjacency(sel_mask, chrom)
        for mac, tname in ((12, "tag12"), (24, "tag24"), (46, "tag46")):
            tag = H[tname].to_numpy()
            # residualize tag on founder-count bin (remove frequency confound)
            resid = tag - pd.Series(tag).groupby(kbin).transform("mean").to_numpy()
            obs, nm, sd, p = rot_null_resid(resid, chrom, sel_idx)
            z = (obs - nm) / sd if sd > 0 else np.nan
            rows.append(dict(axis=axis, top_q=q, sv_mac=mac, n_sel=len(sel_idx),
                             adj_obs=round(adj, 3), adj_exp=round(q, 3),
                             resid_obs=round(obs, 5), rot_null=round(nm, 5),
                             z=round(z, 2), p_rot=round(p, 4)))
    out = pd.DataFrame(rows)
    out.to_csv(f"{G}/hap_sv_rotation2{SFX}.csv", index=False)
    json.dump(rows, open(f"{G}/hap_sv_rotation2{SFX}.json", "w"), indent=2)
    print("=== k-RESIDUALIZED tag under SPATIAL rotation null (controls BOTH freq & space) ===")
    print("resid_obs = winning-hap mean tag ABOVE its founder-count expectation;")
    print("rot_null ~ 0 if no spatial coincidence.  z = (obs-null)/null_sd.\n")
    print(out.to_string(index=False))
    print("\nadj_obs vs adj_exp: if adj_obs ~ adj_exp => selected set is spatially DISPERSED")
    print("(ecotype sorting, spatial confound does NOT bite); >> => contiguous (it does).")


if __name__ == "__main__":
    main()

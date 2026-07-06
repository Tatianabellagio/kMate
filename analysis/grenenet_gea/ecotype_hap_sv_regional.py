"""Last nail: does the winning-haplotype SV-enrichment (~1.4x on the kinship-corrected fitness
axis, exact-k null) survive REGIONAL matching too? Adds pericentromere distance, TE density and
gene density to the null — mirroring the other agent's confound controls, but on the kinship-
corrected garden-fitness axis (not JOINT) and the hap-cluster unit.

Reuses results/.../gwas/hap_sv_table.csv (chrom,start,end,k,tag12/24/46,fit_zone,fit_glob).
Strata scheme copied from _sv_hap_context.py: control draws come from the same
qbin(k,8) x qbin(covariate,5) stratum. arms_only drops clusters within 3 Mb of a centromere.
Output -> results/.../gwas/hap_sv_regional.csv. Env: kmate.
"""
from __future__ import annotations
import os, sys, json
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

G = f"{lib.GEA}/ecotype_fitness/gwas"
GFF = os.path.expanduser("~/ara_key_files/TAIR10_GFF3_genes_transposons.gff")
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
CEN = {"Chr1": 15_086_000, "Chr2": 3_607_000, "Chr3": 13_588_000,
       "Chr4": 3_956_000, "Chr5": 11_726_000}
CLENS = {"Chr1": 30_427_671, "Chr2": 19_698_289, "Chr3": 23_459_830,
         "Chr4": 18_585_056, "Chr5": 26_975_502}
NPERM = 5000
CELLS = [("zone_union", "fit_zone", 0.005), ("global", "fit_glob", 0.01)]


def qbin(x, n):
    r = pd.Series(x).rank(method="first")
    return pd.qcut(r, min(n, r.nunique()), labels=False, duplicates="drop").to_numpy()


def annotate(units):
    gff = pd.read_csv(GFF, sep="\t", header=None, comment="#",
                      names=["chrom", "src", "feat", "start", "end", "sc", "st", "fr", "attr"])
    gff = gff[gff.chrom.isin(CHROMS)]
    genes = gff[gff.feat == "gene"]; tes = gff[gff.feat == "transposable_element"]
    te_frac = np.zeros(len(units)); gene_n = np.zeros(len(units)); pcen = np.zeros(len(units))
    for ch in CHROMS:
        clen = CLENS[ch]
        sel = np.where(units.chrom.to_numpy() == ch)[0]
        if len(sel) == 0:
            continue
        covte = np.zeros(clen + 2, np.int8)
        for a, b in zip(tes[tes.chrom == ch].start, tes[tes.chrom == ch].end):
            covte[max(a, 0):min(b, clen) + 1] = 1
        cte = np.concatenate([[0], np.cumsum(covte)])
        ge = genes[genes.chrom == ch]
        gstart = np.sort(ge.start.to_numpy()); gend = np.sort(ge.end.to_numpy())
        for i in sel:
            s, e = int(units.start.iat[i]), int(units.end.iat[i]); ln = max(e - s, 1)
            te_frac[i] = (cte[min(e, clen) + 1] - cte[max(s, 0)]) / ln
            gene_n[i] = np.searchsorted(gstart, e, "right") - np.searchsorted(gend, s, "left")
            pcen[i] = abs((s + e) / 2 - CEN[ch])
    units = units.copy()
    units["te_frac"] = te_frac; units["gene_n"] = gene_n; units["pcen"] = pcen
    return units


def strata_null(tag, sel, strata, rng, nperm=NPERM):
    pools = {}
    for s in np.unique(strata):
        if s < 0:
            continue
        pools[s] = np.where(strata == s)[0]
    obs = tag[sel].mean()
    null = np.empty((len(sel), nperm))
    for i, j in enumerate(sel):
        pool = pools.get(strata[j])
        null[i] = tag[pool][rng.integers(0, len(pool), nperm)] if pool is not None and len(pool) else 0.0
    nm = null.mean(0)
    fold = obs / nm.mean() if nm.mean() > 0 else np.nan
    p = (1 + (nm >= obs).sum()) / (nperm + 1)
    return float(obs), float(nm.mean()), float(fold), float(p)


def main():
    H = pd.read_csv(f"{G}/hap_sv_table.csv")
    uni = H.drop_duplicates(["chrom", "start", "end"])[["chrom", "start", "end"]].reset_index(drop=True)
    uni = annotate(uni)
    key = uni.set_index(["chrom", "start", "end"])
    idx = list(zip(H.chrom, H.start, H.end))
    for c in ("te_frac", "gene_n", "pcen"):
        H[c] = key[c].reindex(idx).to_numpy()
    k = H.k.to_numpy()
    kbin = qbin(k, 8)
    arm = H.pcen.to_numpy() >= 3e6
    rng = np.random.default_rng(0)
    rows = []
    for axis, fcol, q in CELLS:
        fc = H[fcol].to_numpy()
        sel_all = np.where(fc >= np.quantile(fc, 1 - q))[0]
        sel_arm = sel_all[arm[sel_all]]
        for mac, tname in ((12, "tag12"), (24, "tag24"), (46, "tag46")):
            tag = H[tname].to_numpy()
            controls = {
                "k_only": kbin,
                "k+pericentro": kbin * 10 + qbin(H.pcen.to_numpy(), 5),
                "k+TE_frac": kbin * 10 + qbin(H.te_frac.to_numpy(), 5),
                "k+gene_n": kbin * 10 + qbin(H.gene_n.to_numpy(), 5),
            }
            for name, strata in controls.items():
                obs, nm, fold, p = strata_null(tag, sel_all, strata, rng)
                rows.append(dict(axis=axis, top_q=q, sv_mac=mac, control=name,
                                 n_sel=len(sel_all), fold=round(fold, 3), p=round(p, 4)))
            # arms-only
            strata_a = np.where(arm, kbin, -1)
            obs, nm, fold, p = strata_null(tag, sel_arm, strata_a, rng)
            rows.append(dict(axis=axis, top_q=q, sv_mac=mac, control="arms_only",
                             n_sel=len(sel_arm), fold=round(fold, 3), p=round(p, 4)))
    out = pd.DataFrame(rows)
    out.to_csv(f"{G}/hap_sv_regional.csv", index=False)
    json.dump(rows, open(f"{G}/hap_sv_regional.json", "w"), indent=2)
    print("covariate (top-0.5% zone winners vs all):")
    z = H.iloc[np.where(H.fit_zone.to_numpy() >= np.quantile(H.fit_zone, 0.995))[0]]
    for c in ("te_frac", "gene_n", "pcen"):
        print(f"  {c}: {z[c].mean():.3g} vs {H[c].mean():.3g}")
    print(f"\narms fraction of top-0.5% zone winners: {arm[np.where(H.fit_zone.to_numpy()>=np.quantile(H.fit_zone,0.995))[0]].mean()*100:.0f}%")
    print("\n=== SV-enrichment of winning haplotypes vs REGION+freq-matched null (kinship axis) ===")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()

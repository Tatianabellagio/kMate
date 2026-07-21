#!/usr/bin/env python
"""GO / stress-term over-representation of the NON-SNP-only (kMate-unique) candidate genes.

Question: is the set of genes under blocks the non-SNP (indel+SV) scan flags but the SNP scan
misses enriched for stress / cold / heat / ABA / flowering / circadian function -- BEYOND what
the SNP scan already flags, and beyond what you'd expect from the size of the gene list?

Two things make eyeballing a 969-gene list misleading and are corrected here:
  1. UNIVERSE. GWAS only interrogates gene-dense testable clq0.9 blocks, so the background is NOT
     all ~27k genes -- it is the genes overlapping any SCANNED block (+/-FLANK). We also run the
     sharper contrast with the SNP-HIT genes as the universe (does the kMate-unique layer differ
     functionally from what SNPs catch?).
  2. MULTIPLE TESTING over GO terms: BH-FDR across all tested terms.

Reuses nonsnp_only_genes.py's EXACT significance + block->gene logic (imported), so the foreground
sets match nonsnp_only_genes.csv. GO annotations pulled from the public GO Consortium GAF and the
go-basic OBO (propagated up is_a/part_of), cached under CACHE. Rolls its own hypergeometric
(scipy) -- no goatools/gseapy needed. Env: kmate (needs internet; compute node reaches
current.geneontology.org + purl.obolibrary.org).

Writes analysis/grenenet_gea/varexp/go_enrichment_{ora,themes}.csv + go_enrichment_summary.json.
"""
from __future__ import annotations
import os, sys, re, gzip, json
import numpy as np, pandas as pd
import requests
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
import nonsnp_only_genes as N   # reuse load_class/sig_blocks/climate_p/block_interval_map/FLANK

OUT = f"{lib.GEA}/varexp"
CACHE = f"{lib.GEA}/varexp/go_cache"
GAF_URL = "http://current.geneontology.org/annotations/tair.gaf.gz"
OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"
AGI_RE = re.compile(r"AT[1-5CM]G\d{5}", re.I)
MIN_TERM, MAX_FRAC = 5, 0.5     # ORA: term must have >=5 genes in universe and <50% of it

# ---- theme term-name regexes (matched against GO term NAMES, then genes via propagated annot) ----
THEMES = {
    "cold/freezing":        r"\bcold\b|freezing|chilling|cold acclimation",
    "heat":                 r"response to heat|heat acclimation|thermotoler|high light",
    "ABA/water/drought":    r"abscisic|water deprivation|drought|osmotic",
    "salt/ionic":           r"\bsalt\b|salt stress|sodium|ionic",
    "oxidative":            r"oxidative stress|reactive oxygen|hydrogen peroxide",
    "flowering/photoperiod":r"flower development|flowering|photoperiodism|vernaliz|floral",
    "circadian/clock":      r"circadian|rhythmic",
    "defense/biotic":       r"defense response|response to biotic|immune|to bacterium|to fungus",
    "abiotic (broad)":      r"response to abiotic stimulus",
    "stress (broad)":       r"response to stress",
    "temperature (broad)":  r"response to temperature stimulus",
}


# --------------------------------------------------------------------- GO data
def _fetch(url, path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"  downloading {url}", flush=True)
        hdr = {"User-Agent": "Mozilla/5.0 (kMate GEA enrichment; contact tatiana.bellagio@hhmi.org)"}
        with requests.get(url, headers=hdr, stream=True, timeout=120, allow_redirects=True) as r:
            r.raise_for_status()
            with open(path, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
    return path


def load_go():
    os.makedirs(CACHE, exist_ok=True)
    gaf = _fetch(GAF_URL, f"{CACHE}/tair.gaf.gz")
    obo = _fetch(OBO_URL, f"{CACHE}/go-basic.obo")

    # obo: id -> (name, namespace, parents[is_a/part_of])
    name, ns, parents = {}, {}, {}
    cur = None
    with open(obo) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = {}
            elif line == "" and cur is not None and "id" in cur:
                gid = cur["id"]
                name[gid] = cur.get("name", "")
                ns[gid] = cur.get("namespace", "")
                parents[gid] = cur.get("parents", [])
                cur = None
            elif cur is not None and ":" in line:
                k, _, v = line.partition(": ")
                if k == "id" and v.startswith("GO:"):
                    cur["id"] = v
                elif k == "name":
                    cur["name"] = v
                elif k == "namespace":
                    cur["namespace"] = v
                elif k == "is_a":
                    cur.setdefault("parents", []).append(v.split(" ! ")[0].strip())
                elif k == "relationship" and v.startswith("part_of "):
                    cur.setdefault("parents", []).append(v.split()[1])

    # ancestor closure (memoized DFS)
    anc_cache = {}
    def ancestors(gid):
        if gid in anc_cache:
            return anc_cache[gid]
        acc = set()
        for p in parents.get(gid, []):
            if p in name:
                acc.add(p); acc |= ancestors(p)
        anc_cache[gid] = acc
        return acc

    # gaf: AGI -> direct GO (skip NOT-qualified); then propagate
    direct = {}
    opener = gzip.open if gaf.endswith(".gz") else open
    with opener(gaf, "rt") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 15:
                continue
            qual, go_id, aspect, syn = f[3], f[4], f[8], f[10]
            if "NOT" in qual or not go_id.startswith("GO:"):
                continue
            m = AGI_RE.findall(" ".join([f[2], f[9], syn]))
            if not m:
                continue
            agi = m[0].upper()
            direct.setdefault(agi, set()).add(go_id)

    gene2go = {}
    for agi, gos in direct.items():
        allg = set()
        for g in gos:
            if g in name:
                allg.add(g); allg |= ancestors(g)
        gene2go[agi] = allg
    return gene2go, name, ns


# ------------------------------------------------------- sig genes + universe
def blocks_to_genes(blockset, bmap, genes, flank):
    ids = set()
    for b in blockset:
        if b not in bmap:
            continue
        c, s, e = bmap[b]
        inb = genes[(genes.chrom == c) & (genes.start <= e) & (genes.end >= s)]
        flk = genes[(genes.chrom == c) & (genes.start <= e + flank) & (genes.end >= s - flank)]
        ids.update(inb.gene); ids.update(flk.gene)
    return ids


def build_sets():
    """Returns dict of gene-id sets: nonsnp_only_{fdr,bonf}, snp_hit_{fdr,bonf}, universe."""
    data = {c: N.load_class(c) for c in N.CLASSES}
    sites = data["snp"]["sites"]
    blk = {c: N.block_ids(data[c]["chrom"], data[c]["pos"]) for c in N.CLASSES}
    Cmat = {c: np.corrcoef(data[c]["Z"].T) for c in N.CLASSES}
    clim = lib.load_climate().reindex(sites)[lib.BIO_COLS]
    Bz = (clim - clim.mean()) / clim.std()
    Uu, Ss, _ = np.linalg.svd(Bz.to_numpy(), full_matrices=False)
    axes = {f"bio{i}": clim[f"bio{i}"].to_numpy(float) for i in range(1, 20)}
    axes["PC1_allbio"] = Uu[:, 0] * Ss[0]

    contrasts = {}
    for nm, key in [("JOINT", "p_joint"), ("GLOBAL", "p_global")]:
        contrasts[f"multitrait_{nm}"] = {c: data[c][key] for c in N.CLASSES}
    for ax, bvec in axes.items():
        contrasts[f"multitrait_CLIMATE_{ax}"] = {c: N.climate_p(data[c]["Z"], Cmat[c], bvec) for c in N.CLASSES}
    for si, sid in enumerate(sites):
        contrasts[f"persite_site{int(sid)}"] = {c: 2 * stats.norm.sf(np.abs(data[c]["Z"][:, si])) for c in N.CLASSES}

    nsonly_fdr, nsonly_bonf = set(), set()   # blocks non-SNP-only
    snp_fdr_all, snp_bonf_all = set(), set()  # blocks SNP-significant (any contrast)
    for pv in contrasts.values():
        sf, sb = N.sig_blocks(blk["snp"], pv["snp"])
        nf, nb = N.sig_blocks(blk["nonsnp"], pv["nonsnp"])
        snp_fdr_all |= sf; snp_bonf_all |= sb
        nsonly_fdr |= (nf - sf); nsonly_bonf |= (nb - sb)

    bmap = N.block_interval_map()
    genes = lib.load_genes()
    universe_blocks = set(b for b in np.unique(np.concatenate([blk["snp"], blk["nonsnp"]])) if b)
    S = dict(
        nonsnp_only_fdr = blocks_to_genes(nsonly_fdr, bmap, genes, N.FLANK),
        nonsnp_only_bonf= blocks_to_genes(nsonly_bonf, bmap, genes, N.FLANK),
        snp_hit_fdr     = blocks_to_genes(snp_fdr_all, bmap, genes, N.FLANK),
        snp_hit_bonf    = blocks_to_genes(snp_bonf_all, bmap, genes, N.FLANK),
        universe        = blocks_to_genes(universe_blocks, bmap, genes, N.FLANK),
    )
    return S


# --------------------------------------------------------------- enrichment
def hyper(fg, universe, term_genes):
    """one-sided over-representation: P(X>=k) with X~Hypergeom(N,K,n)."""
    N_ = len(universe); n = len(fg)
    K = len(term_genes & universe)
    k = len(fg & term_genes)
    if k == 0 or K == 0:
        return k, K, n, N_, 1.0
    p = stats.hypergeom.sf(k - 1, N_, K, n)
    return k, K, n, N_, float(p)


def ora(fg, universe, gene2go, go_name, go_ns, label):
    fg = fg & universe                         # restrict fg to universe
    term2genes = {}
    for g in universe:
        for t in gene2go.get(g, ()):
            term2genes.setdefault(t, set()).add(g)
    rows = []
    for t, tg in term2genes.items():
        if go_ns.get(t) != "biological_process":
            continue
        if not (MIN_TERM <= len(tg) <= MAX_FRAC * len(universe)):
            continue
        k, K, n, Nn, p = hyper(fg, universe, tg)
        if k >= 2:
            exp = n * K / Nn
            rows.append(dict(GO=t, name=go_name.get(t, ""), k=k, K=K, n=n, N=Nn,
                             expected=round(exp, 2), fold=round(k / exp, 2) if exp else np.nan, p=p))
    df = pd.DataFrame(rows)
    if len(df):
        df["q_BH"] = lib.bh(df["p"].to_numpy())
        df = df.sort_values("p").reset_index(drop=True)
        df.insert(0, "test", label)
    return df


def theme_test(fg, universe, gene2go, go_name, go_ns, label):
    fg = fg & universe
    # theme -> set of BP terms whose name matches, then genes annotated to any of them
    rows = []
    for theme, rx in THEMES.items():
        terms = {t for t, nm in go_name.items() if go_ns.get(t) == "biological_process"
                 and re.search(rx, nm, re.I)}
        tg = {g for g in universe if gene2go.get(g, set()) & terms}
        k, K, n, Nn, p = hyper(fg, universe, tg)
        exp = n * K / Nn if Nn else np.nan
        rows.append(dict(test=label, theme=theme, n_terms=len(terms), k=k, K=K, n=n, N=Nn,
                         expected=round(exp, 2), fold=round(k / exp, 2) if exp else np.nan, p=p))
    df = pd.DataFrame(rows)
    df["q_BH"] = lib.bh(df["p"].to_numpy())
    return df.sort_values("p").reset_index(drop=True)


def main():
    print("[1/3] gene sets ...", flush=True)
    S = build_sets()
    for k, v in S.items():
        print(f"    {k:18s} {len(v)} genes", flush=True)
    print("[2/3] GO annotations ...", flush=True)
    gene2go, go_name, go_ns = load_go()
    annot_in_univ = sum(1 for g in S["universe"] if g in gene2go)
    print(f"    {len(gene2go)} genes annotated; {annot_in_univ}/{len(S['universe'])} of universe has GO", flush=True)

    print("[3/3] enrichment ...", flush=True)
    jobs = [   # (foreground, universe, label)
        (S["nonsnp_only_fdr"],  S["universe"],    "nsFDR_vs_scanned"),
        (S["nonsnp_only_bonf"], S["universe"],    "nsBONF_vs_scanned"),
        (S["nonsnp_only_fdr"],  S["snp_hit_fdr"], "nsFDR_vs_snpFDR"),
        (S["nonsnp_only_bonf"], S["snp_hit_bonf"],"nsBONF_vs_snpBONF"),
    ]
    ora_all, theme_all = [], []
    for fg, uni, lab in jobs:
        uni = uni | fg                         # ensure fg subset of universe
        ora_all.append(ora(fg, uni, gene2go, go_name, go_ns, lab))
        theme_all.append(theme_test(fg, uni, gene2go, go_name, go_ns, lab))
    ora_df = pd.concat([d for d in ora_all if len(d)], ignore_index=True)
    theme_df = pd.concat(theme_all, ignore_index=True)
    ora_df.to_csv(f"{OUT}/go_enrichment_ora.csv", index=False)
    theme_df.to_csv(f"{OUT}/go_enrichment_themes.csv", index=False)

    print("\n===== THEME TESTS (targeted) =====")
    for lab in [j[2] for j in jobs]:
        sub = theme_df[theme_df.test == lab]
        print(f"\n-- {lab} (fg n={sub.n.iloc[0]}, universe N={sub.N.iloc[0]}) --")
        print(sub[["theme", "k", "K", "expected", "fold", "p", "q_BH"]].to_string(index=False))
    print("\n===== TOP UNBIASED GO-BP TERMS (q_BH<0.1) =====")
    for lab in [j[2] for j in jobs]:
        sub = ora_df[(ora_df.test == lab) & (ora_df.q_BH < 0.1)].head(15)
        print(f"\n-- {lab}: {len(ora_df[(ora_df.test==lab)&(ora_df.q_BH<0.1)])} terms at q<0.1 --")
        if len(sub):
            print(sub[["GO", "name", "k", "K", "fold", "p", "q_BH"]].to_string(index=False))

    summary = {lab: {
        "fg_n": int(theme_df[theme_df.test == lab].n.iloc[0]),
        "universe_N": int(theme_df[theme_df.test == lab].N.iloc[0]),
        "themes_sig_q05": theme_df[(theme_df.test == lab) & (theme_df.q_BH < 0.05)][["theme", "fold", "q_BH"]].to_dict("records"),
        "n_ora_terms_q10": int(len(ora_df[(ora_df.test == lab) & (ora_df.q_BH < 0.1)])),
    } for lab in [j[2] for j in jobs]}
    json.dump(summary, open(f"{OUT}/go_enrichment_summary.json", "w"), indent=2, default=str)
    print(f"\nwrote go_enrichment_{{ora,themes}}.csv + go_enrichment_summary.json")


if __name__ == "__main__":
    main()

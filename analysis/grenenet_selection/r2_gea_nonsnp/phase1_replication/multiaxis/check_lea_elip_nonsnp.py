#!/usr/bin/env python
"""Presence check: LEA (Late Embryogenesis Abundant / dehydrins) and ELIP (Early
Light-Inducible Proteins + LIL/OHP relatives) protein families among the *non-SNP*
(SV + smallindel = pooled non-SNP) LFMM new-peak genes, clq0.9 tiling, all 20 axes.

Source gene lists (results/multiaxis/):
  raw_manhattan_newpeak_genes_nonsnp_tile.csv  -- pooled non-SNP ("all non-SNP")
  raw_manhattan_newpeak_genes_sv_tile.csv      -- SV subset

Family membership is decided two ways and unioned:
  (1) keyword scan of the project annotator's text (TAIR GO + UniProt protein name /
      FUNCTION / keywords) -- the authoritative, annotation-driven route (CLAUDE.md);
  (2) a small curated AGI-ID set of canonical members (Hundertmark & Hincha 2008 dehydrins;
      ELIP1/2 + OHP1/2 + LIL3) as an independent cross-check so nothing is missed if an
      entry lacks a reviewed UniProt name.
"""
import importlib.util as ilu
import os, re, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
import lib  # noqa: E402
RES = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis"

_sp = ilu.spec_from_file_location("ann_tu", f"{os.path.dirname(HERE)}/annotate_genes_tair_uniprot.py")
ann_tu = ilu.module_from_spec(_sp); _sp.loader.exec_module(ann_tu)

# --- curated canonical AGI IDs (independent cross-check) ---------------------
LEA_CURATED = {  # dehydrins (LEA group 2) + a few well-known other-group LEAs
    "AT1G20440": "COR47", "AT1G20450": "ERD10/LTI45", "AT1G76180": "ERD14",
    "AT5G66400": "RAB18", "AT3G50970": "LTI30/XERO2", "AT3G50980": "XERO1",
    "AT1G54410": "dehydrin family", "AT2G42540": "COR15A", "AT2G42530": "COR15B",
    "AT1G01470": "LEA14", "AT3G15670": "LEA (group1)", "AT2G36640": "ECP63/LEA",
    "AT3G17520": "LEA-D29", "AT1G52690": "LEA7", "AT4G36600": "LEA",
    "AT5G06760": "LEA4-5", "AT1G02820": "LEA3-family",
}
ELIP_CURATED = {
    "AT3G22840": "ELIP1", "AT4G14690": "ELIP2",
    "AT5G02120": "OHP1/LIL2", "AT1G34000": "OHP2",
    "AT4G17600": "LIL3:1", "AT5G47110": "LIL3:2",
}

LEA_RE = re.compile(
    r"late embryogenesis abundant|dehydrin|hydrophilin|\bLEA\b|LEA[_ -]?\d|"
    r"\bERD1[04]\b|\bCOR47\b|\bRAB18\b|\bXERO\d\b|\bLTI(?:29|30|45)\b|seed maturation protein",
    re.I)
ELIP_RE = re.compile(
    r"early light[- ]?induc|\bELIP\d?\b|one[- ]?helix protein|\bOHP\d?\b|"
    r"light[- ]?harvesting[- ]?like|\bLIL\d",
    re.I)
# weaker, reported separately (bulk LHCs are ELIP relatives but not the photoprotective set)
CAB_RE = re.compile(r"chlorophyll a[/-]?b[- ]?binding|light-harvesting chlorophyll", re.I)


def load_lists():
    sv = pd.read_csv(f"{RES}/raw_manhattan_newpeak_genes_sv_tile.csv")
    ns = pd.read_csv(f"{RES}/raw_manhattan_newpeak_genes_nonsnp_tile.csv")
    sv["lists"] = "sv"; ns["lists"] = "nonsnp"
    both = pd.concat([sv, ns], ignore_index=True)
    # per gene: which lists, best (max nlp) axis, n_axes max, span
    g = (both.sort_values("nlp", ascending=False)
         .groupby("gene")
         .agg(lists=("lists", lambda s: "+".join(sorted(set(s)))),
              best_axis=("axis", "first"), best_nlp=("nlp", "first"),
              n_axes=("n_axes", "max"), chrom=("chrom", "first"), pos=("pos", "first"))
         .reset_index())
    return g


def main():
    g = load_lists()
    genes = sorted(g.gene.unique())
    print(f"non-SNP new-peak genes (sv + pooled non-SNP union): {len(genes)}", flush=True)

    A = ann_tu.annotate(genes)
    A = A.set_index("gene")
    blob = (A.get("protein_name", "").fillna("") + " || " + A.get("tair_name", "").fillna("")
            + " || " + A.get("mygene_name", "").fillna("") + " || "
            + A.get("uniprot_function", "").fillna("") + " || "
            + A.get("uniprot_keywords", "").fillna(""))

    rows = []
    for gene in genes:
        txt = blob.get(gene, "")
        agi = gene.upper()
        lea = bool(LEA_RE.search(txt)) or agi in LEA_CURATED
        elip = bool(ELIP_RE.search(txt)) or agi in ELIP_CURATED
        cab = bool(CAB_RE.search(txt)) and not elip
        if lea or elip or cab:
            fam = []
            if lea: fam.append("LEA/dehydrin")
            if elip: fam.append("ELIP/LIL/OHP")
            if cab: fam.append("Chl-a/b (LHC, weak ELIP relative)")
            rows.append(dict(
                gene=gene, family=";".join(fam),
                curated=LEA_CURATED.get(agi, "") or ELIP_CURATED.get(agi, ""),
                symbol=A.get("symbol", pd.Series()).get(gene, ""),
                protein_name=A.get("protein_name", pd.Series()).get(gene, ""),
                keywords=A.get("uniprot_keywords", pd.Series()).get(gene, "")))
    hits = pd.DataFrame(rows)
    if len(hits):
        hits = hits.merge(g, on="gene", how="left").sort_values(["family", "best_nlp"],
                                                                ascending=[True, False])
    out = f"{RES}/lea_elip_nonsnp_newpeak_hits.csv"
    hits.to_csv(out, index=False)

    print("\n" + "=" * 78)
    for fam, lab in [("LEA/dehydrin", "LEA (Late Embryogenesis Abundant / dehydrins)"),
                     ("ELIP/LIL/OHP", "ELIP (Early Light-Inducible + LIL/OHP)"),
                     ("Chl-a/b", "Chl a/b-binding (LHC; weak ELIP relatives)")]:
        sub = hits[hits.family.str.contains(fam)] if len(hits) else hits
        print(f"\n### {lab}: {len(sub)} gene(s)")
        for _, r in sub.iterrows():
            print(f"  {r.gene}  {r.symbol or r.curated or '':14s}  lists={r.lists:9s} "
                  f"axes={int(r.n_axes):2d} bestp(-log10)={r.best_nlp:5.1f}@{r.best_axis:5s} "
                  f"{r.chrom}:{int(r.pos):>9d}  | {str(r.protein_name)[:60]}")
    print("\n" + "=" * 78)
    print(f"wrote {out}  ({len(hits)} family hits)")


if __name__ == "__main__":
    main()

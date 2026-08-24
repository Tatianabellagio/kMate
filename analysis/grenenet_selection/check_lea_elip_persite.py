#!/usr/bin/env python
"""Presence check: LEA (Late Embryogenesis Abundant / dehydrins) and ELIP (Early
Light-Inducible Proteins + LIL/OHP relatives) among the PER-SITE (single-garden) new-peak
genes, clq0.9 TILING partition.

Reads the per-site gene tables written by _build_persite_new_peaks_nb.py:
  varexp/persite_new_peaks_nonsnp_genes.csv
  varexp/persite_new_peaks_sv_genes.csv
Family membership: curated canonical AGI IDs UNION annotator-text regex (TAIR GO + UniProt
protein name / FUNCTION / keywords). Same curated sets + regexes as
phase1_replication/multiaxis/check_lea_elip_nonsnp.py.
"""
import importlib.util as ilu
import os, re, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib  # noqa: E402
DIR = f"{lib.GEA}/r3_persite_gwas/results/varexp"

_sp = ilu.spec_from_file_location(
    "ann_tu", f"{lib.GEA}/phase1_replication/annotate_genes_tair_uniprot.py")
ann_tu = ilu.module_from_spec(_sp); _sp.loader.exec_module(ann_tu)

LEA_CURATED = {
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
CAB_RE = re.compile(r"chlorophyll a[/-]?b[- ]?binding|light-harvesting chlorophyll", re.I)


def main():
    frames = {}
    for cls in ["nonsnp", "sv"]:
        f = f"{DIR}/persite_new_peaks_{cls}_genes.csv"
        if os.path.exists(f):
            d = pd.read_csv(f); d["cls"] = cls; frames[cls] = d
    if not frames:
        sys.exit("no per-site gene CSVs found -- run _build_persite_new_peaks_nb.py first")
    allg = pd.concat(frames.values(), ignore_index=True)
    # per gene: which classes, max recurrence, max non-SNP-only recurrence, strongest nlp
    g = (allg.sort_values("nlp", ascending=False).groupby("gene")
         .agg(classes=("cls", lambda s: "+".join(sorted(set(s)))),
              n_sites=("n_sites", "max"), n_sites_nonsnp_only=("n_sites_nonsnp_only", "max"),
              best_nlp=("nlp", "first"), chrom=("chrom", "first"), pos=("pos", "first"))
         .reset_index())
    genes = sorted(g.gene.unique())
    print(f"per-site new-peak genes (nonsnp + sv union): {len(genes)}", flush=True)

    A = ann_tu.annotate(genes).set_index("gene")
    blob = (A.get("protein_name", "").fillna("") + " || " + A.get("tair_name", "").fillna("")
            + " || " + A.get("mygene_name", "").fillna("") + " || "
            + A.get("uniprot_function", "").fillna("") + " || "
            + A.get("uniprot_keywords", "").fillna(""))

    rows = []
    for gene in genes:
        txt = blob.get(gene, ""); agi = gene.upper()
        lea = bool(LEA_RE.search(txt)) or agi in LEA_CURATED
        elip = bool(ELIP_RE.search(txt)) or agi in ELIP_CURATED
        cab = bool(CAB_RE.search(txt)) and not elip
        if lea or elip or cab:
            fam = []
            if lea: fam.append("LEA/dehydrin")
            if elip: fam.append("ELIP/LIL/OHP")
            if cab: fam.append("Chl-a/b (LHC, weak ELIP relative)")
            rows.append(dict(gene=gene, family=";".join(fam),
                             curated=LEA_CURATED.get(agi, "") or ELIP_CURATED.get(agi, ""),
                             symbol=A.get("symbol", pd.Series()).get(gene, ""),
                             protein_name=A.get("protein_name", pd.Series()).get(gene, ""),
                             keywords=A.get("uniprot_keywords", pd.Series()).get(gene, "")))
    hits = pd.DataFrame(rows)
    if len(hits):
        hits = hits.merge(g, on="gene", how="left").sort_values(
            ["family", "n_sites"], ascending=[True, False])
    out = f"{DIR}/lea_elip_persite_hits.csv"
    hits.to_csv(out, index=False)

    print("\n" + "=" * 84)
    for fam, lab in [("LEA/dehydrin", "LEA (Late Embryogenesis Abundant / dehydrins)"),
                     ("ELIP/LIL/OHP", "ELIP (Early Light-Inducible + LIL/OHP)"),
                     ("Chl-a/b", "Chl a/b-binding (LHC; weak ELIP relatives)")]:
        sub = hits[hits.family.str.contains(fam)] if len(hits) else hits
        print(f"\n### {lab}: {len(sub)} gene(s)")
        for _, r in sub.iterrows():
            print(f"  {r.gene}  {str(r.symbol or r.curated or ''):14s} classes={r.classes:9s} "
                  f"n_sites={int(r.n_sites):2d} (nonSNP-only {int(r.n_sites_nonsnp_only):2d}) "
                  f"bestp={r.best_nlp:4.1f} {r.chrom}:{int(r.pos):>9d} | {str(r.protein_name)[:52]}")
    print("\n" + "=" * 84)
    print(f"wrote {out}  ({len(hits)} family hits)")


if __name__ == "__main__":
    main()

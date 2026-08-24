#!/usr/bin/env python
"""Arabidopsis-specific gene annotation: TAIR GO (GO Consortium GAF) + UniProt function.

WHY THIS EXISTS -- the gap in annotate_gene_function.py
------------------------------------------------------
`annotate_gene_function.py` pulls name / summary / GO from mygene.info (NCBI Entrez).
Measured on the 73 genes of the non-SNP WZA hit list:

    name                73/73
    GO BP terms         62/73
    NCBI summary         0/73   <-- NCBI has essentially no free-text summaries for
    categories assigned 25/73        Arabidopsis loci

Because `classify()` keyword-matches against name + summary + GO, an empty summary column
means classification runs on name + GO only, and for the 11 genes with neither it runs on
the gene NAME STRING alone. So the category tags were an undercount of unknown size.

Two Arabidopsis-specific sources fix this, both public and unauthenticated:

  1. TAIR GO annotation -- the GO Consortium's canonical TAIR GAF
     (http://current.geneontology.org/annotations/tair.gaf.gz). Complete curated GO for
     every TAIR locus; strictly better coverage than mygene's subset. Cached locally.
     NB TAIR's own download URLs (arabidopsis.org/download_files/...) now return HTTP 403,
     so the GO Consortium mirror is the usable route.

  2. UniProt REST -- reviewed (Swiss-Prot) entries for taxid 3702, giving
     `protein_name` (far more informative than the GFF `Name=`: e.g. "26.5 kDa heat shock
     protein, mitochondrial (AtHsp26.5)" vs "HSP20-like chaperones superfamily protein"),
     `cc_function` (curated free-text FUNCTION), and `keyword` (a curated controlled
     vocabulary that includes exactly the terms we classify on -- "Stress response",
     "Flowering", "Abscisic acid signaling pathway", ...).
     AGI locus IDs come back in the `gene_oln` (ordered locus names) field, so results map
     back unambiguously.

Ensembl was considered and rejected: Ensembl Plants serves the same TAIR10-assembly /
Araport11 gene models, and its `description` is the same string as the TAIR10 GFF `Name=`.
It adds coordinates we already have locally and no functional text we lack.

Classification reuses annotate_gene_function.py's curated CATEGORIES so the output column
is comparable to the existing tables, but now matched against
  name + protein_name + uniprot_function + uniprot_keywords + GO(TAIR) + GO(mygene)

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python   # needs outbound HTTPS
  $PY annotate_genes_tair_uniprot.py --genes-from <csv with a `gene` column> --out <csv>
  $PY annotate_genes_tair_uniprot.py --genes AT1G52560 AT1G13440
"""
from __future__ import annotations
import argparse, gzip, importlib.util, os, sys, time, urllib.parse, urllib.request
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import lib  # noqa: E402

GAF_URL = "http://current.geneontology.org/annotations/tair.gaf.gz"
GAF_CACHE = f"{lib.GEA}/phase1_replication/results/cache/tair.gaf.gz"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"

# reuse the curated keyword sets + classifier so the `categories` column stays comparable
_spec = importlib.util.spec_from_file_location("agf", f"{HERE}/annotate_gene_function.py")
agf = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(agf)


UA = "Mozilla/5.0 (kmate GEA annotation; contact tbellg)"


def fetch_gaf(force=False) -> str:
    os.makedirs(os.path.dirname(GAF_CACHE), exist_ok=True)
    if force or not os.path.exists(GAF_CACHE) or os.path.getsize(GAF_CACHE) < 1_000_000:
        print(f"downloading {GAF_URL} ...", flush=True)
        # NB current.geneontology.org 403s on urllib's default User-Agent; curl works.
        # Send an explicit UA rather than urlretrieve (which cannot set headers).
        req = urllib.request.Request(GAF_URL, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=180) as r, open(GAF_CACHE, "wb") as fh:
            while True:
                buf = r.read(1 << 20)
                if not buf:
                    break
                fh.write(buf)
    print(f"GAF cache: {GAF_CACHE} ({os.path.getsize(GAF_CACHE)/1e6:.1f} MB)")
    return GAF_CACHE


def tair_go(genes: set[str]) -> dict[str, dict]:
    """AGI locus -> {go_bp, go_mf, go_cc} from the TAIR GAF.

    GAF columns (0-based): 1 DB_Object_ID, 2 DB_Object_Symbol, 4 GO_ID, 8 Aspect(P/F/C),
    9 DB_Object_Name, 10 Synonyms. TAIR puts the AGI locus in the synonyms field, so match
    on both symbol and synonyms, upper-cased.
    """
    want = {g.upper() for g in genes}
    out: dict[str, dict] = {}
    # GO_ID -> term name is not in the GAF; keep the GO id plus the object name, and use
    # the aspect to split BP/MF/CC. Term NAMES come from mygene where available; for
    # classification the DB_Object_Name + UniProt text carry the signal.
    with gzip.open(fetch_gaf(), "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            keys = {f[2].upper()} | {s.strip().upper() for s in f[10].split("|") if s.strip()}
            hit = keys & want
            if not hit:
                continue
            aspect = {"P": "go_bp_ids", "F": "go_mf_ids", "C": "go_cc_ids"}.get(f[8])
            for g in hit:
                r = out.setdefault(g, dict(go_bp_ids=set(), go_mf_ids=set(),
                                           go_cc_ids=set(), tair_name=""))
                if aspect:
                    r[aspect].add(f[4])
                if f[9] and not r["tair_name"]:
                    r["tair_name"] = f[9]
    for g, r in out.items():
        for k in ("go_bp_ids", "go_mf_ids", "go_cc_ids"):
            r[k] = ";".join(sorted(r[k]))
    return out


def uniprot_batch(genes: list[str], chunk=40, retries=3, reviewed=True) -> dict[str, dict]:
    """AGI locus -> UniProt protein name / FUNCTION / keywords.

    `reviewed=True` restricts to Swiss-Prot (curated). Many Arabidopsis loci have no
    reviewed entry, so `annotate()` runs a second pass with reviewed=False over whatever
    is still unmapped and records which tier each gene came from in `uniprot_reviewed`.
    """
    out: dict[str, dict] = {}
    for i in range(0, len(genes), chunk):
        sub = genes[i:i + chunk]
        q = ("(" + " OR ".join(f"gene:{g}" for g in sub) + ")"
             " AND (organism_id:3702)" + (" AND (reviewed:true)" if reviewed else ""))
        url = UNIPROT + "?" + urllib.parse.urlencode({
            "query": q,
            "fields": "accession,gene_oln,gene_primary,protein_name,cc_function,keyword",
            "format": "tsv", "size": "500"})
        for a in range(retries):
            try:
                rq = urllib.request.Request(url, headers={"User-Agent": UA})
                with urllib.request.urlopen(rq, timeout=90) as r:
                    txt = r.read().decode("utf-8", "replace")
                rows = [l.split("\t") for l in txt.strip().split("\n")[1:] if l.strip()]
                for f in rows:
                    if len(f) < 6:
                        continue
                    acc, oln, prim, pname, fn, kw = f[:6]
                    # gene_oln holds the AGI locus (e.g. "At1g52560"), possibly several
                    for tok in [t.strip().upper() for t in oln.replace(",", " ").split()]:
                        if tok in {g.upper() for g in sub}:
                            out[tok] = dict(
                                uniprot=acc, uniprot_symbol=prim, protein_name=pname,
                                uniprot_function=fn.replace("FUNCTION: ", "").strip(),
                                uniprot_keywords=kw, uniprot_reviewed=reviewed)
                break
            except Exception as e:
                print(f"  [uniprot] chunk {i} attempt {a+1} failed: {e}", flush=True)
                time.sleep(3)
        print(f"  uniprot {min(i+chunk, len(genes))}/{len(genes)} -> {len(out)} mapped",
              flush=True)
        time.sleep(0.3)
    return out


def annotate(genes: list[str]) -> pd.DataFrame:
    genes = sorted({g.strip() for g in genes if isinstance(g, str) and g.strip()})
    print(f"annotating {len(genes)} genes\n")
    print("-- TAIR GAF (GO Consortium) --", flush=True)
    go = tair_go(set(genes))
    print(f"   GO for {len(go)}/{len(genes)} genes")
    print("-- UniProt REST pass 1: reviewed / Swiss-Prot (taxid 3702) --", flush=True)
    up = uniprot_batch(genes, reviewed=True)
    print(f"   reviewed: {len(up)}/{len(genes)} genes")
    missing = [g for g in genes if g.upper() not in up]
    if missing:
        print(f"-- UniProt REST pass 2: unreviewed / TrEMBL for the remaining {len(missing)} --",
              flush=True)
        up.update(uniprot_batch(missing, reviewed=False))
    print(f"   UniProt total: {len(up)}/{len(genes)} genes")
    print("-- mygene.info (kept for GO term NAMES + entrez) --", flush=True)
    mg = agf.mygene_batch(genes)

    rows = []
    for g in genes:
        G, U, M = go.get(g.upper(), {}), up.get(g.upper(), {}), (mg.get(g) or {})
        rec = dict(
            gene=g,
            symbol=U.get("uniprot_symbol") or M.get("symbol", "") or "",
            protein_name=U.get("protein_name", ""),
            tair_name=G.get("tair_name", ""),
            mygene_name=M.get("name", ""),
            uniprot=U.get("uniprot", ""),
            uniprot_function=U.get("uniprot_function", ""),
            uniprot_keywords=U.get("uniprot_keywords", ""),
            uniprot_reviewed=U.get("uniprot_reviewed", ""),
            go_bp=";".join(M.get("go_bp") or []),
            go_mf=";".join(M.get("go_mf") or []),
            go_bp_ids=G.get("go_bp_ids", ""),
            n_go_tair=len([x for x in G.get("go_bp_ids", "").split(";") if x]),
        )
        # classify on EVERYTHING now available, not just name+summary+GO
        blob = dict(
            symbol=rec["symbol"], name=" ".join([rec["protein_name"], rec["tair_name"],
                                                 rec["mygene_name"]]),
            summary=" ".join([rec["uniprot_function"], rec["uniprot_keywords"]]),
            go_bp=(M.get("go_bp") or []), go_mf=(M.get("go_mf") or []))
        cats, ev = agf.classify(blob)
        rec["categories"] = ",".join(cats)
        rec["evidence"] = " | ".join(ev)[:400]
        rec["climate_stress_flowering"] = bool(cats)
        rows.append(rec)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", nargs="*", default=None)
    ap.add_argument("--genes-from", default=None,
                    help="CSV with a `gene` column (e.g. *_bonf_genes.csv)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    if a.genes_from:
        src = pd.read_csv(a.genes_from)
        genes = sorted(set(src["gene"].dropna()))
        out = a.out or a.genes_from.replace(".csv", "_tairuniprot.csv")
    elif a.genes:
        genes, out = a.genes, (a.out or "gene_annotation_tairuniprot.csv")
    else:
        sys.exit("give --genes or --genes-from")

    D = annotate(genes)
    D.to_csv(out, index=False)
    print(f"\nwrote {out}  ({len(D)} genes)\n")
    print("COVERAGE")
    for c, lab in [("protein_name", "UniProt protein name"),
                   ("uniprot_function", "UniProt FUNCTION text"),
                   ("uniprot_keywords", "UniProt keywords"),
                   ("go_bp_ids", "TAIR GO (BP)"),
                   ("go_bp", "mygene GO (BP)"),
                   ("uniprot_reviewed", "of which reviewed/unrev."),
                   ("categories", "categories assigned")]:
        print(f"  {lab:24s} {int((D[c].fillna('').astype(str).str.len() > 0).sum()):3d}/{len(D)}")


if __name__ == "__main__":
    main()

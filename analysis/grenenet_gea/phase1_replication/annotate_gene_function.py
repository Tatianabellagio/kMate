#!/usr/bin/env python
"""Functional annotation of the significant-block genes: flag climate / stress / flowering.

For every gene in significant_genes_*.csv, pull NCBI/Entrez name + summary + GO
biological-process & molecular-function terms from the mygene.info API (batch POST;
no auth; independent of TAIR10 GFF and Ensembl). Classify each gene into functional
categories with curated keyword sets matched against name+summary+GO, recording the
evidence that triggered each hit (GO hits are stronger than name-only).

Categories: flowering · temperature(cold/heat) · water(drought/osmotic/ABA/salt) ·
light/UV · oxidative/abiotic · defense(biotic) · calcium/signaling. The combined flag
`climate_stress_flowering` = any of {flowering, temperature, water, light, oxidative,
defense} (the environment/stress/phenology axes the GrENE-Net selection scan targets).

Outputs (analysis/grenenet_gea/phase1_replication/results/):
  gene_function_gen{g}_{clim}_{regime}.csv          per gene: symbol,entrez,name,summary,
        go_bp, go_mf, categories, evidence, climate_stress_flowering
  significant_genes_annotated_gen{g}_{clim}_{regime}.csv   the full (block,gene) sig table
        joined with the function columns (the deliverable table, now functionally tagged)

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python   # needs outbound HTTPS
  $PY analysis/grenenet_gea/phase1_replication/annotate_gene_function.py
"""
from __future__ import annotations
import argparse, json, os, sys, time, urllib.parse, urllib.request
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

MYGENE = "https://mygene.info/v3/query"

# curated keyword sets (lowercase); matched against name + summary + GO terms.
CATEGORIES = {
    "flowering":   ["flower", "floral", "flowering time", "vernaliz", "photoperiod",
                    "circadian", "inflorescence", "meristem identity", "reproductive develop",
                    "flc", "flowering locus", "transition to flowering"],
    "temperature": ["cold", "freezing", "chilling", "frost", "heat", "thermo",
                    "temperature", "heat shock", "cbf", "dreb", " cor ", "cold acclimation"],
    "water":       ["drought", "water deprivation", "dehydrat", "desiccation", "osmotic",
                    "salt stress", "salinity", "abscisic", " aba ", "stomat", "hyperosmotic"],
    "light":       ["response to light", "photomorphogenesis", "red light", "blue light",
                    "far-red", "ultraviolet", "response to uv", "uv-b", "shade avoidance",
                    "photoreceptor", "photosynthe"],
    "oxidative":   ["oxidative stress", "reactive oxygen", "abiotic stress", "hypoxia",
                    "wounding", "response to stress", "hydrogen peroxide"],
    "defense":     ["defense response", "immune", "pathogen", "bacteri", "fungal", "fungus",
                    "virus", "viral", "salicylic", "jasmon", "disease resistance",
                    "hypersensitive", "innate immun", "r protein", "effector"],
    "calcium":     ["calcium", "calmodulin", "calcium-dependent", "calcium ion",
                    "ca2+", "cbl", "cdpk"],
}
# climate/stress/flowering = environmental + stress + phenology axes (exclude pure calcium/signaling)
CLIMATE_SET = ["flowering", "temperature", "water", "light", "oxidative", "defense"]


def _terms(go_field):
    """mygene go.BP / go.MF -> list of term strings (handles dict or list-of-dicts)."""
    if not go_field:
        return []
    if isinstance(go_field, dict):
        go_field = [go_field]
    return [g.get("term", "") for g in go_field if isinstance(g, dict) and g.get("term")]


def mygene_batch(ids, chunk=100, retries=3):
    out = {}
    for i in range(0, len(ids), chunk):
        sub = ids[i:i + chunk]
        data = urllib.parse.urlencode({
            "q": ",".join(sub),
            "scopes": "symbol,alias,locus_tag,ensembl.gene,other_names",
            "species": "3702",
            "fields": "symbol,name,summary,go.BP.term,go.MF.term,entrezgene",
        }).encode()
        req = urllib.request.Request(MYGENE, data=data,
                                     headers={"Accept": "application/json"})
        for a in range(retries):
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    res = json.load(r)
                for h in res:
                    q = h.get("query")
                    if h.get("notfound") or q in out:
                        out.setdefault(q, None)
                        continue
                    out[q] = dict(
                        symbol=h.get("symbol", ""), name=h.get("name", ""),
                        summary=(h.get("summary", "") or ""),
                        entrez=h.get("entrezgene", ""),
                        go_bp=_terms(h.get("go", {}).get("BP")),
                        go_mf=_terms(h.get("go", {}).get("MF")))
                break
            except Exception as e:
                print(f"  [mygene] chunk {i} attempt {a+1} failed: {e}", flush=True)
                time.sleep(2)
        time.sleep(0.3)
    return out


def classify(info):
    """Return (categories list, evidence list) from name+summary+GO text."""
    if not info:
        return [], []
    bp = info["go_bp"]; mf = info["go_mf"]
    hay = " ".join([info["name"].lower(), info["summary"].lower()]
                   + [t.lower() for t in bp + mf])
    cats, ev = [], []
    for cat, kws in CATEGORIES.items():
        for kw in kws:
            if kw.strip() in hay:
                cats.append(cat)
                # tag whether the hit came from a GO term (stronger) or name/summary
                src = "GO" if any(kw.strip() in t.lower() for t in bp + mf) else "desc"
                ev.append(f"{cat}:{kw.strip()}({src})")
                break
    return cats, ev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regime", default="deg7cap2000")
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--gen", type=int, default=9)
    args = ap.parse_args()
    base = f"{lib.GEA}/phase1_replication/results"
    sig = pd.read_csv(f"{base}/significant_genes_gen{args.gen}_{args.climate}_{args.regime}.csv")
    genes = sorted(sig.loc[sig.gene.notna() & (sig.gene.astype(str) != ""), "gene"].unique())
    print(f"annotating {len(genes)} genes via mygene.info…", flush=True)
    info = mygene_batch(genes)
    got = sum(1 for g in genes if info.get(g))
    print(f"  resolved {got}/{len(genes)}", flush=True)

    rows = []
    for g in genes:
        ic = info.get(g)
        cats, ev = classify(ic)
        rows.append(dict(
            gene=g,
            symbol=(ic["symbol"] if ic else ""),
            entrez=(ic["entrez"] if ic else ""),
            fn_name=(ic["name"] if ic else ""),
            summary=((ic["summary"][:300]) if ic else ""),
            go_bp=(";".join(ic["go_bp"]) if ic else ""),
            go_mf=(";".join(ic["go_mf"]) if ic else ""),
            categories=";".join(sorted(set(cats))),
            evidence=";".join(ev),
            climate_stress_flowering=any(c in CLIMATE_SET for c in cats)))
    fn = pd.DataFrame(rows)
    fpath = f"{base}/gene_function_gen{args.gen}_{args.climate}_{args.regime}.csv"
    fn.to_csv(fpath, index=False)

    # join function onto the full (block,gene) significant table
    ann = sig.merge(fn.drop(columns=["symbol"]), on="gene", how="left")
    apath = f"{base}/significant_genes_annotated_gen{args.gen}_{args.climate}_{args.regime}.csv"
    ann.to_csv(apath, index=False)

    # summary
    nrel = int(fn.climate_stress_flowering.sum())
    print(f"\n=== function summary ===")
    print(f"genes annotated: {len(fn)} | climate/stress/flowering-related: {nrel}")
    from collections import Counter
    cc = Counter(c for cs in fn.categories for c in cs.split(";") if c)
    print("category counts:", dict(cc.most_common()))
    rel = fn[fn.climate_stress_flowering].sort_values("symbol")
    print(f"\nclimate/stress/flowering genes ({nrel}):")
    print(rel[["symbol", "gene", "categories", "fn_name"]].to_string(index=False))
    print(f"\n-> {fpath}\n-> {apath}")


if __name__ == "__main__":
    main()

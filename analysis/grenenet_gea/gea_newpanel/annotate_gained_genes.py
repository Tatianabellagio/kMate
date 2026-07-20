#!/usr/bin/env python
"""Fill symbol + Ensembl-Plants function for the non-SNP-only / SV-only 'gained' gene
tables. Reuses/updates the shared cache gene_descriptions.csv; fetches missing AT-IDs
from rest.ensembl.org (compute node has outbound access).

Usage: annotate_gained_genes.py FILE1.csv [FILE2.csv ...]
Each CSV must have a `gene` column; symbol/function are (re)written in place.
"""
import sys, os, json, time
import urllib.request as u
import pandas as pd

DESC = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/varexp/gene_descriptions.csv"


def load_cache():
    if os.path.exists(DESC):
        d = pd.read_csv(DESC).fillna("")
        return {r.gene: (r.symbol, r.function) for r in d.itertuples()}
    return {}


def fetch(ids, chunk=45):
    got = {}
    for i in range(0, len(ids), chunk):
        sub = ids[i:i + chunk]
        for att in range(4):
            try:
                req = u.Request("https://rest.ensembl.org/lookup/id",
                                data=json.dumps({"ids": sub}).encode(),
                                headers={"Content-Type": "application/json", "Accept": "application/json"})
                for g, v in json.load(u.urlopen(req, timeout=45)).items():
                    v = v or {}
                    got[g] = (v.get("display_name", "") or "",
                              (v.get("description", "") or "").split(" [Source")[0])
                break
            except Exception as e:
                time.sleep(2 * (att + 1))
                if att == 3:
                    print(f"[warn] ensembl chunk {i} failed: {type(e).__name__}")
    return got


def main(files):
    cache = load_cache()
    need = sorted({g for f in files for g in pd.read_csv(f)["gene"]})
    missing = [g for g in need if g not in cache or not cache[g][1]]
    if missing:
        print(f"fetching {len(missing)} genes from Ensembl…")
        cache.update(fetch(missing))
        pd.DataFrame([{"gene": g, "symbol": s, "function": fn} for g, (s, fn) in sorted(cache.items())]
                     ).to_csv(DESC, index=False)
    for f in files:
        df = pd.read_csv(f)
        df["symbol"] = df.gene.map(lambda g: cache.get(g, ("", ""))[0])
        df["function"] = df.gene.map(lambda g: cache.get(g, ("", ""))[1])
        df.to_csv(f, index=False)
        print(f"annotated {f}: {len(df)} rows, {df.function.astype(bool).sum()} with function")


if __name__ == "__main__":
    main(sys.argv[1:])

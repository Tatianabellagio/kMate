#!/usr/bin/env python
"""Describe the non-SNP-only candidate genes (nonsnp_only_genes.csv) via public gene APIs.

  - Ensembl Plants REST (batch POST /lookup/id): one-line gene description + biotype (TAIR/NCBI
    annotation). Primary description, robust batch call.
  - UniProt REST (per-gene search, reviewed-preferred): protein name + curated function comment
    (cc_function) where a Swiss-Prot/TrEMBL entry exists -- a fuller functional sentence than the
    Ensembl one-liner. Best-effort; genes with no UniProt hit just keep the Ensembl description.

Writes results/grenenet_gea/varexp/nonsnp_only_genes_described.csv (nonsnp_only_genes.csv + columns
ensembl_description, biotype, uniprot_protein, uniprot_function). Env: kmate (needs internet;
compute node reaches rest.ensembl.org + rest.uniprot.org).
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np, pandas as pd
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

OUT = f"{lib.GEA}/varexp"
ENS = "https://rest.ensembl.org/lookup/id"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


def ensembl_batch(ids, chunk=40):
    """POST /lookup/id in chunks -> {gene: (description, biotype, display_name/symbol)}."""
    out = {}
    for i in range(0, len(ids), chunk):
        sub = ids[i:i + chunk]
        for attempt in range(4):
            try:
                r = requests.post(ENS, headers={"Content-Type": "application/json",
                                                "Accept": "application/json"},
                                  data=json.dumps({"ids": sub}), timeout=40)
                r.raise_for_status()
                for gid, v in r.json().items():
                    out[gid] = ((v or {}).get("description", ""), (v or {}).get("biotype", ""),
                                (v or {}).get("display_name", "")) if v else ("", "", "")
                break
            except Exception as e:
                print(f"  ensembl chunk {i} attempt {attempt} failed: {e}", flush=True)
                time.sleep(2 * (attempt + 1))
        print(f"  ensembl: {min(i+chunk, len(ids))}/{len(ids)}", flush=True)
    return out


def uniprot_one(gid):
    """Reviewed-preferred UniProt entry for an AT locus -> (protein_name, function)."""
    for rev in ("+AND+reviewed:true", ""):
        try:
            url = (f"{UNIPROT}?query=gene:{gid}+AND+organism_id:3702{rev}"
                   f"&fields=protein_name,cc_function&format=tsv&size=1")
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            lines = r.text.strip().split("\n")
            if len(lines) >= 2:
                cols = lines[1].split("\t")
                prot = cols[0] if len(cols) > 0 else ""
                func = cols[1] if len(cols) > 1 else ""
                func = func.replace("FUNCTION: ", "").strip()
                if prot or func:
                    return prot, func
        except Exception as e:
            print(f"  uniprot {gid} ({'rev' if rev else 'all'}) failed: {e}", flush=True)
    return "", ""


def main():
    gdf = pd.read_csv(f"{OUT}/nonsnp_only_genes.csv")
    ids = gdf["gene"].tolist()
    print(f"{len(ids)} genes to describe", flush=True)

    ens = ensembl_batch(ids)
    gdf["symbol"] = gdf["gene"].map(lambda g: ens.get(g, ("", "", ""))[2])
    gdf["ensembl_description"] = gdf["gene"].map(lambda g: ens.get(g, ("", "", ""))[0])
    gdf["biotype"] = gdf["gene"].map(lambda g: ens.get(g, ("", "", ""))[1])

    prot, func = [], []
    for i, g in enumerate(ids):
        p, f = uniprot_one(g)
        prot.append(p); func.append(f)
        time.sleep(0.15)
        if (i + 1) % 20 == 0:
            print(f"  uniprot: {i+1}/{len(ids)}", flush=True)
    gdf["uniprot_protein"] = prot
    gdf["uniprot_function"] = func

    gdf.to_csv(f"{OUT}/nonsnp_only_genes_described.csv", index=False)
    n_ens = (gdf["ensembl_description"].fillna("") != "").sum()
    n_uni = (gdf["uniprot_function"].fillna("") != "").sum()
    print(f"\nwrote nonsnp_only_genes_described.csv: {len(gdf)} genes "
          f"({n_ens} with Ensembl desc, {n_uni} with UniProt function)", flush=True)
    with pd.option_context("display.max_colwidth", 70, "display.width", 200):
        print(gdf.head(15)[["gene", "symbol", "n_contrasts", "ensembl_description"]].to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Resolve AT-IDs -> gene symbols + description via Ensembl Plants REST (batch POST).

arabidopsis.org is 403-gated on the cluster; Ensembl Plants REST needs no auth and
carries TAIR10 symbols. Reads top_block_atids.txt from PB_DIR, writes atid_symbols.csv.
"""
import os, json, time, urllib.request

H = "results/grenenet_gea/hapfreq"
B = os.environ.get("PB_DIR", f"{H}/pipelineB_varlen")
URL = "https://rest.ensembl.org/lookup/id"


def post(ids):
    req = urllib.request.Request(
        URL, data=json.dumps({"ids": ids}).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def main():
    ids = [l.strip() for l in open(f"{B}/top_block_atids.txt") if l.strip()]
    out = {}
    for i in range(0, len(ids), 40):
        chunk = ids[i:i + 40]
        for attempt in range(3):
            try:
                out.update(post(chunk)); break
            except Exception as e:
                print(f"  retry {attempt} ({e})"); time.sleep(3)
    rows = ["atid,symbol,biotype,description"]
    for a in ids:
        d = out.get(a) or {}
        sym = d.get("display_name", "") or ""
        bt = d.get("biotype", "") or ""
        desc = (d.get("description", "") or "").replace(",", ";")
        rows.append(f"{a},{sym},{bt},{desc}")
    open(f"{B}/atid_symbols.csv", "w").write("\n".join(rows) + "\n")
    named = sum(1 for a in ids if (out.get(a) or {}).get("display_name"))
    print(f"[done] {named}/{len(ids)} AT-IDs have a symbol -> {B}/atid_symbols.csv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Independent cross-check of the block->gene region mapping.

Our `build_significant_genes.py` assigned genes to each significant block's span using
the LOCAL TAIR10 GFF (lib.load_genes overlap). Here we re-query each block's genomic
span against the **Ensembl REST `overlap/region`** endpoint (Araport11 annotation — a
DIFFERENT annotation source and a different code path) and reconcile the gene sets, so
any overlap error or annotation-version difference is surfaced explicitly.

Reads:  significant_genes_gen{gen}_{climate}_{regime}.csv  (block, chrom, start, end, gene)
Writes: gene_mapping_verification_gen{gen}_{climate}_{regime}.csv  per block:
        n_local, n_ensembl, n_agree, local_only, ensembl_only, status

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python   # needs outbound HTTPS
  $PY analysis/grenenet_gea/phase1_replication/verify_gene_mapping.py
"""
from __future__ import annotations
import argparse, json, os, sys, time, urllib.request
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OVERLAP = "https://rest.ensembl.org/overlap/region/arabidopsis_thaliana"


def ensembl_region_genes(chrom, start, end, retries=4):
    """Set of protein-coding gene AT-ids overlapping chrom:start-end (Ensembl/Araport11)."""
    c = str(chrom).replace("Chr", "")
    url = f"{OVERLAP}/{c}:{int(start)}-{int(end)}?feature=gene;content-type=application/json"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    for a in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                d = json.load(r)
            # keep AT-ids; note biotype so non-coding don't masquerade as misses
            return {g["id"]: g.get("biotype", "") for g in d if g.get("id", "").startswith("AT")}
        except Exception as e:
            if a == retries - 1:
                print(f"  [ensembl] {c}:{start}-{end} failed: {e}", flush=True)
            time.sleep(2)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regime", default="deg7cap2000")
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--gen", type=int, default=9)
    args = ap.parse_args()
    base = f"{lib.GEA}/phase1_replication"
    src = f"{base}/significant_genes_gen{args.gen}_{args.climate}_{args.regime}.csv"
    d = pd.read_csv(src)

    # local gene set per block (drop the intergenic gene="" placeholder)
    blk = (d[d.gene.notna() & (d.gene.astype(str) != "")]
           .groupby("block").agg(chrom=("chrom", "first"), start=("start", "first"),
                                  end=("end", "first"), local=("gene", lambda s: set(s))))
    # blocks with NO local gene (intergenic) still get checked against Ensembl
    empty = d[~d.block.isin(blk.index)].drop_duplicates("block").set_index("block")
    for b, r in empty.iterrows():
        if pd.notna(r.chrom):
            blk.loc[b] = [r.chrom, r.start, r.end, set()]

    print(f"verifying {len(blk)} significant blocks against Ensembl overlap/region…", flush=True)
    rows = []
    for b, r in blk.iterrows():
        ens = ensembl_region_genes(r.chrom, r.start, r.end)
        time.sleep(0.2)                              # be polite to the API
        if ens is None:
            rows.append(dict(block=b, status="API_FAIL")); continue
        ens_pc = {g for g, bt in ens.items() if bt == "protein_coding"}
        loc = set(r.local)
        agree = loc & ens_pc
        local_only = loc - ens_pc
        ens_only = ens_pc - loc
        status = "OK" if not local_only and not ens_only else "DIFF"
        rows.append(dict(block=b, chrom=r.chrom, start=int(r.start), end=int(r.end),
                         n_local=len(loc), n_ensembl_pc=len(ens_pc), n_agree=len(agree),
                         local_only=";".join(sorted(local_only)),
                         ensembl_only=";".join(sorted(ens_only)), status=status))
    v = pd.DataFrame(rows)
    out = f"{base}/gene_mapping_verification_gen{args.gen}_{args.climate}_{args.regime}.csv"
    v.to_csv(out, index=False)

    nok = (v.status == "OK").sum(); ndiff = (v.status == "DIFF").sum()
    nfail = (v.status == "API_FAIL").sum()
    tot_local = v.n_local.sum(); tot_agree = v.n_agree.sum()
    print(f"\n=== reconciliation ===")
    print(f"blocks: {len(v)} | exact-match: {nok} | differ: {ndiff} | api-fail: {nfail}")
    print(f"local genes total: {tot_local} | confirmed by Ensembl: {tot_agree} "
          f"({100*tot_agree/max(tot_local,1):.1f}%)")
    print(f"local-only (in our TAIR10 map, not in Ensembl/Araport11 region): "
          f"{sum(len(x.split(';')) for x in v.local_only.dropna() if x)}")
    print(f"ensembl-only (Araport11 gene in region we did NOT list): "
          f"{sum(len(x.split(';')) for x in v.ensembl_only.dropna() if x)}")
    if ndiff:
        print("\nblocks that differ:")
        print(v[v.status == "DIFF"][["block", "chrom", "start", "end", "n_local",
              "n_ensembl_pc", "local_only", "ensembl_only"]].to_string(index=False))
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()

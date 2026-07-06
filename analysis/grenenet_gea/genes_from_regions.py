#!/usr/bin/env python
"""Genes overlapping genomic regions, fetched DIRECTLY from the Ensembl Plants REST API.

One step: region -> genes (id, symbol, biotype, description), no local GFF, no separate
symbol-resolution pass. Uses the `overlap/region` endpoint for A. thaliana (TAIR10),
which carries gene symbols natively. arabidopsis.org is 403-gated on the cluster; Ensembl
Plants REST needs no auth and Savio compute nodes have outbound HTTPS.

Usage
-----
  # from a block CSV with chrom/unit_start/unit_end (the GEA/WZA outputs):
  python genes_from_regions.py --csv <blocks.csv> [--n 15] [--sort block_p] [--out genes.csv]
  # or ad-hoc regions:
  python genes_from_regions.py --region Chr3:5823895-5830455 --region Chr1:611128-622074

Importable: `genes_in_region("Chr3", 5823895, 5830455)` -> DataFrame.
"""
import argparse, json, sys, time, urllib.request
import pandas as pd

SERVER = "https://rest.ensembl.org"
SPECIES = "arabidopsis_thaliana"
_CACHE = {}


def _get(path):
    req = urllib.request.Request(SERVER + path,
                                 headers={"Content-Type": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 429:                       # rate-limited: honour Retry-After
                time.sleep(float(e.headers.get("Retry-After", 2)) + 0.5); continue
            if attempt == 3:
                raise
            time.sleep(2)
        except Exception:
            if attempt == 3:
                raise
            time.sleep(2)


def genes_in_region(chrom, start, end) -> pd.DataFrame:
    """All protein-coding/ncRNA genes overlapping chrom:start-end (Ensembl Plants)."""
    chrom = str(chrom).replace("Chr", "").replace("chr", "")  # Ensembl uses 1..5
    key = (chrom, int(start), int(end))
    if key in _CACHE:
        return _CACHE[key]
    feats = _get(f"/overlap/region/{SPECIES}/{chrom}:{int(start)}-{int(end)}"
                 f"?feature=gene;content-type=application/json") or []
    time.sleep(0.08)                                # stay under ~15 req/s
    df = pd.DataFrame([{
        "gene_id": f.get("gene_id") or f.get("id"),
        "symbol": f.get("external_name") or "",
        "biotype": f.get("biotype") or "",
        "g_start": f.get("start"), "g_end": f.get("end"),
        "strand": f.get("strand"),
        "description": (f.get("description") or "").split(" [")[0],
    } for f in feats]).sort_values("g_start").reset_index(drop=True) \
        if feats else pd.DataFrame(
            columns=["gene_id", "symbol", "biotype", "g_start", "g_end", "strand", "description"])
    _CACHE[key] = df
    return df


def label(df):
    """Compact 'symbol-or-id' list for a region's gene table."""
    return ";".join((df.symbol.where(df.symbol != "", df.gene_id)).tolist())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", help="block CSV with chrom/unit_start/unit_end")
    ap.add_argument("--n", type=int, default=15, help="top N rows from --csv")
    ap.add_argument("--sort", default="block_p", help="column to sort --csv ascending by")
    ap.add_argument("--region", action="append", default=[], help="Chr:start-end (repeatable)")
    ap.add_argument("--out", help="write per-block gene table here")
    args = ap.parse_args()

    regions = []  # (unit, chrom, start, end, extra-dict)
    if args.csv:
        d = pd.read_csv(args.csv)
        if args.sort in d:
            d = d.sort_values(args.sort)
        d = d.head(args.n)
        keep = [c for c in ("block_p", "q", "wza", "best_beta1", "lead_beta1") if c in d]
        for _, r in d.iterrows():
            regions.append((r.get("unit", ""), r.chrom, int(r.unit_start),
                            int(r.unit_end), {c: r[c] for c in keep}))
    for rg in args.region:
        ch, se = rg.split(":"); s, e = se.split("-")
        regions.append((rg, ch, int(s), int(e), {}))

    rows = []
    pd.set_option("display.width", 200, "display.max_colwidth", 90)
    for unit, ch, s, e, extra in regions:
        g = genes_in_region(ch, s, e)
        named = g[g.symbol != ""]
        print(f"\n{unit or f'{ch}:{s}-{e}'}  ({(e-s)/1e3:.2f} kb)  "
              + "  ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                          for k, v in extra.items()))
        if len(g):
            for _, r in g.iterrows():
                print(f"    {r.gene_id:14s} {r.symbol:12s} {r.biotype:16s} {r.description}")
        else:
            print("    (no genes overlap)")
        rows.append({"unit": unit, "chrom": ch, "start": s, "end": e,
                     **extra, "n_genes": len(g), "n_named": len(named),
                     "genes": label(g),
                     "symbols": ";".join(named.symbol.tolist())})
    out = pd.DataFrame(rows)
    if args.out:
        out.to_csv(args.out, index=False)
        print(f"\n[done] {len(out)} regions -> {args.out}")


if __name__ == "__main__":
    main()

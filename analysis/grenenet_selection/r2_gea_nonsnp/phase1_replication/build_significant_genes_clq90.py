#!/usr/bin/env python
"""BH-significant clq0.9 LD blocks -> gene table, for the clq90 phase-1 GEA replication.

For every WZA output in {kendall,lfmm,binomial} x {snp,nonsnp} at gen9 / bio1 /
deg2 (canonical, primary regime for the clq90 run — see STATUS_clq90.md), compute
BH-FDR on the block p (Z_pVal) and take the blocks with q < --fdr. For each
significant block:
  - reconstruct its genomic SPAN from the clq0.9 BigLD interval TSVs
    (analysis/grenenet_selection/blocks/results/blocks_mcf90/chr{N}_clq0.9_blocks_clq0.9.tsv),
  - list every TAIR10 gene overlapping that span (local GFF, lib.load_genes),
  - enrich each gene AT-id with its SYMBOL + DESCRIPTION via the Ensembl Plants
    REST API (TAIR10 assembly; no auth, batch POST).

This is the clq0.9-block counterpart of build_significant_genes.py (which uses
the coarse phase-1 hapFIRE blocks / snp+smallindel+sv / deg7cap2000).

Output (analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/clq90/):
  significant_genes_clq90_gen9_bio1_deg2.csv  one row per (block, gene):
    gene, symbol, description, block, region, chrom, start, end, span_bp,
    n_variants, n_genes_in_block, n_sig_combos, sig_models, sig_classes,
    sig_combos, + 6 indicator cols {model}_{class} = BH q-value where significant.

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python   # needs outbound HTTPS
  $PY analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/build_significant_genes_clq90.py --fdr 0.05
"""
from __future__ import annotations
import argparse, json, os, sys, time, urllib.request
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

WDIR = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/clq90/wza"
OUTDIR = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/clq90"
MODELS = ["kendall", "lfmm", "binomial"]
CLASSES = ["snp", "nonsnp"]
ENSEMBL = "https://rest.ensembl.org/lookup/id"


def bh_fdr(p):
    p = np.asarray(p, float); n = np.isfinite(p).sum()
    q = np.full(len(p), np.nan)
    ok = np.where(np.isfinite(p))[0]
    o = ok[np.argsort(p[ok])]
    ranked = p[o] * n / (np.arange(1, len(o) + 1))
    q[o] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    return q


def block_spans(r2: float = 0.9):
    """chrom('Chr#'), start, end, n_variants for every clq{r2} BigLD block.

    Block id = 'Chr{n}_{idx}', idx = the block's rank by start_pos within its
    chrom (0-based) — must match lib.assign_clq_blocks exactly.
    """
    tag = f"clq{r2}"
    rows = []
    for ci in range(1, 6):
        c = f"Chr{ci}"
        f = f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv"
        g = pd.read_csv(f, sep="\t").sort_values("start_pos").reset_index(drop=True)
        for idx, row in g.iterrows():
            rows.append((f"{c}_{idx}", c, int(row.start_pos), int(row.end_pos),
                        int(row.n_variants)))
    return pd.DataFrame(rows, columns=["block", "chrom", "start", "end", "n_variants"]
                        ).set_index("block")


def ensembl_symbols(ids, chunk=900, retries=3):
    """AT-id -> (symbol, description) via Ensembl REST batch POST. Robust to failures."""
    out = {}
    for i in range(0, len(ids), chunk):
        sub = ids[i:i + chunk]
        body = json.dumps({"ids": sub}).encode()
        req = urllib.request.Request(
            ENSEMBL, data=body,
            headers={"Content-Type": "application/json", "Accept": "application/json"})
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    d = json.load(r)
                for k, v in d.items():
                    if v:
                        out[k] = (v.get("display_name", "") or "",
                                  (v.get("description", "") or "").split(" [Source")[0])
                break
            except Exception as e:
                print(f"  [ensembl] chunk {i} attempt {attempt+1} failed: {e}", flush=True)
                time.sleep(2)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--regime", default="deg2")
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--r2", type=float, default=0.9)
    ap.add_argument("--no-api", action="store_true", help="skip Ensembl symbol enrichment")
    ap.add_argument("--q-floor", type=float, default=1e-16,
                    help="floor for reported BH q (WZA Z_pVal underflows float64 to 0 / rounds "
                         "to ~0 for extreme high-LD blocks); values <=floor are clipped and flagged")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    FLOOR = args.q_floor

    # 1) BH-significant blocks per (model, class): {block: q}
    sig = {}                       # (model, class) -> {block: q}
    for m in MODELS:
        for c in CLASSES:
            f = f"{WDIR}/wza_{m}_{c}_gen{args.gen}_{args.climate}_{args.regime}.csv"
            w = pd.read_csv(f).rename(columns={"index": "block"})
            w["block"] = w["block"].astype(str)
            w = w[w["Z_pVal"].notna()].copy()
            w["q"] = bh_fdr(w["Z_pVal"].to_numpy())
            s = w[w["q"] < args.fdr]
            sig[(m, c)] = dict(zip(s["block"], s["q"]))
            print(f"{m:9s} {c:8s}: {len(w):,} blocks | BH q<{args.fdr}: {len(s)}", flush=True)

    # 2) union of significant blocks, with which combos flagged each
    all_blocks = sorted({b for d in sig.values() for b in d})
    print(f"\nUnion of BH-significant blocks across all {len(MODELS)*len(CLASSES)} combos: "
          f"{len(all_blocks)}", flush=True)
    if not all_blocks:
        sys.exit("no significant blocks at this FDR")

    spans = block_spans(args.r2)
    genes = lib.load_genes()

    # 3) span + overlapping genes per block; collect all AT-ids
    block_genes = {}               # block -> list of AT-ids
    block_meta = {}
    for blk in all_blocks:
        if blk not in spans.index:
            block_meta[blk] = (None, None, None, None); block_genes[blk] = []
            continue
        sp = spans.loc[blk]
        gc = genes[(genes.chrom == sp.chrom) & (genes.end >= sp.start) & (genes.start <= sp.end)]
        block_meta[blk] = (sp.chrom, int(sp.start), int(sp.end), int(sp.n_variants))
        block_genes[blk] = list(gc.gene)

    uniq_genes = sorted({g for gs in block_genes.values() for g in gs})
    print(f"Unique TAIR10 genes overlapping significant blocks: {len(uniq_genes)}", flush=True)

    # 4) Ensembl symbol/description enrichment
    sym = {}
    if uniq_genes and not args.no_api:
        print(f"Querying Ensembl REST for {len(uniq_genes)} gene symbols…", flush=True)
        sym = ensembl_symbols(uniq_genes)
        print(f"  got {len(sym)}/{len(uniq_genes)} symbols", flush=True)

    # 5) build the gene x combo table
    combo_cols = [f"{m}_{c}" for m in MODELS for c in CLASSES]
    rows = []
    for blk in all_blocks:
        ch, st, en, nvar = block_meta[blk]
        flagged = {(m, c): sig[(m, c)][blk] for m in MODELS for c in CLASSES if blk in sig[(m, c)]}
        sig_models = sorted({m for (m, c) in flagged})
        sig_classes = sorted({c for (m, c) in flagged})
        sig_combos = ";".join(f"{m}:{c}" for (m, c) in sorted(flagged))
        gids = block_genes[blk] or [""]      # keep block even if intergenic (gene="")
        for gid in gids:
            symbol, desc = sym.get(gid, ("", ""))
            region = f"{ch}:{st}-{en}" if st is not None else ""   # browser-pastable (TAIR10)
            # floor the BH q (WZA Z_pVal underflows to 0 / is untrustworthy in magnitude for
            # extreme high-LD blocks); flag any combo whose raw q was at/below the floor.
            q_floored = any(flagged[(m, c)] <= FLOOR for (m, c) in flagged)
            row = dict(gene=gid, symbol=symbol, description=desc, block=blk,
                       region=region, chrom=ch, start=st, end=en,
                       span_bp=(en - st) if st is not None else None,
                       n_variants=nvar, n_genes_in_block=len(block_genes[blk]),
                       n_sig_combos=len(flagged), sig_models=";".join(sig_models),
                       sig_classes=";".join(sig_classes), sig_combos=sig_combos,
                       q_floored=q_floored)
            for m in MODELS:
                for c in CLASSES:
                    row[f"{m}_{c}"] = max(flagged[(m, c)], FLOOR) if (m, c) in flagged else ""
            rows.append(row)

    out = pd.DataFrame(rows)
    # sort: most-supported blocks first, then by block, gene
    out = out.sort_values(["n_sig_combos", "block", "gene"], ascending=[False, True, True])
    path = args.out or (f"{OUTDIR}/significant_genes_clq90_gen{args.gen}_{args.climate}_"
                        f"{args.regime}.csv")
    out.to_csv(path, index=False)
    named = out[out.symbol != ""]
    print(f"\n-> {len(out)} (block,gene) rows | {out.block.nunique()} blocks | "
          f"{out[out.gene!=''].gene.nunique()} genes ({named.symbol.nunique()} with symbols)")
    print(f"   written: {path}")
    # quick peek: the multi-supported hits
    top = out[out.n_sig_combos >= 2][["symbol", "gene", "block", "n_sig_combos", "sig_combos"]].head(20)
    if len(top):
        print("\n   blocks significant in >=2 combos (top 20 rows):")
        print(top.to_string(index=False))


if __name__ == "__main__":
    main()

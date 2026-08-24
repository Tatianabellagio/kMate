#!/usr/bin/env python
"""SUPERSEDED (2026-07-27) — reads deg-2 WZA output + the retired pooled 2-class
(snp/nonsnp) split. deg-2 fabricates significance on sparse large blocks (see
`../STATUS_clq90.md`'s 2026-07-21/07-27 banners); the pipeline moved to isotonic
+ the 3-class (snp/sv/smallindel) split. Current equivalent: the `class_specific()`
helper in `_build_snp_vs_nonsnp_peaks_nb.py` -> `{cls}_specific_peaks_{bonf,fdr}.csv`
+ `snp_vs_nonsnp_new_peaks.ipynb`. Kept for provenance/archaeology; do not cite its
output (`nonsnp_specific_genes.csv`, the "437 blocks -> 548 genes, heat-stress
dominates" list) as current.

Genes on nonSNP-specific hit blocks (BH-sig in non-SNP but NOT in SNP), all axes.

A block is "nonSNP-specific" for a given (axis, model) if it is WZA BH-FDR q<0.05 in
the NON-SNP class but NOT in the SNP class for that same axis+model. We collect the
UNION of such blocks across all 20 axes x 3 models, record which (axis,model) flagged
each (and whether the same block is EVER a SNP hit anywhere — a stricter "never a SNP
hit" flag), map each block to its clq0.9 genomic span, list overlapping TAIR10 genes,
and enrich AT-ids with SYMBOL+DESCRIPTION via the Ensembl Plants REST API.

Output: multiaxis/nonsnp_specific_genes.csv  (one row per block,gene)

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python   # needs outbound HTTPS
  $PY nonsnp_specific_genes.py --fdr 0.05
"""
from __future__ import annotations
import argparse, json, os, sys, time, urllib.request
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

WD = f"{lib.GEA}/phase1_replication/results/multiaxis/wza"
OUTDIR = f"{lib.GEA}/phase1_replication/results/multiaxis"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
MODELS = ["kendall", "lfmm", "binomial"]
ENSEMBL = "https://rest.ensembl.org/lookup/id"


def bh_sig(model, cls, axis, q=0.05):
    f = f"{WD}/wza_{model}_{cls}_gen9_{axis}_deg2.csv"
    if not os.path.exists(f):
        return None
    w = pd.read_csv(f).rename(columns={"index": "block"})
    w["block"] = w["block"].astype(str)
    w = w[w["Z_pVal"].notna()].copy()
    p = w["Z_pVal"].to_numpy(); n = len(p); o = np.argsort(p); qv = np.empty(n)
    qv[o] = (p[o] * n) / (np.arange(n) + 1); qv[o] = np.minimum.accumulate(qv[o][::-1])[::-1]
    return set(w.loc[np.clip(qv, 0, 1) < q, "block"])


def block_spans(r2=0.9):
    tag = f"clq{r2}"; rows = []
    for ci in range(1, 6):
        g = pd.read_csv(f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv", sep="\t"
                        ).sort_values("start_pos").reset_index(drop=True)
        for idx, r in g.iterrows():
            rows.append((f"Chr{ci}_{idx}", f"Chr{ci}", int(r.start_pos), int(r.end_pos), int(r.n_variants)))
    return pd.DataFrame(rows, columns=["block", "chrom", "start", "end", "n_variants"]).set_index("block")


def ensembl_symbols(ids, chunk=900, retries=3):
    out = {}
    for i in range(0, len(ids), chunk):
        sub = ids[i:i + chunk]
        req = urllib.request.Request(ENSEMBL, data=json.dumps({"ids": sub}).encode(),
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
                print(f"  [ensembl] chunk {i} attempt {attempt+1}: {e}", flush=True); time.sleep(2)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--no-api", action="store_true")
    args = ap.parse_args()

    # 1) collect nonSNP-specific flags per (axis, model); and every-SNP-hit set
    ns_flags = {}          # block -> list of "axis:model" where nonsnp-sig & snp-not
    snp_ever = set()       # block BH-sig in SNP for ANY axis/model
    nonsnp_ever = set()
    for axis in AXES:
        for model in MODELS:
            s = bh_sig(model, "snp", axis, args.fdr); n = bh_sig(model, "nonsnp", axis, args.fdr)
            if s is None or n is None:
                continue
            snp_ever |= s; nonsnp_ever |= n
            for b in (n - s):
                ns_flags.setdefault(b, []).append(f"{axis}:{model}")
    blocks = sorted(ns_flags)
    print(f"nonSNP-specific hit blocks (nonsnp-sig & snp-not, some axis+model): {len(blocks)}")
    print(f"  of which NEVER a SNP hit anywhere: {len([b for b in blocks if b not in snp_ever])}")
    if not blocks:
        sys.exit("no nonSNP-specific blocks")

    spans = block_spans(); genes = lib.load_genes()
    block_genes, meta = {}, {}
    for b in blocks:
        if b not in spans.index:
            meta[b] = (None, None, None, None); block_genes[b] = []; continue
        sp = spans.loc[b]
        gc = genes[(genes.chrom == sp.chrom) & (genes.end >= sp.start) & (genes.start <= sp.end)]
        meta[b] = (sp.chrom, int(sp.start), int(sp.end), int(sp.n_variants)); block_genes[b] = list(gc.gene)

    uniq = sorted({g for gs in block_genes.values() for g in gs})
    sym = {} if args.no_api else ensembl_symbols(uniq)
    print(f"unique genes on nonSNP-specific blocks: {len(uniq)} ({len(sym)} with symbols)")

    rows = []
    for b in blocks:
        ch, st, en, nv = meta[b]
        flags = ns_flags[b]
        axes_hit = sorted({f.split(':')[0] for f in flags})
        models_hit = sorted({f.split(':')[1] for f in flags})
        for gid in (block_genes[b] or [""]):
            symbol, desc = sym.get(gid, ("", ""))
            rows.append(dict(block=b, region=f"{ch}:{st}-{en}" if st else "", chrom=ch, start=st, end=en,
                             n_variants=nv, gene=gid, symbol=symbol, description=desc,
                             n_flags=len(flags), axes_nonsnp_specific=";".join(axes_hit),
                             models=";".join(models_hit), flags=";".join(flags),
                             never_snp_hit=(b not in snp_ever)))
    out = pd.DataFrame(rows).sort_values(["n_flags", "block", "gene"], ascending=[False, True, True])
    path = f"{OUTDIR}/nonsnp_specific_genes.csv"
    out.to_csv(path, index=False)
    print(f"\n-> {len(out)} (block,gene) rows | {out.block.nunique()} blocks | "
          f"{out[out.gene!=''].gene.nunique()} genes -> {path}")
    top = out[out.symbol != ""].drop_duplicates("gene").head(25)
    print("\ntop nonSNP-specific genes (by #axis-model flags):")
    print(top[["block", "symbol", "gene", "n_flags", "axes_nonsnp_specific", "never_snp_hit"]].to_string(index=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Corrected (2026-07-27) twin of `nonsnp_specific_genes.py` — symbol/description/
functional-category annotation of the class-specific candidate genes, but on the
CURRENT isotonic + 3-class (snp/sv/smallindel) outputs instead of the retired
deg-2 + pooled-nonsnp ones. Same Ensembl Plants REST lookup + same CATS/CURATED
keyword scheme as the old script, so the two runs are directly comparable.

Source: `{nonsnp,sv,smallindel}_specific_peaks_{bonf,fdr}.csv` (from
`_build_snp_vs_nonsnp_peaks_nb.py`'s `class_specific()`, already on disk).
Output: `multiaxis/candidate_genes_corrected.csv` (one row per gene, union over
classes at the FDR tier, with per-class/tier provenance).

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY candidate_genes_corrected.py
"""
from __future__ import annotations
import json, os, sys, time, urllib.request
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import lib

MA = f"{lib.GEA}/phase1_replication/results/multiaxis"
CLASSES = ["nonsnp", "sv", "smallindel"]
ENSEMBL = "https://rest.ensembl.org/lookup/id"


def block_spans(r2=0.9):
    tag = f"clq{r2}"; rows = []
    for ci in range(1, 6):
        g = pd.read_csv(f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv", sep="\t"
                        ).sort_values("start_pos").reset_index(drop=True)
        for idx, r in g.iterrows():
            rows.append((f"Chr{ci}_{idx}", f"Chr{ci}", int(r.start_pos), int(r.end_pos)))
    return pd.DataFrame(rows, columns=["block", "chrom", "start", "end"]).set_index("block")


def genes_on(block, spans, genes):
    if block not in spans.index:
        return []
    s = spans.loc[block]
    gc = genes[(genes.chrom == s.chrom) & (genes.end >= s.start) & (genes.start <= s.end)]
    return list(gc.gene)


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


CATS = {
    "heat":        r"heat|thermo|high temperature|chaperone|heat shock|hsp|hsf|dnaj",
    "drought_ABA": r"drought|abscisic|\baba\b|dehydrat|desicc|osmotic|water stress|stomat|dehydrin|late embryogenesis|proline|salt",
    "cold":        r"\bcold\b|freezing|chilling|\bcbf\b|frost|cor15",
    "flower_circ": r"flower|floral|photoperiod|vernaliz|infloresc|florigen|circadian|clock|rhythm|oscillat",
}
CURATED = {
    "HSFA2": "heat", "HSBP": "heat", "HSP17.4": "heat", "MBF1C": "heat",
    "RAS1": "drought_ABA", "DREB2A": "drought_ABA", "RD29A": "drought_ABA", "NCED3": "drought_ABA",
    "COR15A": "cold", "VRN2": "flower_circ", "FT": "flower_circ", "FLC": "flower_circ",
    "GI": "flower_circ", "CCA1": "flower_circ", "TOC1": "flower_circ", "PIF4": "flower_circ",
}


def main():
    spans = block_spans(); genes_df = lib.load_genes()

    # 1) collect block-level rows (union over classes, FDR tier = superset) with provenance
    block_rows = {}  # block -> dict(classes:set, tiers:dict(cls->tier), n_flags(max), axes:set, models:set, never_snp_hit)
    for cls in CLASSES:
        for tier in ["fdr", "bonf"]:
            f = f"{MA}/{cls}_specific_peaks_{tier}.csv"
            if not os.path.exists(f):
                print(f"  SKIP missing {f}"); continue
            df = pd.read_csv(f)
            for _, r in df.iterrows():
                b = r.block
                e = block_rows.setdefault(b, dict(classes=set(), tiers={}, n_flags=0,
                                                   axes=set(), models=set(), never_snp_hit=False))
                e["classes"].add(cls)
                e["tiers"][cls] = tier if e["tiers"].get(cls) != "bonf" else "bonf"  # bonf wins if seen at all
                e["n_flags"] = max(e["n_flags"], int(r.n_flags))
                rax, rmo = r["axes"], r["models"]  # NB: `.axes` shadows the pandas Series builtin
                e["axes"] |= set(str(rax).split(";")) if pd.notna(rax) else set()
                e["models"] |= set(str(rmo).split(";")) if pd.notna(rmo) else set()
                e["never_snp_hit"] = e["never_snp_hit"] or bool(r.never_snp_hit)

    print(f"union blocks (any class, FDR or Bonf): {len(block_rows)}")

    # 2) expand blocks -> genes (full overlap, not truncated)
    gene_rows = []
    for b, e in block_rows.items():
        gs = genes_on(b, spans, genes_df)
        for gid in (gs or [""]):
            gene_rows.append(dict(block=b, gene=gid, classes=";".join(sorted(e["classes"])),
                                   tiers=";".join(f"{c}:{t}" for c, t in sorted(e["tiers"].items())),
                                   n_flags=e["n_flags"], axes=";".join(sorted(e["axes"])),
                                   models=";".join(sorted(e["models"])), never_snp_hit=e["never_snp_hit"]))
    d = pd.DataFrame(gene_rows)
    d = d[d.gene != ""]
    print(f"union blocks -> {d.block.nunique()} blocks -> {d.gene.nunique()} unique genes")

    # 3) one row per gene (keep the block with the most flags per gene)
    g = (d.sort_values("n_flags", ascending=False).drop_duplicates("gene").reset_index(drop=True))

    # 4) Ensembl symbol + description
    uniq = sorted(g.gene.unique())
    sym = ensembl_symbols(uniq)
    g["symbol"] = g.gene.map(lambda x: sym.get(x, ("", ""))[0])
    g["description"] = g.gene.map(lambda x: sym.get(x, ("", ""))[1])
    print(f"{len(g)} unique genes ({(g.symbol != '').sum()} with symbols)")

    # 5) functional-category tagging (identical scheme to the old script)
    txt = (g.symbol.fillna("") + " " + g.description.fillna("")).str.lower()
    g["category"] = ""
    for cat, pat in CATS.items():
        g.loc[(g.category == "") & txt.str.contains(pat, regex=True, na=False), "category"] = cat
    for symn, cat in CURATED.items():
        g.loc[g.symbol == symn, "category"] = cat
    print("\ncategory counts:\n", g["category"].replace("", "other").value_counts().to_string())

    out = f"{MA}/candidate_genes_corrected.csv"
    g[["block", "gene", "symbol", "description", "category", "n_flags", "classes", "tiers",
       "axes", "models", "never_snp_hit"]].to_csv(out, index=False)
    print(f"\n-> {out}")

    for cat in ["heat", "drought_ABA", "cold", "flower_circ"]:
        sub = g[g.category == cat].sort_values("n_flags", ascending=False)
        print(f"\n===== {cat}: {len(sub)} genes =====")
        if len(sub):
            print(sub[["block", "symbol", "gene", "description", "n_flags", "classes", "never_snp_hit"]]
                  .to_string(index=False))


if __name__ == "__main__":
    main()

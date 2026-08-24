#!/usr/bin/env python
"""Compare the four block-partition arms: significant blocks + genes per arm.

Arms (all same current data, same isotonic WZA; only the partition differs):
  a09s  clq0.9 strict   (production baseline -- discards 40-49% of records)
  a09t  clq0.9 tiling   (same LD boundaries, HapFM gap-free rule, 0% discarded)
  a05s  clq0.5 strict
  a05t  clq0.5 tiling

a09s->a09t isolates the DATA-LOSS effect; a09t->a05t isolates BLOCK COARSENESS.

Block spans are computed per-arm (tiling blocks are extended leftward to absorb the
preceding gap, so their gene overlap differs from the strict interval).

Output: multiaxis/block_arms_{summary,genes}.csv
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
import lib

MA = f"{lib.GEA}/phase1_replication/results/multiaxis"
WD = f"{MA}/wza_arms"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
MODELS = ["kendall", "lfmm", "binomial"]
CLASSES = ["snp", "sv", "smallindel", "nonsnp"]
ARMS = {"a09s": (0.9, "strict"), "a09t": (0.9, "tiling"),
        "a05s": (0.5, "strict"), "a05t": (0.5, "tiling")}


def bh(p):
    p = np.asarray(p, float); n = len(p)
    if n == 0: return p
    o = np.argsort(p); q = np.empty(n)
    q[o] = (p[o] * n) / (np.arange(n) + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)


def spans_for(r2, how):
    """Block id -> (chrom, start, end). Tiling blocks absorb the preceding gap."""
    out = {}
    for ci in range(1, 6):
        ch = f"Chr{ci}"; tag = f"clq{r2}"
        g = pd.read_csv(f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv", sep="\t"
                        ).sort_values("start_pos").reset_index(drop=True)
        ends = g.end_pos.to_numpy(np.int64); starts = g.start_pos.to_numpy(np.int64)
        for i in range(len(g)):
            if how == "tiling":
                lo = 1 if i == 0 else int(ends[i-1]) + 1
                hi = int(ends[i]) if i < len(g)-1 else 10**9
            else:
                lo, hi = int(starts[i]), int(ends[i])
            out[f"{ch}_{i}"] = (ch, lo, hi)
    return out


def main():
    GENES = lib.load_genes()
    span_cache = {arm: spans_for(*ARMS[arm]) for arm in ARMS}
    gene_cache = {}

    def genes_on(arm, blk):
        key = (arm, blk)
        if key in gene_cache: return gene_cache[key]
        sp = span_cache[arm].get(blk)
        if sp is None:
            gene_cache[key] = []; return []
        ch, lo, hi = sp
        g = GENES[(GENES.chrom == ch) & (GENES.end >= lo) & (GENES.start <= hi)]
        gene_cache[key] = list(g.gene)
        return gene_cache[key]

    summary, generows = [], []
    for arm in ARMS:
        for model in MODELS:
            for cls in CLASSES:
                for axis in AXES:
                    f = f"{WD}/wza_{model}_{cls}_gen9_{axis}_{arm}.csv"
                    if not os.path.exists(f):
                        continue
                    w = pd.read_csv(f).rename(columns={"index": "block"})
                    w["block"] = w["block"].astype(str)
                    w = w[w.Z_pVal.notna() & (w.Z_pVal > 0)]
                    n = len(w)
                    if n == 0: continue
                    q = bh(w.Z_pVal.to_numpy())
                    bonf = w.Z_pVal < 0.05 / n
                    fdr = q < 0.05
                    summary.append(dict(arm=arm, model=model, cls=cls, axis=axis,
                                        n_blocks=n, n_bonf=int(bonf.sum()), n_fdr=int(fdr.sum())))
                    for tier, mask in [("bonf", bonf), ("fdr", fdr)]:
                        for _, r in w[mask].iterrows():
                            for g in genes_on(arm, r["block"]):
                                generows.append(dict(arm=arm, model=model, cls=cls, axis=axis,
                                                     tier=tier, block=r["block"], gene=g,
                                                     Z=r["Z"], p=r["Z_pVal"]))
    S = pd.DataFrame(summary); G = pd.DataFrame(generows)
    S.to_csv(f"{MA}/block_arms_summary.csv", index=False)
    G.to_csv(f"{MA}/block_arms_genes.csv", index=False)

    print("=== significant BLOCKS summed over 20 axes, per arm x model x class ===\n")
    for tier, col in [("Bonferroni", "n_bonf"), ("FDR q<0.05", "n_fdr")]:
        print(f"--- {tier} ---")
        piv = S.pivot_table(index=["model", "cls"], columns="arm", values=col, aggfunc="sum")
        piv = piv.reindex(columns=["a09s", "a09t", "a05s", "a05t"])
        print(piv.fillna(0).astype(int).to_string()); print()

    print("\n=== unique GENES implicated per arm (union over axes/models) ===\n")
    print(f"{'arm':6s} {'tier':5s} | " + " ".join(f"{c:>12s}" for c in CLASSES))
    print("-" * 70)
    for arm in ["a09s", "a09t", "a05s", "a05t"]:
        for tier in ["bonf", "fdr"]:
            sub = G[(G.arm == arm) & (G.tier == tier)]
            cnt = [sub[sub.cls == c].gene.nunique() for c in CLASSES]
            print(f"{arm:6s} {tier:5s} | " + " ".join(f"{v:12,}" for v in cnt))

    print(f"\n-> {MA}/block_arms_summary.csv\n-> {MA}/block_arms_genes.csv")


if __name__ == "__main__":
    main()

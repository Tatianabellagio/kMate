#!/usr/bin/env python
"""One tidy table behind notebooks/gea_hit_landscape.ipynb: every record the non-SNP GEA
tested, what it is, where it lands, and on which axes the GEA called it.

Nothing is re-estimated. Inputs are the GEA's own outputs and annotations:

  wza_in_clq09_tile/lfmm_{nonsnp,sv,smallindel}_gen9_<axis>.csv   per-record LFMM p, MAF, block
  results/gea_universe_regions.csv.gz   region per record, the hits' classifier on every
                                        tested record (annotate_gea_universe.py)

The hit call is build_gea_pool.py's, reproduced exactly: per class x axis, MAF > 0.05 and a
block, Bonferroni at 0.05 / n_records; a record is a hit on an axis if it clears that bar in
any of the three class scans (sv, smallindel, pooled nonsnp). Records are identified by
(chrom, pos, ref_len, alt_len, MAF) -- one position can carry several records with the same
lengths (the multiallelic trap), and MAF is what tells them apart.

Class: `sv` if the record is in the sv scan (>= 50 bp in the panel's class split), else
`smallindel`. Type: deletion if REF is longer, insertion if ALT is; equal-length MNPs are kept
as `mnp` so nothing is silently dropped, and the notebook sets them aside explicitly.

Signed distance to the nearest TSS (strand-aware; negative = upstream) comes from the same
gene models as the classifier (lib.load_genes).

Writes results/gea_hit_landscape.csv.gz. env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, GEA)
import lib                                                       # noqa: E402

LFMM = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile"
OUT = f"{HERE}/results"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1", "pc2", "pc3"]
ID = ["chrom", "pos", "ref_len", "alt_len", "mafr"]


def tested(cls, ax):
    d = pd.read_csv(f"{LFMM}/lfmm_{cls}_gen9_{ax}.csv")
    d = d[(d.MAF > 0.05) & d.block.notna() & (d.block != "")].copy()
    d["mafr"] = d.MAF.round(6)
    d["sig"] = d.pval < 0.05 / len(d)
    return d


def tss_distance(U):
    """Signed distance from the variant start to the nearest TSS, in gene orientation."""
    G = lib.load_genes().drop_duplicates("gene")
    G["tss"] = np.where(G.strand == "+", G.start, G.end)
    out = np.full(len(U), np.nan); near = np.empty(len(U), dtype=object)
    for ch, g in G.groupby("chrom"):
        m = (U.chrom == ch).values
        if not m.any():
            continue
        g = g.sort_values("tss"); t = g.tss.values; s = g.strand.values; n = g.gene.values
        p = U.pos.values[m]
        i = np.clip(np.searchsorted(t, p), 1, len(t) - 1)
        j = np.where(np.abs(p - t[i - 1]) <= np.abs(p - t[i]), i - 1, i)
        d = p - t[j]
        out[m] = np.where(s[j] == "+", d, -d)      # negative = upstream of the TSS
        near[m] = n[j]
    return out, near


def main():
    U = tested("nonsnp", AXES[0])[["chrom", "pos", "ref_len", "alt_len", "MAF", "mafr", "block"]]
    sv = tested("sv", AXES[0])[ID].drop_duplicates().assign(_sv=True)
    U = U.merge(sv, on=ID, how="left")
    U["cls"] = np.where(U._sv.fillna(False).astype(bool), "sv", "smallindel")
    U = U.drop(columns="_sv")
    U["kind"] = np.select([U.ref_len > U.alt_len, U.alt_len > U.ref_len],
                          ["deletion", "insertion"], "mnp")
    U["size"] = (U.ref_len - U.alt_len).abs()

    hits = {k: [] for k in map(tuple, U[ID].values)}
    minp = pd.Series(1.0, index=pd.MultiIndex.from_frame(U[ID]))
    for ax in AXES:
        for cls in ("nonsnp", "sv", "smallindel"):
            d = tested(cls, ax)
            s = d[d.sig]
            for k in map(tuple, s[ID].values):
                if k in hits and ax not in hits[k]:
                    hits[k].append(ax)
            if cls == "nonsnp":
                pm = d.groupby(ID).pval.min()
                minp = np.minimum(minp, pm.reindex(minp.index).fillna(1.0))
        print(f"  {ax}: running union {sum(bool(v) for v in hits.values()):,} hit records", flush=True)
    U["hit_axes"] = [",".join(hits[k]) for k in map(tuple, U[ID].values)]
    U["n_hit_axes"] = U.hit_axes.str.split(",").map(lambda x: len([a for a in x if a]))
    U["hit"] = U.n_hit_axes > 0
    U["min_p"] = minp.reindex(pd.MultiIndex.from_frame(U[ID])).values

    R = pd.read_csv(f"{OUT}/gea_universe_regions.csv.gz")
    U = U.merge(R, on=["chrom", "pos", "ref_len"], how="left")
    U["dist_tss"], U["nearest_tss_gene"] = tss_distance(U)

    U.to_csv(f"{OUT}/gea_hit_landscape.csv.gz", index=False)
    h = U[U.hit]
    print(f"\n{len(U):,} tested records; {int(U.hit.sum()):,} hits "
          f"({dict(h.cls.value_counts())}; types {dict(h.kind.value_counts())}); "
          f"{h.block.nunique():,} hit blocks")
    print(f"wrote {OUT}/gea_hit_landscape.csv.gz")


if __name__ == "__main__":
    main()

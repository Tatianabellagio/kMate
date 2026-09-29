#!/usr/bin/env python
"""Genomic context for every variant in `sv_purging_classify_variants{TAG}.csv`.

One class per variant, first match wins:  CDS > UTR > intron > TE > promoter > intergenic
  CDS / UTR   overlaps a TAIR10 CDS / 5' or 3' UTR
  intron      inside a gene but not in an exon
  TE          inside a reference transposable element (TAIR10 transposable_element)
  promoter    within 2 kb upstream of a gene TSS (strand-aware)
  intergenic  none of the above
Also records the overlapping gene, or the gene whose promoter it sits in -- reused for the candidate
list. Span = [pos, pos + ref_len) (an insertion is effectively the point `pos`).

Env: basic. Writes results/sv_adaptive/sv_purging_annotated.csv.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
PROM = 2000
TAG = os.environ.get("SV_TAG", "")


def merged(df):
    """non-overlapping sorted intervals per chrom -> {chrom: (starts, ends)}"""
    out = {}
    for c, d in df.groupby("chrom"):
        d = d.sort_values("start"); s, e = [], []
        for a, b in zip(d.start.to_numpy(), d.end.to_numpy()):
            if s and a <= e[-1]: e[-1] = max(e[-1], b)
            else: s.append(a); e.append(b)
        out[c] = (np.array(s), np.array(e))
    return out


def hit(iv, chrom, lo, hi):
    """bool: does [lo, hi] overlap any merged interval on chrom"""
    res = np.zeros(lo.size, bool)
    for c in np.unique(chrom):
        if c not in iv: continue
        m = chrom == c; s, e = iv[c]
        j = np.searchsorted(s, hi[m], side="right") - 1
        ok = j >= 0; r = np.zeros(m.sum(), bool)
        r[ok] = e[j[ok]] >= lo[m][ok]
        res[m] = r
    return res


def main():
    V = pd.read_csv(f"{OUT}/sv_purging_classify_variants{TAG}.csv")
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz")
    V["ref_len"] = idx["ref_len"].astype(np.int64)[V.afrow.to_numpy()]
    g = pd.read_csv(lib.TAIR10_GENES_TE, sep="\t", header=None, comment="#",
                    names=["chrom", "src", "feat", "start", "end", "score", "strand", "frame", "attr"])
    g = g[g.chrom.isin([f"Chr{i}" for i in range(1, 6)])]
    genes = g[g.feat == "gene"].copy(); genes["gid"] = genes.attr.str.extract(r"ID=([^;]+)")
    prom = genes.copy()
    prom["start"] = np.where(genes.strand == "+", genes.start - PROM, genes.end + 1)
    prom["end"] = np.where(genes.strand == "+", genes.start - 1, genes.end + PROM)
    IV = {"CDS": merged(g[g.feat == "CDS"]),
          "UTR": merged(g[g.feat.isin(["five_prime_UTR", "three_prime_UTR"])]),
          "exon": merged(g[g.feat == "exon"]),
          "gene": merged(genes),
          "TE": merged(g[g.feat == "transposable_element"]),
          "promoter": merged(prom)}
    ch = V.chrom.to_numpy(); lo = V.pos.to_numpy(); hi = lo + np.maximum(V.ref_len.to_numpy() - 1, 0)
    H = {k: hit(v, ch, lo, hi) for k, v in IV.items()}
    cat = np.full(len(V), "intergenic", dtype=object)
    for k, cond in (("promoter", H["promoter"]), ("TE", H["TE"]),
                    ("intron", H["gene"] & ~H["exon"]), ("UTR", H["UTR"]), ("CDS", H["CDS"])):
        cat[cond] = k                                    # applied low->high so the highest priority wins
    V["context"] = cat

    # gene ID: overlapping gene, else the gene whose promoter window contains it
    gid = np.full(len(V), "", dtype=object)
    for c in np.unique(ch):
        m = np.where(ch == c)[0]
        for src in (genes[genes.chrom == c], prom[prom.chrom == c]):
            s = src.start.to_numpy(); e = src.end.to_numpy(); ids = src.gid.to_numpy()
            o = np.argsort(s); s, e, ids = s[o], e[o], ids[o]
            for k in m:
                if gid[k]: continue
                j = np.searchsorted(s, hi[k], side="right")
                cand = np.where(e[max(0, j - 200):j] >= lo[k])[0]
                if cand.size: gid[k] = ids[max(0, j - 200) + cand[-1]]
    V["gene"] = gid
    V.to_csv(f"{OUT}/sv_purging_annotated{TAG}.csv", index=False)
    print(pd.crosstab([V.cls, V.origin], V.context).loc[["insertion", "deletion", "indel"]]
          [["CDS", "UTR", "intron", "TE", "promoter", "intergenic"]].to_string())
    print(f"\n[wrote] {OUT}/sv_purging_annotated{TAG}.csv")


if __name__ == "__main__":
    main()

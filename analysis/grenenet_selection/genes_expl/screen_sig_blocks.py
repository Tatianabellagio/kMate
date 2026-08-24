#!/usr/bin/env python
"""Functional-context screen of every significant SV / small-indel / non-SNP GEA
block-lead variant, judged AT THE VARIANT'S OWN POSITION (not the block's
attributed genes).

Motivation (the CARK lesson)
----------------------------
A clq0.9 tiling+merge block can be attributed to a gene its lead variant does not
touch: for the sparse SV class, `merge_small_blocks` chains empty blocks backward,
so the WZA unit balloons to tens of kb and the causal-candidate SV can sit in a
gene desert far from the block's overlapping genes (see genes_expl/README.md,
CARK8/9 <- a 5.1 kb insertion ~29 kb away, r²≤0.05 to CARK). A distal, unlinked SV
has no obvious functional handle. So to find REAL candidates we ask, per
significant lead variant: does IT physically fall in a gene body (CDS / UTR /
intron), a promoter (≤1 kb upstream, strand-aware), a TE, or intergenic space?

Significance set
----------------
Union over cls in {sv, smallindel, nonsnp} x axis in {bio1..19, pc1} of the
per-class genome-wide Bonferroni (0.05/n, MAF>0.05) block-LEAD hits from the raw
(pre-WZA) LFMM p in wza_in_clq09_tile. NB: raw uncalibrated p (matches the
newpeak dotgrid); treat as a candidate net, not a calibrated hit list.

Outputs (genes_expl/):
  sig_block_functional_screen.csv   one row per unique lead variant: position,
     class(es), best axis/nlp, n_axes, MAF, genomic context (tier + gene +
     sub-region + strand + distance), promoter/TE flags, block tagging-artifact
     flag (block-span genes != variant gene), + gene annotation for genic/promoter.
env: kmate.  Run on a compute node (GFF parse + optional API annotation).
"""
from __future__ import annotations
import os, sys, importlib.util as ilu
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import lib
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "phase1_replication", "multiaxis"))

WZAIN = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
CLASSES = ["sv", "smallindel", "nonsnp"]
GFF = lib.TAIR10_GENES_TE                      # full hierarchy incl. TEs
PROM_BP = 1000                                 # promoter window upstream of TSS
PROXIMAL_BP = 5000                             # "near a gene" cutoff


# ------------------------------------------------------------------ significance
def sig_leads():
    rows = []
    for cls in CLASSES:
        for ax in AXES:
            f = f"{WZAIN}/lfmm_{cls}_gen9_{ax}.csv"
            if not os.path.exists(f):
                continue
            d = pd.read_csv(f)
            d = d[(d.MAF > 0.05) & d.block.notna() & (d.block != "")].copy()
            if not len(d):
                continue
            d["nlp"] = -np.log10(d.pval.clip(lower=1e-300))
            bonf = -np.log10(0.05 / len(d))
            L = d.loc[d.groupby("block")["nlp"].idxmax()]      # block lead
            sig = L[L.nlp > bonf]
            for _, r in sig.iterrows():
                rows.append((cls, ax, r.block, r.chrom, int(r.pos), int(r.ref_len),
                             int(r.alt_len), float(r.MAF), float(r.nlp)))
    C = pd.DataFrame(rows, columns=["cls", "axis", "block", "chrom", "pos",
                                    "ref_len", "alt_len", "MAF", "nlp"])
    # dedup to unique variant
    key = ["chrom", "pos", "ref_len", "alt_len"]
    g = C.sort_values("nlp", ascending=False).groupby(key)
    out = g.first().reset_index()[key + ["MAF"]]
    out["size"] = (out.alt_len - out.ref_len).abs()
    out["vclass"] = np.where(out["size"] == 0, "snp",
                             np.where(out["size"] > 50, "sv", "smallindel"))
    out["best_axis"] = g["axis"].first().values
    out["best_nlp"] = g["nlp"].max().values
    out["n_axes"] = g["axis"].nunique().values
    out["sig_classes"] = g["cls"].apply(lambda s: ",".join(sorted(set(s)))).values
    out["block"] = g["block"].first().values                  # block of strongest hit
    return out.sort_values("best_nlp", ascending=False).reset_index(drop=True)


# ------------------------------------------------------------------ GFF parse
def load_gff():
    genes, cds, u5, u3, exon, te = [], [], [], [], [], []
    with open(GFF) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            ch, typ, s, e, strand, attr = p[0], p[2], int(p[3]), int(p[4]), p[6], p[8]
            if ch not in lib.CHROMS:
                continue
            if typ == "gene":
                gid = attr.split("ID=")[1].split(";")[0]
                genes.append((ch, s, e, strand, gid))
            elif typ in ("CDS", "exon", "five_prime_UTR", "three_prime_UTR"):
                par = attr.split("Parent=")[1].split(";")[0].split(",")[0]
                gid = par.split(".")[0]
                tgt = {"CDS": cds, "exon": exon, "five_prime_UTR": u5,
                       "three_prime_UTR": u3}[typ]
                tgt.append((ch, s, e, gid))
            elif "transposable_element" in typ:
                te.append((ch, s, e))
    G = pd.DataFrame(genes, columns=["chrom", "start", "end", "strand", "gene"])
    sub = {k: pd.DataFrame(v, columns=["chrom", "start", "end", "gene"])
           for k, v in (("cds", cds), ("u5", u5), ("u3", u3), ("exon", exon))}
    TE = pd.DataFrame(te, columns=["chrom", "start", "end"])
    return G, sub, TE


def _feat_index(df):
    """per-chrom dict of (start[], end[], gene[]?) sorted by start, for overlap."""
    out = {}
    for ch, d in df.groupby("chrom"):
        d = d.sort_values("start")
        out[ch] = (d.start.to_numpy(), d.end.to_numpy(),
                   d.gene.to_numpy() if "gene" in d else None)
    return out


def overlaps_gene_of(idx, ch, gid, vs, ve):
    """True if [vs,ve] overlaps any interval in idx[ch] belonging to gene gid."""
    if ch not in idx:
        return False
    st, en, gn = idx[ch]
    m = (gn == gid) & (st <= ve) & (en >= vs)
    return bool(m.any())


# ------------------------------------------------------------------ classify
def classify(C):
    G, sub, TE = load_gff()
    # per-chrom numpy arrays for genes (sorted by start) + feature/gene overlap sets
    gnp = {}
    for ch, d in G.groupby("chrom"):
        d = d.sort_values("start")
        gnp[ch] = (d.start.to_numpy(np.int64), d.end.to_numpy(np.int64),
                   d.strand.to_numpy(str), d.gene.to_numpy(str))
    cds_i, u5_i, u3_i, exon_i = (_feat_index(sub[k]) for k in ("cds", "u5", "u3", "exon"))
    tenp = {ch: (d.start.to_numpy(np.int64), d.end.to_numpy(np.int64))
            for ch, d in TE.groupby("chrom")}

    tiers, genes, subreg, strands, dists, te_flags, nearest = [], [], [], [], [], [], []
    for _, r in C.iterrows():
        ch, vs, ve = r.chrom, int(r.pos), int(r.pos) + int(r.ref_len) - 1
        gene = subr = strand = ""
        tier = None
        st = en = sd = gn = None
        if ch in gnp:
            st, en, sd, gn = gnp[ch]
            body = (st <= ve) & (en >= vs)
            if body.any():
                k = int(np.where(body)[0][0])
                gene, strand = gn[k], sd[k]
                if overlaps_gene_of(cds_i, ch, gene, vs, ve):
                    tier, subr = "1_CDS", "CDS"
                elif overlaps_gene_of(u5_i, ch, gene, vs, ve):
                    tier, subr = "2_UTR", "5'UTR"
                elif overlaps_gene_of(u3_i, ch, gene, vs, ve):
                    tier, subr = "2_UTR", "3'UTR"
                elif overlaps_gene_of(exon_i, ch, gene, vs, ve):
                    tier, subr = "4_exon_noncoding", "exon(noncoding)"
                else:
                    tier, subr = "5_intron", "intron"
            elif True:                                          # promoter (vectorized)
                lo = np.where(sd == "+", st - PROM_BP, en + 1)
                hi = np.where(sd == "+", st - 1, en + PROM_BP)
                pm = (vs <= hi) & (ve >= lo)
                if pm.any():
                    k = int(np.where(pm)[0][0])
                    gene, strand, tier, subr = gn[k], sd[k], "3_promoter", "promoter"
        # TE overlap
        te_hit = False
        if ch in tenp:
            tst, ten = tenp[ch]
            te_hit = bool(((tst <= ve) & (ten >= vs)).any())
        # nearest gene + distance
        ngene, ndist = "", np.nan
        if st is not None:
            d0 = np.maximum.reduce([st - ve, vs - en, np.zeros(len(st), np.int64)])
            j = int(np.argmin(d0)); ngene, ndist = gn[j], int(d0[j])
        if tier is None:
            if te_hit:
                tier, subr = "6_TE", "TE"
            elif not np.isnan(ndist) and ndist <= PROXIMAL_BP:
                tier, subr = "7_proximal_intergenic", "intergenic"
            else:
                tier, subr = "8_gene_desert", "intergenic"
        tiers.append(tier); genes.append(gene); subreg.append(subr)
        strands.append(strand); te_flags.append(te_hit)
        nearest.append(ngene); dists.append(ndist)
    C = C.copy()
    C["tier"] = tiers; C["region"] = subreg; C["gene"] = genes
    C["strand"] = strands; C["te_overlap"] = te_flags
    C["nearest_gene"] = nearest; C["dist_to_gene"] = dists
    return C, G


# ------------------------------------------------------------ block tagging flag
def block_span_genes(C, G):
    """Genes overlapping each block's TILING span (the dotgrid attribution) -> flag
    when the variant's own gene is NOT among them (CARK-type merge artifact)."""
    import blocks_tiling as bt
    spans = {}
    for ch in lib.CHROMS:
        b = bt.load_blocks(0.9, ch)
        ends = b.end_pos.to_numpy(np.int64)
        for i in range(len(b)):
            lo = 1 if i == 0 else int(ends[i - 1]) + 1
            hi = int(ends[i]) if i < len(b) - 1 else 10**9
            spans[f"{ch}_{i}"] = (ch, lo, hi)
    battr, flag = [], []
    for _, r in C.iterrows():
        sp = spans.get(r.block)
        if sp is None:
            battr.append(""); flag.append(False); continue
        ch, lo, hi = sp
        gs = G[(G.chrom == ch) & (G.end >= lo) & (G.start <= hi)].gene.tolist()
        battr.append(";".join(gs))
        # artifact: block attributes gene(s) but the variant's own gene is not among
        # them (CARK-type — the lead variant sits outside the block's genes)
        flag.append(bool(gs) and (r.gene not in gs))
    C = C.copy()
    C["block_span_genes"] = battr
    C["block_gene_mismatch"] = flag
    return C


# ------------------------------------------------------------ gene annotation
def annotate(C):
    annp = os.path.join(os.path.dirname(HERE), "phase1_replication",
                        "annotate_genes_tair_uniprot.py")
    spec = ilu.spec_from_file_location("ann_tu", annp)
    mod = ilu.module_from_spec(spec); spec.loader.exec_module(mod)
    want = sorted(set(C.loc[C.tier.isin(["1_CDS", "2_UTR", "3_promoter",
                                         "4_exon_noncoding", "5_intron"]), "gene"])
                  - {""})
    print(f"annotating {len(want)} candidate genes ...", flush=True)
    A = mod.annotate(want).set_index("gene")
    for col in ("symbol", "protein_name", "categories"):
        C[col] = C["gene"].map(A[col]) if col in A else ""
    return C


def main():
    C = sig_leads()
    print(f"unique significant lead variants: {len(C)} "
          f"({(C.vclass=='sv').sum()} sv, {(C.vclass=='smallindel').sum()} indel)")
    C, G = classify(C)
    C = block_span_genes(C, G)
    C = annotate(C)
    C = C.sort_values(["tier", "best_nlp"], ascending=[True, False])
    cols = ["chrom", "pos", "ref_len", "alt_len", "size", "vclass", "sig_classes",
            "best_axis", "best_nlp", "n_axes", "MAF", "tier", "region", "gene",
            "symbol", "protein_name", "categories", "strand", "te_overlap",
            "nearest_gene", "dist_to_gene", "block", "block_span_genes",
            "block_gene_mismatch"]
    out = f"{HERE}/sig_block_functional_screen.csv"
    C[cols].to_csv(out, index=False)
    print(f"wrote {out}")
    print("\ntier counts (all classes):")
    print(C.tier.value_counts().sort_index().to_string())
    print("\ntier x SV only:")
    print(C[C.vclass == "sv"].tier.value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()

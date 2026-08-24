#!/usr/bin/env python
"""Genes under blocks the NON-SNP scan flags that the SNP scan misses -- the candidate
kMate-unique (indel+SV) signals, across every contrast we ran.

For each contrast, take the clq0.9 blocks BH-FDR-significant in the non-SNP GWAS but NOT in the
SNP GWAS (block-level, same convention as class_gwas_multitrait.ipynb / class_gwas_persite.ipynb),
across:
  - multitrait JOINT, GLOBAL
  - multitrait CLIMATE x (bio1..19 + PC1 of all bioclim)
  - per-site single-site scans (31 gardens)
Union the non-SNP-only blocks, map each to its clq0.9 interval -> overlapping TAIR10 genes.

Writes analysis/grenenet_selection/r3_persite_gwas/results/varexp/nonsnp_only_blocks.csv (one row per block) and
nonsnp_only_genes.csv (one row per gene; UniProt descriptions filled in by
nonsnp_only_genes_describe.py). Env: kmate.
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

OUT = f"{lib.GEA}/r3_persite_gwas/results/varexp"
CLASSES = ["snp", "nonsnp"]
FLANK = 2000          # promoter/regulatory flank (bp) each side of the block span for gene overlap


def load_class(name):
    d = np.load(f"{OUT}/class_gwas_{name}.npz", allow_pickle=True)
    return dict(chrom=d["chrom"], pos=d["pos"], Z=d["Z"], sites=d["sites"], bio1=d["bio1"],
                p_joint=d["p_joint"], p_global=d["p_global"])


def block_ids(chrom_lower, pos):
    chrom_cap = np.array([c.replace("chr", "Chr") for c in chrom_lower])
    return lib.assign_clq_blocks(chrom_cap, pos, r2=0.9)


def sig_blocks(blk, p, q=0.05):
    """(fdr_set, bonf_set) of clq0.9 blocks significant at BH-FDR q<0.05 and Bonferroni 0.05,
    block-level (lead = min-p marker per block), matching the notebook overlap tables."""
    keep = blk != ""
    nlp = -np.log10(np.clip(p[keep], 1e-300, 1))
    df = pd.DataFrame({"block": blk[keep], "nlp": nlp})
    lead = lib.collapse_to_blocks(df, stat="nlp", block_col="block")
    p_lead = 10.0 ** (-lead["nlp"].to_numpy())
    m = len(p_lead)
    qv = lib.bh(p_lead)
    fdr = set(lead.loc[qv < q, "block"])
    bonf = set(lead.loc[p_lead < 0.05 / m, "block"])
    return fdr, bonf


def climate_p(Z, C, bvec):
    S = Z.shape[1]
    Cinv = np.linalg.pinv(C); one = np.ones(S)
    dg = float(one @ Cinv @ one)
    c0 = (bvec - bvec.mean()) / bvec.std()
    c = c0 - (float(one @ Cinv @ c0) / dg) * one
    dc = float(c @ Cinv @ c)
    z = (Z @ Cinv @ c) / np.sqrt(dc)
    return 2 * stats.norm.sf(np.abs(z))


def block_interval_map():
    """block id 'Chr{n}_{idx}' -> (chrom, start_pos, end_pos). idx = rank by start_pos within
    chrom, matching lib.assign_clq_blocks."""
    m = {}
    for ci in range(1, 6):
        c = f"Chr{ci}"
        f = f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_clq0.9_blocks_clq0.9.tsv"
        g = pd.read_csv(f, sep="\t").sort_values("start_pos").reset_index(drop=True)
        for idx, r in g.iterrows():
            m[f"{c}_{idx}"] = (c, int(r.start_pos), int(r.end_pos))
    return m


def main():
    data = {c: load_class(c) for c in CLASSES}
    sites = data["snp"]["sites"]
    blk = {c: block_ids(data[c]["chrom"], data[c]["pos"]) for c in CLASSES}
    Cmat = {c: np.corrcoef(data[c]["Z"].T) for c in CLASSES}

    clim = lib.load_climate().reindex(sites)[lib.BIO_COLS]
    Bz = (clim - clim.mean()) / clim.std()
    Uu, Ss, _ = np.linalg.svd(Bz.to_numpy(), full_matrices=False)
    axes = {f"bio{i}": clim[f"bio{i}"].to_numpy(float) for i in range(1, 20)}
    axes["PC1_allbio"] = Uu[:, 0] * Ss[0]

    # contrast -> {class -> p-vector}
    contrasts = {}
    for name, key in [("JOINT", "p_joint"), ("GLOBAL", "p_global")]:
        contrasts[f"multitrait_{name}"] = {c: data[c][key] for c in CLASSES}
    for ax, bvec in axes.items():
        contrasts[f"multitrait_CLIMATE_{ax}"] = {c: climate_p(data[c]["Z"], Cmat[c], bvec) for c in CLASSES}
    for si, sid in enumerate(sites):
        contrasts[f"persite_site{int(sid)}"] = {
            c: 2 * stats.norm.sf(np.abs(data[c]["Z"][:, si])) for c in CLASSES}

    # non-SNP-only blocks per contrast, at BOTH thresholds:
    #   FDR  : block nonsnp-FDR-sig  AND NOT snp-FDR-sig  in that contrast
    #   Bonf : block nonsnp-Bonf-sig AND NOT snp-Bonf-sig in that contrast (stricter subset)
    block_fdr = {}    # block -> set(contrasts) where non-SNP-only at FDR
    block_bonf = {}   # block -> set(contrasts) where non-SNP-only at Bonferroni
    for cname, pv in contrasts.items():
        snp_fdr, snp_bonf = sig_blocks(blk["snp"], pv["snp"])
        ns_fdr, ns_bonf = sig_blocks(blk["nonsnp"], pv["nonsnp"])
        for b in (ns_fdr - snp_fdr):
            block_fdr.setdefault(b, set()).add(cname)
        for b in (ns_bonf - snp_bonf):
            block_bonf.setdefault(b, set()).add(cname)
    print(f"{len(contrasts)} contrasts scanned; {len(block_fdr)} non-SNP-only blocks at FDR, "
          f"{len(block_bonf)} at Bonferroni", flush=True)

    # map blocks -> intervals -> genes (in-block AND within +/-FLANK bp)
    bmap = block_interval_map()
    genes = lib.load_genes()
    rows = []
    for b in sorted(block_fdr):
        if b not in bmap:
            print(f"  WARN block {b} has no interval", flush=True); continue
        c, s, e = bmap[b]
        inb = genes[(genes.chrom == c) & (genes.start <= e) & (genes.end >= s)]
        flk = genes[(genes.chrom == c) & (genes.start <= e + FLANK) & (genes.end >= s - FLANK)]
        flk = flk[~flk.gene.isin(inb.gene)]                     # flank-only (not already in-block)
        cs = block_fdr[b]
        rows.append(dict(block=b, chrom=c, start=s, end=e, span_bp=e - s,
                         bonferroni=b in block_bonf,
                         n_contrasts=len(cs), contrasts=";".join(sorted(cs)),
                         contrasts_bonf=";".join(sorted(block_bonf.get(b, []))),
                         n_genes_in=len(inb), genes_in_block=";".join(inb.gene),
                         n_genes_flank=len(flk), genes_flank2kb=";".join(flk.gene),
                         gene_names_in=";".join(inb.name.fillna(""))))
    bdf = pd.DataFrame(rows).sort_values(["bonferroni", "n_contrasts", "block"],
                                         ascending=[False, False, True])
    bdf.to_csv(f"{OUT}/nonsnp_only_blocks.csv", index=False)
    print(f"wrote nonsnp_only_blocks.csv ({len(bdf)} blocks; "
          f"{bdf.bonferroni.sum()} Bonferroni; {(bdf.n_genes_in+bdf.n_genes_flank).sum()} gene-hits incl flank)",
          flush=True)

    # explode to unique genes (in-block + flank), tagging overlap_type & threshold
    grows = {}
    for _, r in bdf.iterrows():
        for gid, otype in ([(g, "in_block") for g in (r.genes_in_block.split(";") if r.genes_in_block else [])]
                           + [(g, "flank2kb") for g in (r.genes_flank2kb.split(";") if r.genes_flank2kb else [])]):
            if gid not in grows:
                grows[gid] = dict(gene=gid, chrom=r.chrom, blocks=set(), contrasts=set(),
                                  overlap="flank2kb", bonferroni=False)
            grows[gid]["blocks"].add(r.block)
            grows[gid]["contrasts"].update(r.contrasts.split(";"))
            if otype == "in_block":
                grows[gid]["overlap"] = "in_block"     # in_block wins over flank2kb
            if r.bonferroni:
                grows[gid]["bonferroni"] = True
    gout = [dict(gene=g, chrom=d["chrom"], overlap_type=d["overlap"], bonferroni=d["bonferroni"],
                 n_blocks=len(d["blocks"]), blocks=";".join(sorted(d["blocks"])),
                 n_contrasts=len(d["contrasts"]), contrasts=";".join(sorted(d["contrasts"])))
            for g, d in grows.items()]
    gdf = pd.DataFrame(gout).sort_values(["bonferroni", "n_contrasts", "gene"],
                                         ascending=[False, False, True])
    gdf.to_csv(f"{OUT}/nonsnp_only_genes.csv", index=False)
    print(f"wrote nonsnp_only_genes.csv ({len(gdf)} unique genes; "
          f"{(gdf.overlap_type=='in_block').sum()} in-block, {(gdf.overlap_type=='flank2kb').sum()} flank-only)",
          flush=True)


if __name__ == "__main__":
    main()

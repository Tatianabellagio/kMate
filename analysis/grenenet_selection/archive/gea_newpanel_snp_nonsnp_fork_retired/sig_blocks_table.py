#!/usr/bin/env python
"""Site-level significant-BLOCK table across all climate axes (bio1..19 + pc1),
SNP vs non-SNP, with the SNP<->nonSNP block overlap per axis.

Corrections applied (per this session's findings):
  * site-level MAF >= 0.05 filter (computed on the 31 site means), BEFORE BH/Bonferroni
  * deflate-only 'fair' p: use GIF-calibrated p only when GIF>=1, else raw p (never let
    an under-dispersed axis with GIF<1 inflate significance)
  * significant records -> clq0.9 haploblocks (blocks_clq09); a block is 'significant'
    if it contains >=1 Bonferroni-significant record.

Writes:
  power/sig_blocks_by_axis.csv     cls,axis,GIF,credible,n_bonf_rec,n_bonf_blocks,n_fdr05_blocks,min_p
  power/sig_blocks_overlap.csv     axis,snp_blocks,nonsnp_blocks,shared,shared_ids
  power/sig_blocks_detail.csv      every Bonferroni-sig record (cls,axis,chrom,pos,block,...,gene)
  lfmm_{cls}_sitepc1_clq09.csv     per-record site PC1 fair p + clq0.9 block (for the notebook)
"""
from __future__ import annotations
import os, sys, glob, re
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blocks_clq09
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

POW = f"{lib.GEA}/gea_newpanel/results/power"
SITE = f"{lib.GEA}/gea_newpanel/results/lfmm_site"
GNP = f"{lib.GEA}/gea_newpanel/results"
CM = f"{lib.GEA}/phase1_replication/results/class_matrices"
_P0 = {"snp": "p0_snp.npy", "nonsnp": "p0_nonsnp.npy", "sv": "p0_nonsnp.npy"}
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
MAF_MIN = 0.05


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = (p[o] * n) / (np.arange(n) + 1); q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)


def gifval(f):
    g = f.replace(".csv", ".gif.txt"); return float(open(g).read().strip())


def main():
    genes = lib.load_genes()
    by_axis, detail, overlap = [], [], []
    sig_blocks = {}  # (axis) -> {cls: set(blocks)}
    for cls in ["snp", "nonsnp", "sv"]:
        recs = pd.read_csv(f"{CM}/{cls}_gen9.records.csv")
        dims = np.loadtxt(f"{SITE}/lfmm_{cls}_site_dims.txt", dtype=int)
        nsite, nrec = int(dims[0]), int(dims[1])
        dp = np.fromfile(f"{SITE}/lfmm_{cls}_site_Y.f64").reshape(nsite, nrec)
        p0 = np.load(f"{lib.AF_STORE}/{_P0[cls]}")[recs["col"].to_numpy()]
        site_af = dp + p0[None, :]
        site_maf = np.minimum(site_af.mean(0), 1 - site_af.mean(0))
        block = blocks_clq09.assign_clq09_blocks(recs.chrom.to_numpy(), recs.pos.to_numpy())
        recs = recs.assign(site_maf=site_maf, block=block)
        keep = site_maf >= MAF_MIN
        print(f"{cls}: {nrec:,} records, site-MAF>=%.2f -> {int(keep.sum()):,}" % MAF_MIN, flush=True)

        for ax in AXES:
            f = f"{POW}/exp2_scan_{cls}_{ax}_bothp.csv"
            if not os.path.exists(f):
                continue
            g = gifval(f); d = pd.read_csv(f)
            fair = d.pval_gif.values   # pure GIF-calibrated (no deflate-only; inflated axes shown as-is)
            rk = recs.assign(p=fair)[keep].copy()
            n = len(rk); q = bh(rk.p.values); rk["q"] = q
            bonf_thr = 0.05 / n
            sig = rk[rk.p < bonf_thr]
            fdr = rk[rk.q < 0.05]
            fdr10 = rk[rk.q < 0.10]
            sb = set(sig.block[sig.block != ""])
            fb = set(fdr.block[fdr.block != ""])
            sig_blocks.setdefault(ax, {})[cls] = {"bonf": sb, "fdr": fb}
            by_axis.append(dict(cls=cls, axis=ax, GIF=round(g, 2),
                credible=bool(0.85 <= g <= 1.25), n_bonf_rec=len(sig),
                n_bonf_blocks=len(sb), n_fdr05_blocks=len(fb),
                n_fdr05_rec=len(fdr), n_fdr10_blocks=fdr10.block[fdr10.block != ""].nunique(),
                min_p=float(rk.p.min())))
            # detail: keep every FDR.05-significant record (superset of Bonferroni)
            if len(fdr):
                a = lib.annotate_svs(fdr.copy(), flank=2000, genes=genes)
                for (_, r), gn in zip(fdr.iterrows(), a.gene_name.values):
                    detail.append(dict(cls=cls, axis=ax, GIF=round(g, 2), chrom=r.chrom,
                        pos=int(r.pos), block=r.block, sv_size=int(abs(r.alt_len - r.ref_len)),
                        site_maf=round(float(r.site_maf), 3), p=float(r.p), q=float(r.q),
                        bonferroni=bool(r.p < bonf_thr), gene=gn or ""))
        # write per-record PC1 file for the notebook Manhattan (fair p + clq0.9 block)
        fpc = f"{POW}/exp2_scan_{cls}_pc1_bothp.csv"
        if os.path.exists(fpc):
            g = gifval(fpc); d = pd.read_csv(fpc)
            fair = d.pval_gif.values if g >= 1 else d.pval_raw.values
            out = recs.assign(pval=fair)[["chrom", "pos", "ref_len", "alt_len",
                                          "site_maf", "block", "pval"]].rename(columns={"site_maf": "MAF"})
            out.to_csv(f"{GNP}/lfmm_{cls}_sitepc1_clq09.csv", index=False)

    A = pd.DataFrame(by_axis)
    A.to_csv(f"{POW}/sig_blocks_by_axis.csv", index=False)
    # overlap per axis: pairwise SNP/nonSNP/SV block sharing, both thresholds
    def sh(a, b):
        return len(a & b)
    for ax, dd in sig_blocks.items():
        g = lambda c, t: dd.get(c, {}).get(t, set())
        row = dict(axis=ax)
        for t in ("bonf", "fdr"):
            s, n, v = g("snp", t), g("nonsnp", t), g("sv", t)
            row.update({f"snp_{t}": len(s), f"nonsnp_{t}": len(n), f"sv_{t}": len(v),
                        f"snp_nonsnp_{t}": sh(s, n), f"snp_sv_{t}": sh(s, v),
                        f"nonsnp_sv_{t}": sh(n, v)})
        row["sv_fdr_ids"] = ";".join(sorted(g("sv", "fdr")))
        overlap.append(row)
    O = pd.DataFrame(overlap)
    O.to_csv(f"{POW}/sig_blocks_overlap.csv", index=False)
    D = pd.DataFrame(detail); D.to_csv(f"{POW}/sig_blocks_detail.csv", index=False)

    pd.set_option("display.width", 240, "display.max_colwidth", 40)
    for thr, col, t in [("Bonferroni", "n_bonf_blocks", "bonf"),
                        ("FDR q<0.05", "n_fdr05_blocks", "fdr")]:
        print(f"\n===== significant BLOCKS by axis ({thr}, site-MAF>=0.05, GIF-aware) =====")
        print("  columns: snp/nonsnp/sv = #sig blocks; snp∩sv, nonsnp∩sv, snp∩nonsnp = shared blocks")
        piv = A.pivot(index="axis", columns="cls", values=col).fillna(0).astype(int)
        piv = piv.join(O.set_index("axis")[[f"snp_nonsnp_{t}", f"snp_sv_{t}", f"nonsnp_sv_{t}"]])
        gifp = A.pivot(index="axis", columns="cls", values="GIF")
        piv = piv.join(gifp["snp"].round(2).rename("GIF"))   # GIF shown, not used to filter
        print(piv.sort_values(["sv", "nonsnp", "snp"], ascending=False).to_string())
    print("\n===== significant records (FDR q<0.05; bonferroni flag) =====")
    if len(D):
        print(D.sort_values(["cls", "p"]).to_string(index=False))
    print(f"\nwrote {POW}/sig_blocks_by_axis.csv, sig_blocks_overlap.csv, sig_blocks_detail.csv")


if __name__ == "__main__":
    main()

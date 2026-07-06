#!/usr/bin/env python
"""Post-process LF-rank GEA: FDR, SNP-vs-nonSNP block concordance, gene annotation,
and save figure inputs. Run in kmate env (uses lib.annotate_svs / lib.load_genes)."""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/gea_newpanel/lf_rank"


def bh_fdr(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    m = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    q_ranked = ranked * m / (np.arange(1, m + 1))
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    qv = np.empty(m); qv[order] = np.clip(q_ranked, 0, 1)
    q[ok] = qv
    return q


def main():
    genes = lib.load_genes()
    res = {}
    for cls in ["snp", "nonsnp"]:
        d = pd.read_csv(f"{OUT}/lf_rank_{cls}_gen9_bio1.csv")
        d["q"] = bh_fdr(d["pval"].to_numpy())
        d["nlp"] = -np.log10(d["pval"].clip(lower=1e-300))
        m = int(np.isfinite(d["pval"]).sum())
        bonf = 0.05 / m
        d["bonf_hit"] = d["pval"] < bonf
        res[cls] = d
        n_bonf = int(d["bonf_hit"].sum())
        n_fdr = int((d["q"] < 0.05).sum())
        n_fdr01 = int((d["q"] < 0.01).sum())
        print(f"{cls}: tested={m:,} | Bonferroni({bonf:.1e}) hits={n_bonf} | "
              f"FDR q<0.05={n_fdr:,} | q<0.01={n_fdr01:,}")

    # ---- block-peak concordance ----
    def block_peak(d):
        g = d.groupby("block").agg(minp=("pval", "min"),
                                   chrom=("chrom", "first")).reset_index()
        # peak record position per block
        idx = d.groupby("block")["pval"].idxmin()
        pk = d.loc[idx, ["block", "chrom", "pos", "pval", "q", "MAF"]].rename(
            columns={"chrom": "chrom", "pos": "peak_pos", "pval": "peak_p", "q": "peak_q"})
        g = g.merge(pk[["block", "peak_pos", "peak_p", "peak_q"]], on="block")
        g["nlp"] = -np.log10(g["minp"].clip(lower=1e-300))
        return g

    bs, bn = block_peak(res["snp"]), block_peak(res["nonsnp"])
    mrg = bs.merge(bn, on="block", suffixes=("_snp", "_non"))
    from scipy.stats import spearmanr, pearsonr
    rho = spearmanr(mrg["nlp_snp"], mrg["nlp_non"]).correlation
    pear = pearsonr(mrg["nlp_snp"], mrg["nlp_non"])[0]
    print(f"\nBlock-peak concordance (n_blocks_shared={len(mrg):,}): "
          f"Spearman rho(-log10 p) = {rho:.3f} | Pearson = {pear:.3f}")

    # ---- class-specific top blocks (strong in one, weak in other) ----
    mrg["snp_only"] = mrg["nlp_snp"] - mrg["nlp_non"]
    mrg["non_only"] = mrg["nlp_non"] - mrg["nlp_snp"]
    thr_weak = -np.log10(0.05)
    snp_spec = mrg[(mrg["peak_p_snp"] < 1e-5) & (mrg["nlp_non"] < thr_weak)] \
        .sort_values("snp_only", ascending=False).head(15).copy()
    non_spec = mrg[(mrg["peak_p_non"] < 1e-5) & (mrg["nlp_snp"] < thr_weak)] \
        .sort_values("non_only", ascending=False).head(15).copy()

    def annot(df, chromcol, poscol):
        a = df.rename(columns={chromcol: "chrom", poscol: "pos"})[["chrom", "pos"]].copy()
        a = lib.annotate_svs(a, flank=2000, genes=genes)
        out = df.copy()
        out["gene_name"] = a["gene_name"].values
        out["genes_all"] = a["genes_all"].values
        return out

    snp_spec = annot(snp_spec, "chrom_snp", "peak_pos_snp")
    non_spec = annot(non_spec, "chrom_non", "peak_pos_non")

    # ---- top-25 hits per class, gene-annotated ----
    for cls in ["snp", "nonsnp"]:
        d = res[cls].sort_values("pval").head(25).copy()
        da = d.rename(columns={"pos": "pos"})[["chrom", "pos"]].copy()
        da = lib.annotate_svs(da, flank=2000, genes=genes)
        d["gene_name"] = da["gene_name"].values
        d["genes_all"] = da["genes_all"].values
        cols = ["chrom", "pos", "ref_len", "alt_len", "MAF", "block", "r",
                "pval", "pval_unadj", "q", "gene_name", "genes_all"]
        d[cols].to_csv(f"{OUT}/lf_rank_{cls}_gen9_bio1_top25_genes.csv", index=False)
        print(f"  top25 {cls} genes -> lf_rank_{cls}_gen9_bio1_top25_genes.csv "
              f"(lead {d.iloc[0]['chrom']}:{int(d.iloc[0]['pos'])} p={d.iloc[0]['pval']:.2e} "
              f"gene={d.iloc[0]['gene_name']})")

    scols = ["block", "chrom_snp", "peak_pos_snp", "peak_p_snp", "nlp_snp",
             "nlp_non", "gene_name", "genes_all"]
    ncols = ["block", "chrom_non", "peak_pos_non", "peak_p_non", "nlp_non",
             "nlp_snp", "gene_name", "genes_all"]
    snp_spec[scols].to_csv(f"{OUT}/lf_rank_SNPspecific_blocks.csv", index=False)
    non_spec[ncols].to_csv(f"{OUT}/lf_rank_NONSNPspecific_blocks.csv", index=False)
    print(f"\nSNP-specific blocks (strong SNP, weak nonSNP): top genes = "
          f"{list(snp_spec['gene_name'].head(8))}")
    print(f"nonSNP-specific blocks (strong nonSNP, weak SNP): top genes = "
          f"{list(non_spec['gene_name'].head(8))}")

    # ---- cross-check LFMM lead near Chr1:13,979,249 ----
    tgt = 13_979_249
    for cls in ["snp", "nonsnp"]:
        d = res[cls]
        near = d[(d["chrom"] == "Chr1") & (d["pos"].between(tgt - 20000, tgt + 20000))]
        if len(near):
            best = near.loc[near["pval"].idxmin()]
            print(f"  Chr1:{tgt} +/-20kb [{cls}]: n={len(near)} best "
                  f"pos={int(best['pos'])} p={best['pval']:.2e} q={best['q']:.2e} "
                  f"bonf_hit={bool(best['pval'] < 0.05/np.isfinite(d['pval']).sum())}")
        else:
            print(f"  Chr1:{tgt} +/-20kb [{cls}]: no records")

    # ---- save figure inputs ----
    np.savez_compressed(
        f"{OUT}/fig_inputs.npz",
        **{f"{c}_chrom": res[c]["chrom"].to_numpy().astype("U8") for c in res},
        **{f"{c}_pos": res[c]["pos"].to_numpy() for c in res},
        **{f"{c}_p": res[c]["pval"].to_numpy() for c in res},
        **{f"{c}_punadj": res[c]["pval_unadj"].to_numpy() for c in res},
        **{f"{c}_q": res[c]["q"].to_numpy() for c in res},
        block_rho=rho)
    # carry GIF numbers from the per-class npz
    gifs = {}
    for c in res:
        z = np.load(f"{OUT}/lf_rank_{c}_gen9_bio1.npz")
        gifs[c] = (float(z["gif_unadj"]), float(z["gif_adj"]),
                   z["gif_curve_k"], z["gif_curve_v"])
    np.savez_compressed(f"{OUT}/fig_gif.npz",
                        **{f"{c}_unadj": gifs[c][0] for c in gifs},
                        **{f"{c}_adj": gifs[c][1] for c in gifs},
                        **{f"{c}_curve_k": gifs[c][2] for c in gifs},
                        **{f"{c}_curve_v": gifs[c][3] for c in gifs})
    print("\nsaved fig_inputs.npz + fig_gif.npz")


if __name__ == "__main__":
    main()

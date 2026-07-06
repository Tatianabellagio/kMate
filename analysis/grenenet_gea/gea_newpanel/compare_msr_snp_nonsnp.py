#!/usr/bin/env python
"""SNP-vs-non-SNP comparison of the MSR structure-null rank GEA + gene annotation.

Reads msr_kendall_{cls}_site_bio1.csv (from run_msr_kendall.py):
  (i)   gene-annotated top-25 MSR hits per class
  (ii)  block-peak Spearman rho per class (shared LD-block ids) -> SAME vs DIFFERENT
  (iii) class-specific top blocks with genes
  (iv)  hit counts (p_msr<1e-4, FDR q<0.05) per class
  (v)   cross-check the site-LFMM lead near Chr1:13,979,249
Writes comparison.json + top-hit CSVs + block-merge CSV + manhattan npz for plotting.
"""
from __future__ import annotations
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/gea_newpanel/msr_kendall"
CLIM = "bio1"
MAF = 0.05
LFMM_LEAD = ("Chr1", 13_979_249)


def bh_fdr(p):
    p = np.asarray(p, float); n = p.size
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out


def main():
    genes = lib.load_genes()
    d = {}
    for cls in ["snp", "nonsnp"]:
        df = pd.read_csv(f"{OUT}/msr_kendall_{cls}_site_{CLIM}.csv")
        df["q_msr"] = bh_fdr(df.p_msr.to_numpy())
        d[cls] = df

    summary = {"climate": CLIM, "maf_min": MAF, "classes": {}}
    # ---- (i)+(iv) per-class hits + top-25 annotated ----
    for cls, df in d.items():
        m = df.site_maf >= MAF
        sub = df[m].copy()
        n_p1e4 = int((sub.p_msr < 1e-4).sum())
        n_fdr = int((sub.q_msr < 0.05).sum())
        sub["_absrho"] = sub.rho.abs()
        top = sub.sort_values(["p_msr", "_absrho"], ascending=[True, False]).head(25).copy()
        top = lib.annotate_svs(top, flank=2000, genes=genes)
        top.to_csv(f"{OUT}/top25_{cls}_{CLIM}.csv", index=False)
        summary["classes"][cls] = dict(
            M=int(len(df)), M_maf=int(m.sum()),
            n_p1e4=n_p1e4, n_fdr=n_fdr,
            min_p=float(sub.p_msr.min()), min_q=float(sub.q_msr.min()),
            max_absrho=float(sub.rho.abs().max()),
            top_genes=[g for g in top.gene_name.tolist() if g][:12])
        print(f"[{cls}] MAF>=.05 {int(m.sum()):,} | p<1e-4={n_p1e4} FDRq<.05={n_fdr} "
              f"| top genes: {[g for g in top.gene_name.tolist() if g][:8]}")

    # ---- (ii) block-peak rho per class, merged on shared LD-block ids ----
    def block_peak(df):
        s = df[df.site_maf >= MAF].copy()
        s["absrho"] = s.rho.abs()
        idx = s.groupby("block").absrho.idxmax()
        bp = s.loc[idx, ["block", "chrom", "pos", "rho", "absrho", "p_msr"]]
        return bp.set_index("block")
    bs, bn = block_peak(d["snp"]), block_peak(d["nonsnp"])
    shared = bs.join(bn, lsuffix="_snp", rsuffix="_non", how="inner")
    r_pear = float(np.corrcoef(shared.rho_snp, shared.rho_non)[0, 1])
    from scipy.stats import spearmanr
    r_spear = float(spearmanr(shared.rho_snp, shared.rho_non).statistic)
    shared.reset_index().to_csv(f"{OUT}/blockpeak_shared_{CLIM}.csv", index=False)
    summary["block_peak"] = dict(
        n_shared_blocks=int(len(shared)),
        pearson_rho_snp_vs_non=r_pear, spearman_snp_vs_non=r_spear)
    print(f"\n[block-peak] shared blocks={len(shared):,} | peak-rho SNP vs nonSNP: "
          f"Pearson={r_pear:.3f} Spearman={r_spear:.3f}")

    # ---- (iii) class-specific top blocks (strong in one, weak in other) ----
    shared["dominance"] = shared.absrho_non - shared.absrho_snp
    non_specific = shared.sort_values("dominance", ascending=False).head(15).reset_index()
    snp_specific = shared.sort_values("dominance").head(15).reset_index()
    for tagdf, tag in [(non_specific, "nonsnp_specific"), (snp_specific, "snp_specific")]:
        tagdf = tagdf.rename(columns={"pos_snp": "pos"})
        tagdf["ref_len"] = 1
        tagdf = lib.annotate_svs(tagdf.assign(chrom=tagdf.chrom_snp), flank=2000, genes=genes)
        tagdf.to_csv(f"{OUT}/blocks_{tag}_{CLIM}.csv", index=False)
        summary.setdefault("class_specific_blocks", {})[tag] = [
            dict(block=r.block, chrom=r.chrom_snp, pos=int(r.pos),
                 rho_snp=round(r.rho_snp, 3), rho_non=round(r.rho_non, 3),
                 gene=r.gene_name) for _, r in tagdf.head(8).iterrows()]

    # ---- (v) LFMM lead cross-check Chr1:13,979,249 ----
    lc, lp = LFMM_LEAD
    win = {}
    for cls, df in d.items():
        s = df[(df.chrom == lc) & (df.pos.between(lp - 25000, lp + 25000))
               & (df.site_maf >= MAF)]
        if len(s):
            best = s.loc[s.rho.abs().idxmax()]
            win[cls] = dict(n_in_window=int(len(s)),
                            best_pos=int(best.pos), best_rho=float(best.rho),
                            best_p_msr=float(best.p_msr),
                            best_q=float(bh_fdr(df.p_msr.to_numpy())[df.index.get_loc(best.name)])
                            if False else float(best.p_msr))
        else:
            win[cls] = dict(n_in_window=0)
    summary["lfmm_lead_crosscheck"] = dict(locus=f"{lc}:{lp}", window="+-25kb", **{f"{k}": v for k, v in win.items()})
    print(f"\n[LFMM lead {lc}:{lp} +-25kb] {json.dumps(win)}")

    # ---- manhattan npz (subsample non-significant for plotting size) ----
    for cls, df in d.items():
        s = df[df.site_maf >= MAF]
        pos_abs = s.pos.to_numpy()
        chrom = s.chrom.to_numpy().astype(str)
        p = s.p_msr.to_numpy()
        # keep all p<0.01 + a random 3% of the rest to bound figure size
        rng = np.random.default_rng(0)
        keep = (p < 0.01) | (rng.random(len(p)) < 0.03)
        np.savez(f"{OUT}/manhattan_{cls}_{CLIM}.npz",
                 chrom=chrom[keep], pos=pos_abs[keep], p_msr=p[keep],
                 p_naive=s.p_naive.to_numpy()[keep])

    with open(f"{OUT}/comparison.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n=== comparison.json ===")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

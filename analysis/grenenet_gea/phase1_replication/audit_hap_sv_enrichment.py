#!/usr/bin/env python
"""AUDIT of hap_sv_enrichment: data integrity + method correctness.

(A) data integrity: merge, NaN, has_sv sanity, direct SV-overlap spot-check on 3 blocks.
(B) representative choice: min-p (selection-biased) vs stat-independent (first-by-pos).
(C) the over-control problem: has_sv is collinear with block length (SVs live in long
    blocks) -> length is a MEDIATOR, not a clean confounder. Check common support and use a
    SIZE-MATCHED null (the project's established method) instead of a length-covariate logistic:
    match each SV block to the nearest SNP-only blocks by log10(len)+log10(nvar), compare
    %sig and median -log10p, permutation p by shuffling has_sv within length strata.
Env: kmate.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

CM = f"{lib.GEA}/phase1_replication/class_matrices"
KEN = f"{lib.GEA}/phase1_replication/kendall/kendall_hap_gen9_bio1.csv"
rng = np.random.default_rng(0)


def bh(p):
    p = np.asarray(p, float); m = len(p); o = p.argsort(); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def block_table(d, rep="minp"):
    """One row per block. rep='minp' (min-p hap, Sidak by n_hap) or 'firstpos'
    (stat-independent: the hap whose block_start is first, no within-block correction)."""
    g = d.groupby("block")
    if rep == "minp":
        lead = d.loc[g.pval.idxmin()].copy()
        lead["n_hap"] = g.size().reindex(lead.block).to_numpy()
        lead["p_block"] = 1 - np.power(1 - lead.pval.clip(0, 1), lead.n_hap)
    else:
        lead = d.loc[g.pval.apply(lambda s: s.index[0])].copy()  # arbitrary/first hap
        lead["n_hap"] = g.size().reindex(lead.block).to_numpy()
        lead["p_block"] = lead.pval
    lead["q_block"] = bh(lead.p_block.to_numpy())
    return lead.reset_index(drop=True)


def main():
    ken = pd.read_csv(KEN); reg = pd.read_csv(f"{CM}/hap_gen9.registry.csv")
    print(f"[A] ken rows {len(ken):,} | reg rows {len(reg):,} | "
          f"hap_id identical & aligned: {bool((ken.hap_id.values==reg.hap_id.values).all())}")
    d = ken.merge(reg, on="hap_id", suffixes=("", "_reg"))
    print(f"    merged {len(d):,} | block col == reg block: {bool((d.block==d.block_reg).all())} | "
          f"any pval NaN: {int(d.pval.isna().sum())} | has_sv NaN: {int(d.has_sv.isna().sum())}")
    d = d[np.isfinite(d.pval)].copy()

    # --- direct SV-overlap spot-check: recompute has_sv for 3 blocks from the raw panel ---
    samp = lib.list_samples(lib.OUT)[0]
    meta = pd.read_csv(f"{lib.OUT}/{samp}.tsv", sep="\t", usecols=["chrom","pos","ref_len","alt_len"])
    meta["dlen"] = (meta.alt_len - meta.ref_len).abs(); sv = meta[meta.dlen > 50]
    chk = d.drop_duplicates("block").sample(3, random_state=1)
    print("[A] SV-overlap spot-check (block span vs raw panel SVs):")
    for _, r in chk.iterrows():
        s = sv[(sv.chrom==r.chrom) & (sv.pos < r.block_end) & (sv.pos+sv.ref_len > r.block_start)]
        print(f"    {r.block} [{r.block_start}-{r.block_end}] tagged has_sv={r.has_sv} n_sv={r.n_sv} | "
              f"recomputed overlap={len(s)}  {'OK' if (len(s)>0)==bool(r.has_sv) else 'MISMATCH!!'}")

    # --- (B) representative-choice sensitivity ---
    for rep in ["minp", "firstpos"]:
        blk = block_table(d, rep)
        f_sv = blk[blk.has_sv]; f_sn = blk[~blk.has_sv]
        rate_sv = (blk.q_block[blk.has_sv] < 0.05).mean(); rate_sn = (blk.q_block[~blk.has_sv] < 0.05).mean()
        print(f"[B] rep={rep:8s} blocks={len(blk):,} | %BHsig SV={rate_sv:.1%} SNP={rate_sn:.1%} | "
              f"median -log10p SV={np.median(-np.log10(f_sv.p_block.clip(1e-300,1))):.2f} "
              f"SNP={np.median(-np.log10(f_sn.p_block.clip(1e-300,1))):.2f}")

    # --- (C) common support + SIZE-MATCHED null (min-p representative) ---
    blk = block_table(d, "minp")
    blk["logL"] = np.log10((blk.block_end - blk.block_start).clip(lower=1))
    blk["logV"] = np.log10(blk.unit_nvar.clip(lower=1))
    blk["nlp"] = -np.log10(blk.p_block.clip(1e-300, 1))
    blk["sig"] = (blk.q_block < 0.05).astype(int)
    sv_m = blk.has_sv.to_numpy()
    print(f"\n[C] common support in log10(blocklen):")
    for q in [0, 25, 50, 75, 90, 100]:
        print(f"    pctl {q:3d}:  SV={np.percentile(blk.logL[sv_m], q):.2f}  "
              f"SNP-only={np.percentile(blk.logL[~sv_m], q):.2f}")

    # nearest-neighbour matching: each SV block -> nearest SNP-only by (logL,logV)
    A = blk[sv_m]; B = blk[~sv_m].reset_index(drop=True)
    Bxy = np.column_stack([B.logL, B.logV]); used = np.zeros(len(B), bool)
    mi = []
    order = A.nlp.to_numpy().argsort()          # deterministic order
    for idx in A.index.to_numpy()[order]:
        ax = np.array([blk.logL[idx], blk.logV[idx]])
        dd = ((Bxy - ax) ** 2).sum(1); dd[used] = np.inf
        j = int(dd.argmin()); used[j] = True; mi.append((idx, B.index[j]))
    ai = [m[0] for m in mi]; bj = [m[1] for m in mi]
    sv_nlp = blk.loc[ai, "nlp"].to_numpy(); mm_nlp = B.loc[bj, "nlp"].to_numpy()
    sv_sig = blk.loc[ai, "sig"].to_numpy();  mm_sig = B.loc[bj, "sig"].to_numpy()
    dL = np.abs(blk.loc[ai, "logL"].to_numpy() - B.loc[bj, "logL"].to_numpy())
    print(f"\n[C] size-matched (nearest SNP-only by logLen+logNvar), {len(ai):,} pairs, "
          f"median |ΔlogLen|={np.median(dL):.3f}:")
    print(f"    %BHsig   SV={sv_sig.mean():.1%}  matched-SNP={mm_sig.mean():.1%}  Δ={sv_sig.mean()-mm_sig.mean():+.3f}")
    w = stats.wilcoxon(sv_nlp, mm_nlp, zero_method="zsplit")
    print(f"    median -log10p  SV={np.median(sv_nlp):.3f}  matched-SNP={np.median(mm_nlp):.3f} | "
          f"paired Wilcoxon p={w.pvalue:.2g}")
    mb = stats.binomtest((sv_sig > mm_sig).sum(), (sv_sig != mm_sig).sum())
    print(f"    McNemar sig-discordant SV>SNP {(sv_sig>mm_sig).sum()} vs SNP>SV {(mm_sig>sv_sig).sum()} "
          f"-> p={mb.pvalue:.2g}")

    # stratified permutation null: shuffle has_sv WITHIN logL deciles, stat = mean nlp diff
    strata = pd.qcut(blk.logL, 10, labels=False, duplicates="drop").to_numpy()
    def stat(mask):
        return np.mean([blk.nlp[mask & (strata==s)].mean() - blk.nlp[~mask & (strata==s)].mean()
                        for s in np.unique(strata)
                        if (mask & (strata==s)).any() and (~mask & (strata==s)).any()])
    obs = stat(sv_m)
    null = np.empty(2000)
    for k in range(2000):
        perm = np.zeros(len(blk), bool)
        for s in np.unique(strata):
            si = np.where(strata==s)[0]; nsv = int(sv_m[si].sum())
            perm[rng.choice(si, nsv, replace=False)] = True
        null[k] = stat(perm)
    p = (np.abs(null) >= abs(obs)).mean()
    print(f"\n[C] within-length-decile permutation: obs mean Δ(-log10p, SV−SNP)={obs:+.3f} | "
          f"perm p={p:.3f} (null sd {null.std():.3f})")


if __name__ == "__main__":
    main()

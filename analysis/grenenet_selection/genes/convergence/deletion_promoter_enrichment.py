#!/usr/bin/env python
"""Which variant TYPES in which REGIONS carry the GEA signal, and which way do they move?

Prompted by the burden test (another line of work): promoter deletions look advantageous in
warm gardens, many of the strongest candidates here are promoter deletions (GPX6, CYP28,
PGIP1, AT2G14910), and CDS insertions look purifying-selected. Rather than test only the
combinations already suspected, this scores the full grid -- deletions and insertions x
promoter / CDS / UTR / intron / TE / intergenic -- with one statistic and one multiple-testing
bar (Bonferroni over the 12 cells), so promoter deletions and CDS insertions are two cells of
one table rather than two hand-picked tests.

  Q1  are hits enriched for deletions at all?
  Q2  within each type, which regions are over-represented among hits? (12 cells)
  Q3  is a region's enrichment specific to one type? (ratio of the two ORs)
  Q4  direction: per cell, the share of LD blocks whose SHORTER allele rises with bio1

Note on CDS insertions and purifying selection: purifying selection acts mostly BEFORE this
test -- it keeps CDS insertions rare, and the comparison universe is only variants common
enough to test (MAF >= 0.05), which already carries that depletion. Q2 therefore asks which
types carry CLIMATE signal among testable variants, not whether a type is constrained.

Hits: every GEA-significant non-SNP record (Bonferroni on any axis, raw p -- memory
`raw-lfmm-over-wza-decision`), MNPs excluded because an equal-length substitution is neither
a deletion nor an insertion. Background: `atac_enrich_background.csv.gz`, 80,000 records drawn
at random from the same testable universe (MAF >= 0.05, small indels + SVs), annotated with
the SAME region classifier (dissection/screen_sig_blocks.classify); hits are removed from it,
so the comparison is hits against non-hits.

Three things could fake an enrichment, and each is controlled:

  size / class  deletion:insertion ratios differ enormously by size and SV insertions are
                harder to call (memory `sv-panel-support-asymmetry`), so odds ratios are
                Mantel-Haenszel over vclass x size-bin x MAF-bin strata;
  MAF           LFMM has most power at intermediate frequency, so hits skew common -- the
                MAF bin is in the strata;
  clustering    hits in one LD block are not independent (1,375 blocks for 2,404 records;
                the atac_enrich lesson), so confidence intervals come from resampling whole
                BLOCKS of hits; the 80k background is effectively independent.

Direction (the burden-test claim) is asked separately: for each hit, the sign of the
garden-level correlation between the SHORTER allele's frequency and bio1. For a deletion the
shorter allele is ALT. Caveat that applies to every "deletion" here: it is a deletion
RELATIVE TO Col-0, not relative to the ancestral state -- there is no outgroup polarisation.

Writes results/deletion_promoter_{enrichment,interaction,direction}.csv.
env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", ".."))); sys.path.insert(0, HERE)
OUT = f"{HERE}/results"
K4 = ["chrom", "pos", "ref_len", "alt_len"]
SIZE_BINS = [0, 1, 9, 49, 499, np.inf]
SIZE_LAB = ["1", "2-9", "10-49", "50-499", "500+"]
MAF_BINS = [0.05, 0.1, 0.2, 0.3, 0.5]
# the classifier's tiers, grouped; 4_exon_noncoding (6 hits) joins UTR, gene deserts join
# intergenic
REGION6 = {"1_CDS": "CDS", "2_UTR": "UTR", "4_exon_noncoding": "UTR", "3_promoter": "promoter",
           "5_intron": "intron", "6_TE": "TE", "7_proximal_intergenic": "intergenic",
           "8_gene_desert": "intergenic"}
REGIONS = ["promoter", "CDS", "UTR", "intron", "TE", "intergenic"]


def annotate(D, maf_col):
    D = D.copy()
    D["deletion"] = D.ref_len > D.alt_len
    D["promoter"] = D.tier.astype(str).eq("3_promoter")
    D["region6"] = D.tier.astype(str).map(REGION6).fillna("intergenic")
    D["abs_size"] = (D.ref_len - D.alt_len).abs()
    D["size_bin"] = pd.cut(D.abs_size, SIZE_BINS, labels=SIZE_LAB)
    D["maf_bin"] = pd.cut(D[maf_col], MAF_BINS, include_lowest=True).astype(str)
    D["stratum"] = D.vclass.astype(str) + "|" + D.size_bin.astype(str) + "|" + D.maf_bin
    return D


def _counts(case, ctrl, x, strata):
    """Per-block case counts and per-stratum control counts, so a bootstrap is a matmul."""
    S = sorted(set(case[strata]) | set(ctrl[strata])); si = {v: i for i, v in enumerate(S)}
    c = ctrl.groupby([ctrl[strata], ctrl[x]]).size()
    cw = np.array([c.get((v, True), 0) for v in S], float)
    cwo = np.array([c.get((v, False), 0) for v in S], float)
    blocks = case.block.astype(str).values; ub, bi = np.unique(blocks, return_inverse=True)
    aw = np.zeros((len(ub), len(S))); awo = np.zeros((len(ub), len(S)))
    col = case[strata].map(si).values; xv = case[x].values.astype(bool)
    np.add.at(aw, (bi[xv], col[xv]), 1); np.add.at(awo, (bi[~xv], col[~xv]), 1)
    return aw, awo, cw, cwo


def _mh(a, b, c, d):
    n = a + b + c + d; ok = n > 0
    num = (a[ok] * d[ok] / n[ok]).sum(); den = (b[ok] * c[ok] / n[ok]).sum()
    return num / den if den else np.nan


def mh_or(case, ctrl, x, strata):
    """Mantel-Haenszel OR of feature x (bool) in cases vs controls, over strata."""
    aw, awo, cw, cwo = _counts(case, ctrl, x, strata)
    return _mh(aw.sum(0), awo.sum(0), cw, cwo)


def block_boot(case, ctrl, x, strata, reps=2000, seed=1):
    """Resample whole LD blocks of hits; the background is effectively independent."""
    aw, awo, cw, cwo = _counts(case, ctrl, x, strata)
    rng = np.random.default_rng(seed); nb = aw.shape[0]
    out = []
    for _ in range(reps):
        w = np.bincount(rng.integers(0, nb, nb), minlength=nb).astype(float)
        out.append(_mh(w @ aw, w @ awo, cw, cwo))
    out = np.array(out)
    p = 2 * min(np.nanmean(out <= 1), np.nanmean(out >= 1))
    return (*np.nanpercentile(out, [2.5, 97.5]), min(p, 1.0))


def row(label, case, ctrl, x, strata):
    lo, hi, bp = block_boot(case, ctrl, x, strata)
    return dict(test=label, n_hits=len(case), hits_with=int(case[x].sum()),
                hit_frac=case[x].mean(), n_bg=len(ctrl), bg_with=int(ctrl[x].sum()),
                bg_frac=ctrl[x].mean(), raw_or=(case[x].mean() / (1 - case[x].mean()))
                / (ctrl[x].mean() / (1 - ctrl[x].mean())),
                mh_or=mh_or(case, ctrl, x, strata), ci_lo=lo, ci_hi=hi, boot_p=bp)


def direction(H):
    """Sign of the garden-level r between the SHORTER allele's frequency and bio1."""
    import lib, axis_clusters as ac, candidate_evidence as CE
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    K = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                      "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    K["row"] = np.arange(len(K))
    E = pd.read_csv(f"{OUT}/evidence_matrix.csv")[K4 + ["store_row"]]
    H = H.merge(E, on=K4, how="left")
    H = H[H.store_row.notna()].copy(); H["store_row"] = H.store_row.astype(int)
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())
    G = CE.site_means(M, p0, H.store_row.values, site); G.columns = H.index
    r = CE.corr_axes(G, A[["bio1"]])["bio1"]          # r for the ALT allele
    # the SHORTER allele: ALT for a deletion, REF for an insertion (sign flips)
    H["r_short_bio1"] = np.where(H.deletion, r.reindex(H.index).values,
                                 -r.reindex(H.index).values)
    H["short_up_warm"] = H.r_short_bio1 > 0
    return H


def main():
    V = pd.read_csv(f"{OUT}/variants_classified.csv", low_memory=False)
    H = V[V.in_gea.astype(bool) & (V.vclass != "mnp")].copy()
    H = annotate(H, "gea_MAF")
    B = pd.read_csv(f"{OUT}/atac_enrich_background.csv.gz")
    B = B.merge(H[K4], on=K4, how="left", indicator=True)
    B = annotate(B[B._merge == "left_only"].drop(columns="_merge"), "maf")
    B["block"] = np.arange(len(B))
    print(f"hits {len(H)} records in {H.block.nunique()} blocks | background {len(B)} non-hits")

    rows = []
    # Q1: are hits enriched for deletions at all?
    rows.append(row("Q1 deletion | all", H, B, "deletion", "stratum"))
    # Q2 grid: within each variant type, which REGIONS are over-represented among hits?
    # 2 types x 6 regions = 12 tests -> Bonferroni alpha 0.05/12 = 0.0042
    for kind, hk, bk in [("deletion", H[H.deletion], B[B.deletion]),
                         ("insertion", H[~H.deletion], B[~B.deletion])]:
        for reg in REGIONS:
            hk = hk.assign(_x=hk.region6.eq(reg)); bk = bk.assign(_x=bk.region6.eq(reg))
            rows.append(row(f"Q2 {reg} | {kind}s", hk, bk, "_x", "stratum"))
    R = pd.DataFrame(rows)
    R["bonferroni_12"] = np.where(R.test.str.startswith("Q2"), (R.boot_p * 12).clip(upper=1), np.nan)

    # Q3: is each region's enrichment specific to one type? Ratio of the deletion OR to the
    # insertion OR, with the SAME block resample applied to both. (A logit with stratum
    # dummies is singular here: many strata hold no hits, i.e. perfect separation.)
    Hd, Hi, Bd, Bi = H[H.deletion], H[~H.deletion], B[B.deletion], B[~B.deletion]
    ub = np.unique(H.block.astype(str)); pos = {bb: i for i, bb in enumerate(ub)}
    bd = np.array([pos[bb] for bb in np.unique(Hd.block.astype(str))])
    bi_ = np.array([pos[bb] for bb in np.unique(Hi.block.astype(str))])
    irows = []
    for reg in REGIONS:
        cd = _counts(Hd.assign(_x=Hd.region6.eq(reg)), Bd.assign(_x=Bd.region6.eq(reg)), "_x", "stratum")
        ci_ = _counts(Hi.assign(_x=Hi.region6.eq(reg)), Bi.assign(_x=Bi.region6.eq(reg)), "_x", "stratum")
        o_d = _mh(cd[0].sum(0), cd[1].sum(0), cd[2], cd[3])
        o_i = _mh(ci_[0].sum(0), ci_[1].sum(0), ci_[2], ci_[3])
        rng = np.random.default_rng(2); boot = []
        for _ in range(2000):
            w = np.bincount(rng.integers(0, len(ub), len(ub)), minlength=len(ub)).astype(float)
            boot.append(_mh(w[bd] @ cd[0], w[bd] @ cd[1], cd[2], cd[3])
                        / _mh(w[bi_] @ ci_[0], w[bi_] @ ci_[1], ci_[2], ci_[3]))
        boot = np.array(boot)
        irows.append(dict(region=reg, or_deletion=o_d, or_insertion=o_i, ratio=o_d / o_i,
                          ci_lo=np.nanpercentile(boot, 2.5), ci_hi=np.nanpercentile(boot, 97.5),
                          boot_p=min(1.0, 2 * min(np.nanmean(boot <= 1), np.nanmean(boot >= 1)))))
    I = pd.DataFrame(irows)

    # Q4: direction among hits, per cell -- one vote per LD block
    from scipy.stats import binomtest
    D = direction(H)
    D["cell"] = D.region6 + " " + np.where(D.deletion, "deletion", "insertion")
    Db = D.groupby(["cell", "block"]).short_up_warm.mean().reset_index()
    drows = []
    for c, g in Db.groupby("cell"):
        up, tie = int((g.short_up_warm > 0.5).sum()), int((g.short_up_warm == 0.5).sum())
        n = len(g) - tie
        drows.append(dict(cell=c, n_records=int((D.cell == c).sum()), n_blocks=len(g),
                          blocks_shorter_up_warm=up / n if n else np.nan,
                          sign_p=binomtest(up, n).pvalue if n else np.nan))
    Dd = pd.DataFrame(drows).sort_values("n_blocks", ascending=False)

    R.to_csv(f"{OUT}/deletion_promoter_enrichment.csv", index=False)
    I.to_csv(f"{OUT}/deletion_promoter_interaction.csv", index=False)
    Dd.to_csv(f"{OUT}/deletion_promoter_direction.csv", index=False)
    pd.set_option("display.width", 240)
    print(R.round(3).to_string(index=False))
    print("\nis the region enrichment specific to one type? (OR deletions / OR insertions)")
    print(I.round(3).to_string(index=False))
    print("\ndirection among hits -- share of LD blocks where the SHORTER allele rises in warm "
          "(high bio1) gardens; for insertions that means the insertion allele FALLS")
    print(Dd.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

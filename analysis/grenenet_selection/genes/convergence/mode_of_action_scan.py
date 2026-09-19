#!/usr/bin/env python
"""Beyond promoter deletions: score every coding / splice / UTR / CDS-repeat variant in
the non-SNP pool on the same evidence used for GPX6.

The 336-gene review and the mechanistic shortlist were built around promoter structure,
because that is what GPX6 is. Other modes of action are just as real and the pool already
carries transcript-level effects for every record (results/repeat_variants.csv, from
build_sv_hit_dossier.site_effect): truncating, frameshift, in-frame, whole-CDS deletion,
CDS span, splice region, exon-intron span, stop-lost, UTR, and repeat-copy changes.

For each pooled key with a coding/splice/UTR effect:
  1. resolve the af_store row WITHOUT re-joining on chrom:pos
       GEA keys  -> the row whose gen9 pool MAF equals the GEA record's MAF
       GWAS keys -> the row at `key_rank` among rows sharing the key
     (memory `panel-multiallelic-pos-key-trap`; same rules as audit_gea_colliders.py and
     resolve_gwas_alleles.py)
  2. garden-level gradient, best of 22 axes, against a selection-aware background
  3. genome-wide share of its garden profile and the locus-specific (residual) gradient
  4. local rank among non-SNP records and SNP rivals within +-25 kb, on its best axis
  5. founder carriers / called with split records merged

Writes results/mode_of_action_scan.csv.
env: kmate. Compute node, ~20-40 min.
"""
import os, sys
import numpy as np, pandas as pd, pysam
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", ".."))); sys.path.insert(0, HERE)
import lib, axis_clusters as ac                                  # noqa: E402
import candidate_evidence as CE                                  # noqa: E402

OUT = f"{HERE}/results"
PROJ = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
K4 = ["chrom", "pos", "ref_len", "alt_len"]
SEV = {"CDS:whole-CDS-deleted": 9, "CDS:span": 8, "CDS:truncating": 7, "CDS:frameshift": 6,
       "CDS:stop-lost": 6, "CDS:in-frame": 5, "CDS:altered": 4, "splice-region": 4,
       "exon-intron span": 4, "CDS:synonymous": 2, "5'UTR": 3, "3'UTR": 3,
       "exon(non-coding)": 3, "intron": 2, "transcript(other)": 1}
MODE = {"CDS:whole-CDS-deleted": "gene loss", "CDS:span": "coding SV",
        "CDS:truncating": "truncating", "CDS:frameshift": "truncating",
        "CDS:stop-lost": "stop-lost", "CDS:in-frame": "in-frame", "CDS:altered": "in-frame",
        "splice-region": "splice", "exon-intron span": "splice",
        "5'UTR": "5'UTR", "3'UTR": "3'UTR"}


def main():
    R = pd.read_csv(f"{OUT}/repeat_variants.csv")
    R["sev"] = R.site_feature.map(SEV).fillna(0)
    B = R.sort_values("sev", ascending=False).drop_duplicates(K4)
    B["mode"] = B.site_feature.map(MODE)
    cds_rep = B.repeat_copy_change.astype(bool) & B.site_feature.fillna("").str.startswith("CDS")
    B.loc[cds_rep, "mode"] = "CDS repeat copy"
    B = B[B["mode"].notna()].copy()
    print(f"{len(B)} pooled keys with a coding / splice / UTR effect:")
    print(B["mode"].value_counts().to_string(), flush=True)

    # ---- 1. allele identity --------------------------------------------------------------
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    K = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                      "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    K["row"] = np.arange(len(K))
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    GP = pd.read_csv(f"{OUT}/gea_pool.csv")[K4 + ["MAF"]]
    GW = pd.read_csv(f"{OUT}/gwas_hits_allele_resolved.csv")[K4 + ["key_rank"]].drop_duplicates()
    rows, how = [], []
    for b in B.itertuples():
        hit = K[(K.chrom == b.chrom) & (K.pos == b.pos) & (K.ref_len == b.ref_len)
                & (K.alt_len == b.alt_len)].row.values
        if not len(hit):
            rows.append(-1); how.append("no store row"); continue
        if len(hit) == 1:
            rows.append(int(hit[0])); how.append("unique"); continue
        g = GP[(GP.chrom == b.chrom) & (GP.pos == b.pos) & (GP.ref_len == b.ref_len)
               & (GP.alt_len == b.alt_len)]
        if len(g):
            mu = np.array([np.nanmean(np.asarray(M[:, h], float)) for h in hit])
            d = np.abs(np.minimum(mu, 1 - mu) - float(g.MAF.iloc[0]))
            if d.min() < 1e-6:
                rows.append(int(hit[d.argmin()])); how.append("gea MAF"); continue
        w = GW[(GW.chrom == b.chrom) & (GW.pos == b.pos) & (GW.ref_len == b.ref_len)
               & (GW.alt_len == b.alt_len)]
        if len(w) == 1 and int(w.key_rank.iloc[0]) < len(hit):
            rows.append(int(np.sort(hit)[int(w.key_rank.iloc[0])])); how.append("gwas key_rank"); continue
        rows.append(-1); how.append("ambiguous, unresolved")
    B["store_row"], B["row_resolved_by"] = rows, how
    print("\nallele resolution:", B.row_resolved_by.value_counts().to_dict(), flush=True)
    B = B[B.store_row >= 0].reset_index(drop=True)

    # ---- 2-3. gradient, selection-aware, genome-wide share ------------------------------
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())
    rng = np.random.default_rng(7)
    Bg = pd.read_csv(f"{OUT}/repeat_variants_background_v2.csv.gz", usecols=K4)
    brow = np.sort(rng.choice(K.merge(Bg, on=K4)["row"].values, 3000, replace=False))
    Gb = CE.site_means(M, p0, brow, site)
    bg_sel = np.sort(CE.corr_axes(Gb, A).abs().max(axis=1).values)
    z = lambda G: ((G - G.mean()) / G.std()).fillna(0)
    U, S, _ = np.linalg.svd(z(Gb).values, full_matrices=False)
    P = U[:, :3]; H = P @ np.linalg.pinv(P)
    Rb = pd.DataFrame(z(Gb).values - H @ z(Gb).values, index=Gb.index)
    bg_res = np.sort(CE.corr_axes(Rb, A).abs().max(axis=1).values)

    Gc = CE.site_means(M, p0, B.store_row.values, site); Gc.columns = B.index
    RC = CE.corr_axes(Gc, A)
    B["best_axis"] = RC.abs().idxmax(axis=1).values
    B["r_best"] = RC.abs().max(axis=1).values
    B["pct_sel"] = 100 * np.searchsorted(bg_sel, B.r_best.values) / len(bg_sel)
    X = z(Gc).values; fit = H @ X
    B["gw_share"] = (fit ** 2).sum(0) / np.maximum((X ** 2).sum(0), 1e-12)
    Rc = pd.DataFrame(X - fit, index=Gc.index, columns=Gc.columns)
    rr = CE.corr_axes(Rc, A).abs()
    B["r_resid"] = [rr.loc[i, ax] for i, ax in zip(B.index, B.best_axis)]
    B["pct_resid"] = 100 * np.searchsorted(bg_res, B.r_resid.values) / len(bg_res)
    print("gradients done", flush=True)

    # ---- 4. local rank + SNP rivals -------------------------------------------------------
    Ms = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_snp_af.npy", mmap_mode="r")
    p0s = np.load(f"{lib.AF_STORE}/p0_snp.npy")
    isn = np.load(f"{lib.AF_STORE}/index_snp.npz", allow_pickle=True)
    s_ch, s_pos = isn["chrom"].astype(str), isn["pos"]
    lr, sr = [], []
    for c in B.itertuples():
        a = A[[c.best_axis]]
        sel = np.where((K.chrom.values == c.chrom) & (np.abs(K.pos.values - c.pos) <= 25_000))[0]
        sel = np.union1d(sel[np.minimum(p0[sel], 1 - p0[sel]) > 0.05], [c.store_row])
        rn = CE.corr_axes(CE.site_means(M, p0, sel, site), a)[c.best_axis].abs()
        me = rn.iloc[int(np.where(sel == c.store_row)[0][0])]
        lr.append(int((rn > me + 1e-9).sum()) + 1)
        ss = np.where((s_ch == c.chrom) & (np.abs(s_pos - c.pos) <= 25_000))[0]
        ss = ss[np.minimum(p0s[ss], 1 - p0s[ss]) > 0.05]
        rs = CE.corr_axes(CE.site_means(Ms, p0s, ss, site), a)[c.best_axis].abs() \
            if len(ss) else pd.Series([], dtype=float)
        sr.append(int((rs >= me - 1e-9).sum()))
    B["local_rank"], B["snp_rivals"] = lr, sr
    print("local ranks done", flush=True)

    # ---- 5. founder support ---------------------------------------------------------------
    car, cal, vfs = [], [], {}
    for c in B.itertuples():
        f = f"{PROJ}/panel/arch3/{c.chrom.lower()}/merged_231_{c.chrom.lower()}_final.vcf.gz"
        vf = vfs.setdefault(f, pysam.VariantFile(f))
        a, b = CE.founder_support(vf, c.chrom, int(c.pos), int(c.ref_len), int(c.alt_len))
        car.append(a); cal.append(b)
    B["carriers"], B["called"] = car, cal
    B["size"] = B.alt_len - B.ref_len
    keep = ["mode", "site_feature", "chrom", "pos", "ref_len", "alt_len", "size",
            "site_effect_gene", "protein_effect", "protein_pos", "wt_aa", "mut_aa",
            "splice_dist", "repeat_copy_change", "unit", "n_units_changed",
            "in_gea", "in_gwas", "gea_best_axis", "gea_nlp", "gwas_n_gardens",
            "store_row", "row_resolved_by", "best_axis", "r_best", "pct_sel",
            "gw_share", "r_resid", "pct_resid", "local_rank", "snp_rivals",
            "carriers", "called"]
    B[[c for c in keep if c in B.columns]].to_csv(f"{OUT}/mode_of_action_scan.csv", index=False)
    print(f"wrote {OUT}/mode_of_action_scan.csv ({len(B)} variants)")


if __name__ == "__main__":
    main()

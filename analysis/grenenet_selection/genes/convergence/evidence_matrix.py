#!/usr/bin/env python
"""One row per variant in the non-SNP pool, one column per line of evidence.

The candidate lists so far were each built around one shape (GPX6-like promoter deletions,
then coding / splice / UTR effects) and each covered part of the pool. This scores ALL
3,020 pooled keys on every line we now have, so candidates are ranked by how many
INDEPENDENT lines they satisfy -- never by p (README "Why"; raw p, no inflation correction,
memory `raw-lfmm-over-wza-decision`).

Stages (results/evidence_matrix_*.csv):
  a   allele resolution (gen9 MAF for GEA keys, key_rank for GWAS keys -- never chrom:pos),
      site effect / mode of action, target gene and TSS distance, garden-level gradient
      (best of 22 axes, selection-aware; genome-wide share; locus-specific residual),
      ATAC, founder carriers / called (split records merged), shared-AF flag
  b   local rank among non-SNP records and SNP rivals within +-25 kb (slow; --worker/--nworkers)
  c   carrier-vs-non-carrier expression (needs data/eqtl/expr_*.csv from
      export_expression_long.R for the stage-a target genes), with lineage control
  d   merge + lines-of-evidence scoring -> results/evidence_matrix.csv

env: kmate. Compute node.
"""
import os, sys, argparse
import numpy as np, pandas as pd, pysam
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", ".."))); sys.path.insert(0, HERE)
import lib, axis_clusters as ac                                  # noqa: E402
import candidate_evidence as CE                                  # noqa: E402
import atac_overlap as AO                                        # noqa: E402

OUT = f"{HERE}/results"
PROJ = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
K4 = ["chrom", "pos", "ref_len", "alt_len"]
SEV = {"CDS:whole-CDS-deleted": 9, "CDS:span": 8, "CDS:truncating": 7, "CDS:frameshift": 6,
       "CDS:stop-lost": 6, "CDS:in-frame": 5, "CDS:altered": 4, "splice-region": 4,
       "exon-intron span": 4, "CDS:synonymous": 2, "5'UTR": 3, "3'UTR": 3,
       "exon(non-coding)": 3, "intron": 2, "transcript(other)": 1}


def mode_of(r):
    e, f = str(r.protein_effect), str(r.site_feature)
    if f == "CDS:whole-CDS-deleted":
        return "gene loss"
    if f.startswith("CDS"):
        import re
        m = re.search(r"ends at (\d+) of (\d+)", e)
        if ("frameshift" in e and (not m or int(m.group(1)) / int(m.group(2)) < 0.9)) or \
           (m and int(m.group(1)) / int(m.group(2)) < 0.9) or f == "CDS:span":
            return "truncation"
        return "in-frame / minor coding"
    if f in ("splice-region", "exon-intron span"):
        return "splice"
    if f in ("5'UTR", "3'UTR"):
        return f
    t = str(r.tier)
    if t.startswith("3_promoter"):
        return "promoter"
    if f == "intron" or t.startswith("5_intron"):
        return "intron"
    if t.startswith("6_TE"):
        return "TE"
    return "intergenic"


def stage_a():
    V = pd.read_csv(f"{OUT}/variants_classified.csv")
    V["chrom"] = V.chrom.str.replace("^chr", "Chr", regex=True)
    R = pd.read_csv(f"{OUT}/repeat_variants.csv")
    R["sev"] = R.site_feature.map(SEV).fillna(0)
    R = R.sort_values("sev", ascending=False).drop_duplicates(K4)
    V = V.merge(R[K4 + ["site_feature", "protein_effect", "site_effect_gene",
                        "repeat_copy_change"]], on=K4, how="left")
    V["mode"] = V.apply(mode_of, axis=1)
    V.loc[V.vclass == "mnp", "mode"] = "MNP (" + V.loc[V.vclass == "mnp", "mode"] + ")"
    V["target_gene"] = V.site_effect_gene.where(V.site_effect_gene.notna(), V.gene)
    V["target_gene"] = V.target_gene.where(V.target_gene.notna(), V.nearest_gene)

    # ---- allele identity ------------------------------------------------------------------
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    K = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                      "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    K["row"] = np.arange(len(K))
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    GP = pd.read_csv(f"{OUT}/gea_pool.csv")[K4 + ["MAF"]]
    GW = pd.read_csv(f"{OUT}/gwas_hits_allele_resolved.csv")[K4 + ["key_rank"]].drop_duplicates()
    Kg = K.groupby(K4).row.apply(lambda x: np.sort(x.values)).to_dict()
    rows, how = [], []
    for v in V.itertuples():
        hit = Kg.get((v.chrom, v.pos, v.ref_len, v.alt_len))
        if hit is None:
            rows.append(-1); how.append("no store row"); continue
        if len(hit) == 1:
            rows.append(int(hit[0])); how.append("unique"); continue
        g = GP[(GP.chrom == v.chrom) & (GP.pos == v.pos) & (GP.ref_len == v.ref_len)
               & (GP.alt_len == v.alt_len)]
        if len(g):
            mu = np.array([np.nanmean(np.asarray(M[:, h], float)) for h in hit])
            d = np.abs(np.minimum(mu, 1 - mu) - float(g.MAF.iloc[0]))
            if d.min() < 1e-6:
                rows.append(int(hit[d.argmin()])); how.append("gea MAF"); continue
        w = GW[(GW.chrom == v.chrom) & (GW.pos == v.pos) & (GW.ref_len == v.ref_len)
               & (GW.alt_len == v.alt_len)]
        if len(w) == 1 and int(w.key_rank.iloc[0]) < len(hit):
            rows.append(int(hit[int(w.key_rank.iloc[0])])); how.append("gwas key_rank"); continue
        rows.append(-1); how.append("ambiguous, unresolved")
    V["store_row"], V["row_resolved_by"] = rows, how
    print("allele resolution:", V.row_resolved_by.value_counts().to_dict(), flush=True)
    V = V[V.store_row >= 0].reset_index(drop=True)

    # ---- gradients ------------------------------------------------------------------------
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())
    rng = np.random.default_rng(11)
    Bg = pd.read_csv(f"{OUT}/repeat_variants_background_v2.csv.gz", usecols=K4)
    brow = np.sort(rng.choice(K.merge(Bg, on=K4)["row"].values, 4000, replace=False))
    Gb = CE.site_means(M, p0, brow, site)
    bg_sel = np.sort(CE.corr_axes(Gb, A).abs().max(axis=1).values)
    z = lambda G: ((G - G.mean()) / G.std()).fillna(0)
    U, S, _ = np.linalg.svd(z(Gb).values, full_matrices=False)
    P = U[:, :3]; H = P @ np.linalg.pinv(P)
    bg_res = np.sort(CE.corr_axes(pd.DataFrame(z(Gb).values - H @ z(Gb).values, index=Gb.index),
                                  A).abs().max(axis=1).values)
    Gc = CE.site_means(M, p0, V.store_row.values, site); Gc.columns = V.index
    RC = CE.corr_axes(Gc, A)
    V["best_axis"] = RC.abs().idxmax(axis=1).values
    V["r_best"] = RC.abs().max(axis=1).values
    V["pct_sel"] = 100 * np.searchsorted(bg_sel, V.r_best.values) / len(bg_sel)
    X = z(Gc).values; fit = H @ X
    V["gw_share"] = (fit ** 2).sum(0) / np.maximum((X ** 2).sum(0), 1e-12)
    rr = CE.corr_axes(pd.DataFrame(X - fit, index=Gc.index, columns=Gc.columns), A).abs()
    V["r_resid"] = [rr.loc[i, ax] for i, ax in zip(V.index, V.best_axis)]
    V["pct_resid"] = 100 * np.searchsorted(bg_res, V.r_resid.values) / len(bg_res)

    # ---- TSS distance (negative = upstream) -----------------------------------------------
    G = lib.load_genes().drop_duplicates("gene").set_index("gene")
    d = []
    for v in V.itertuples():
        if v.target_gene in G.index:
            g = G.loc[v.target_gene]
            d.append((g.end - v.pos) * -1 if g.strand == "-" and False else
                     ((g.end - v.pos) if g.strand == "-" else (v.pos - g.start)) * 1)
        else:
            d.append(np.nan)
    V["dist_tss"] = d          # >0 = downstream of the TSS in gene orientation, <0 upstream
    V.loc[V.dist_tss.notna(), "dist_tss"] = V.loc[V.dist_tss.notna(), "dist_tss"]

    # ---- ATAC -----------------------------------------------------------------------------
    peaks = AO.load_peaks()
    ov, rem, tis = [], [], []
    for v in V.itertuples():
        ps, pe, num, _ = peaks[v.chrom]
        lo, hi = v.pos, v.pos + max(int(v.ref_len), 1) - 1
        h = np.where((pe >= lo) & (ps <= hi))[0]
        ov.append(int(sum(min(pe[j], hi) - max(ps[j], lo) + 1 for j in h)))
        rem.append(int(sum((ps[j] >= lo) and (pe[j] <= hi) for j in h)))
        tis.append(int(num[h].max()) if len(h) else 0)
    V["atac_bp"], V["atac_removed"], V["atac_tissues"] = ov, rem, tis

    # ---- founder support ------------------------------------------------------------------
    car, cal, vfs = [], [], {}
    for v in V.itertuples():
        f = f"{PROJ}/panel/arch3/{v.chrom.lower()}/merged_231_{v.chrom.lower()}_final.vcf.gz"
        vf = vfs.setdefault(f, pysam.VariantFile(f))
        a, b = CE.founder_support(vf, v.chrom, int(v.pos), int(v.ref_len), int(v.alt_len))
        car.append(a); cal.append(b)
    V["carriers"], V["called"] = car, cal

    # ---- shared pool-AF flag --------------------------------------------------------------
    D = pd.read_csv(f"{OUT}/af_duplicate_groups.csv")
    V["af_shared_vector"] = V.store_row.isin(set(D.row))
    V.to_csv(f"{OUT}/evidence_matrix_A.csv", index=False)
    print(f"stage a: {len(V)} variants; {V.target_gene.nunique()} target genes", flush=True)
    open(f"{OUT}/_evidence_genes.txt", "w").write(" ".join(sorted(V.target_gene.dropna().unique())))


def stage_b(worker, nworkers):
    V = pd.read_csv(f"{OUT}/evidence_matrix_A.csv").iloc[worker::nworkers]
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    kch, kpos = idx["chrom"].astype(str), idx["pos"]
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    Ms = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_snp_af.npy", mmap_mode="r")
    p0s = np.load(f"{lib.AF_STORE}/p0_snp.npy")
    isn = np.load(f"{lib.AF_STORE}/index_snp.npz", allow_pickle=True)
    s_ch, s_pos = isn["chrom"].astype(str), isn["pos"]
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())
    out = []
    for n, v in enumerate(V.itertuples()):
        a = A[[v.best_axis]]
        sel = np.where((kch == v.chrom) & (np.abs(kpos - v.pos) <= 25_000))[0]
        sel = np.union1d(sel[np.minimum(p0[sel], 1 - p0[sel]) > 0.05], [v.store_row])
        rn = CE.corr_axes(CE.site_means(M, p0, sel, site), a)[v.best_axis].abs()
        me = rn.iloc[int(np.where(sel == v.store_row)[0][0])]
        ss = np.where((s_ch == v.chrom) & (np.abs(s_pos - v.pos) <= 25_000))[0]
        ss = ss[np.minimum(p0s[ss], 1 - p0s[ss]) > 0.05]
        rs = CE.corr_axes(CE.site_means(Ms, p0s, ss, site), a)[v.best_axis].abs() \
            if len(ss) else pd.Series([], dtype=float)
        out.append(dict(store_row=v.store_row, local_rank=int((rn > me + 1e-9).sum()) + 1,
                        n_local=len(sel), snp_rivals=int((rs >= me - 1e-9).sum()),
                        n_snp_window=len(ss)))
        if n % 100 == 0:
            print(f"  worker {worker}: {n}/{len(V)}", flush=True)
    pd.DataFrame(out).to_csv(f"{OUT}/evidence_matrix_B_{worker}.csv", index=False)
    print(f"stage b worker {worker} done", flush=True)


def stage_c():
    import carrier_expression_test as CT
    V = pd.read_csv(f"{OUT}/evidence_matrix_A.csv")
    Gx = pd.read_csv(f"{CT.EQ}/expr_genes.csv")
    Cx = pd.read_csv(f"{CT.EQ}/expr_control.csv")
    graw = Gx.pivot_table(index="acc", columns="id", values="raw", aggfunc="mean")
    ctl = Cx.pivot_table(index="acc", columns="id", values="value", aggfunc="mean")
    med = graw.median()
    vfs, out = {}, []
    for v in V.itertuples():
        o = dict(store_row=v.store_row)
        g = v.target_gene
        if not isinstance(g, str) or g not in graw.columns:
            o["expr_status"] = "not in 1001T"; out.append(o); continue
        o["gene_median_raw"] = float(med[g])
        f = f"{PROJ}/panel/arch3/{v.chrom.lower()}/merged_231_{v.chrom.lower()}_final.vcf.gz"
        vf = vfs.setdefault(f, pysam.VariantFile(f))
        car, ref = CT.carriers(vf, v.chrom, int(v.pos), int(v.ref_len), int(v.alt_len))
        r = graw[g]
        cc = [a for a in r.index if a in car and np.isfinite(r[a])]
        rr = [a for a in r.index if a in ref and np.isfinite(r[a])]
        o.update(n_car_expr=len(cc), n_ref_expr=len(rr),
                 med_raw_car=r.reindex(cc).median(), med_raw_ref=r.reindex(rr).median())
        o["expr_fold"] = (o["med_raw_car"] + 1) / (o["med_raw_ref"] + 1) if cc and rr else np.nan
        p = CT.mw(r.reindex(cc).values, r.reindex(rr).values)
        o["p_expr"] = p
        if np.isfinite(p):
            pc = np.array([CT.mw(ctl[c].reindex(cc).values, ctl[c].reindex(rr).values)
                           for c in ctl.columns])
            pc = pc[np.isfinite(pc)]
            o["p_expr_emp"] = float((np.sum(pc <= p) + 1) / (len(pc) + 1))
            o["expr_status"] = "tested"
        else:
            o["expr_status"] = f"too few expressed carriers ({len(cc)})"
        out.append(o)
    pd.DataFrame(out).to_csv(f"{OUT}/evidence_matrix_C.csv", index=False)
    print("stage c done:", pd.DataFrame(out).expr_status.value_counts().to_dict())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--nworkers", type=int, default=1)
    a = ap.parse_args()
    {"a": stage_a, "c": stage_c}.get(a.stage, lambda: stage_b(a.worker, a.nworkers))()

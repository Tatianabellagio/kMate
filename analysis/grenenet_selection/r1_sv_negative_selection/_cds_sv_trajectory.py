#!/usr/bin/env python
"""Per-garden linear fit of frequency vs generation for CDS SVs, against matched SNPs.

The burden model says founders carrying more coding SVs lose fitness, and more so in warm gardens.
This asks the allele-level version directly: per garden, does the frequency of CDS SVs fall over the
experimental generations, relative to SNPs matched on starting frequency AND carrier-origin climate?

  unit      arch3 allele RECORD of a CDS SV event (~4 records per event; the burden model's unit is the
            event, so events here are weighted by how many assembly paths arch3 split them into)
  response  dp = p(t) - p0, p0 = seed-mix mean; dp(0) = 0 is anchored into every fit
  matching  p0 decile x carrier-origin-climate quintile; carrier_bio1 = mean origin bio1 of the founders
            carrying the variant, from var_pa x sv_origin_partition_founders.csv. This is the confound
            that frequency matching alone cannot remove: cold-origin ecotypes decline as gardens warm,
            so any variant they carry falls with them (sv_founder_confound.ipynb, rho = -0.542).
  statistic per variant, the slope of dp on generation within a garden; then the garden's mean slope,
            Wilcoxon vs 0, and Mann-Whitney against the matched SNPs' slopes. A single pooled regression
            per garden would treat ~13.6k correlated variants as independent -- not done.

Stage 1 (matched sets) is cached to results/sv_adaptive/cds_traj_sets.npz; delete it to rebuild.
Writes results/sv_adaptive/cds_traj_pergarden.csv.  env: kmate. Compute node, ~20 min cold.
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import lib                                                                    # noqa: E402

GM = f"{lib.GEA}/common/results/gen_matrices"
OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
CACHE = f"{OUT}/cds_traj_sets.npz"
K = ["chrom", "pos", "ref_len", "alt_len"]
GENS = [1, 2, 3]
N_SNP_POOL, PER_SV, NBP0, NBCB = 400_000, 10, 10, 5


def carrier_bio1(idx, fb):
    """mean origin bio1 of the founders carrying each record in `idx` (cols K)."""
    out = np.full(len(idx), np.nan)
    for c in ("chr1", "chr2", "chr3", "chr4", "chr5"):
        C = "Chr" + c[-1]
        sel = idx.chrom.to_numpy() == C
        if not sel.any():
            continue
        base = f"{lib.PROJ}/panel/arch3/{c}/var_pa_231_arch3_{c}"
        m = np.load(f"{base}.meta.npz", allow_pickle=True)
        key = pd.DataFrame({"chrom": C, "pos": m["pos"].astype(np.int64),
                            "ref_len": m["ref_len"].astype(np.int64),
                            "alt_len": m["alt_len"].astype(np.int64),
                            "j": np.arange(len(m["pos"]))})
        key = key.drop_duplicates(subset=K, keep="first")     # panel carries duplicate 4-keys
        j = idx[sel].merge(key, on=K, how="left", validate="many_to_one").j.to_numpy()
        ok = np.isfinite(j)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
        sub = (vp[:, j[ok].astype(int)] > 0).toarray()
        n = sub.sum(0)
        tmp = np.full(sel.sum(), np.nan)
        tmp[ok] = np.where(n > 0, (sub * fb[:, None]).sum(0) / np.maximum(n, 1), np.nan)
        out[sel] = tmp
        print(f"  carrier_bio1 {C}: {int(ok.sum()):,}/{int(sel.sum()):,}", flush=True)
    return out


def build_sets():
    fb = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str}).bio1.to_numpy()
    z = np.load(f"{lib.GEA}/common/results/group_means.npz")
    P0 = pd.DataFrame({k: z[k] for k in (*K, "p0")})
    P0["chrom"] = P0.chrom.astype(str)
    P0 = P0.drop_duplicates(subset=K, keep="first")           # 4.76% duplicate keys
    A = pd.read_csv(f"{OUT}/sv_purging_annotated_event2_cr70.csv")
    cds = A[(A.cls.isin(["insertion", "deletion"])) & (A.context == "CDS")][["chrom", "pos"]].drop_duplicates()
    SV = pd.read_csv(f"{GM}/sv_column_index.csv").merge(cds, on=["chrom", "pos"], how="inner")
    SV = SV.merge(P0, on=K, how="left", validate="many_to_one").dropna(subset=["p0"])
    SV["cb"] = carrier_bio1(SV[K], fb)
    SV = SV.dropna(subset=["cb"])
    print(f"CDS SV records: {len(SV):,}", flush=True)
    rng = np.random.default_rng(0)
    SNi = pd.read_csv(f"{GM}/snp_column_index.csv")
    SNi = SNi.iloc[rng.choice(len(SNi), N_SNP_POOL, replace=False)].reset_index(drop=True)
    SN = SNi.merge(P0, on=K, how="left", validate="many_to_one").dropna(subset=["p0"])
    SN["cb"] = carrier_bio1(SN[K], fb)
    SN = SN.dropna(subset=["cb"]).reset_index(drop=True)
    print(f"SNP pool: {len(SN):,}", flush=True)
    pe = np.quantile(np.r_[SV.p0, SN.p0], np.linspace(0, 1, NBP0 + 1))[1:-1]
    ce = np.quantile(np.r_[SV.cb, SN.cb], np.linspace(0, 1, NBCB + 1))[1:-1]
    cell = lambda d: np.digitize(d.p0, pe) * 10 + np.digitize(d.cb, ce)
    SV["cell"], SN["cell"] = cell(SV), cell(SN)
    pick = []
    for cl, k in (SV.cell.value_counts() * PER_SV).items():
        pool = SN.index[SN.cell == cl].to_numpy()
        if pool.size:
            pick.append(rng.choice(pool, min(int(k), pool.size), replace=False))
    SNm = SN.loc[np.concatenate(pick)]
    print(f"matched SNPs: {len(SNm):,}  (SV cells covered {SV.cell.isin(SNm.cell).mean():.1%})", flush=True)
    np.savez(CACHE, sv_col=SV.col.to_numpy(), sv_p0=SV.p0.to_numpy(),
             snp_col=SNm.col.to_numpy(), snp_p0=SNm.p0.to_numpy())
    print(f"[cached] {CACHE}", flush=True)


def site_means(kind, cols):
    """(gen -> (sites, site x variant mean AF)); samples -> pools (flower-weighted) -> site mean."""
    out = {}
    for g in GENS:
        M = np.load(f"{GM}/gen{g}_{kind}_af.npy", mmap_mode="r")
        rm = pd.read_csv(f"{GM}/gen{g}.rowmeta.csv")
        X = np.asarray(M[:, cols], dtype=np.float32)
        w = rm["flowerscollected"].to_numpy(float)
        w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        pk = (rm["site"].astype(int).astype(str) + "_" + rm["plot"].astype(str)).to_numpy()
        num = pd.DataFrame(X * w[:, None]).groupby(pk).sum()
        den = pd.DataFrame(np.isfinite(X) * w[:, None]).groupby(pk).sum()
        pool = num / den.replace(0, np.nan)
        sm = pool.groupby([int(i.split("_")[0]) for i in pool.index]).mean()
        out[g] = (sm.index.to_numpy(), sm.to_numpy(dtype=np.float32))
        print(f"  {kind} gen{g}: {X.shape[0]} samples -> {pool.shape[0]} pools -> {sm.shape[0]} sites", flush=True)
        del X, M
    return out


def main():
    if not os.path.exists(CACHE):
        build_sets()
    S = np.load(CACHE)
    sets = {}
    for kind in ("sv", "snp"):
        c = S[f"{kind}_col"]; o = np.argsort(c)
        sets[kind] = (c[o], S[f"{kind}_p0"][o])
    sm = {k: site_means(k, cols) for k, (cols, _) in sets.items()}
    bio = np.load(f"{lib.GEA}/r3_persite_gwas/results/ecotype_fitness/founder_site_dh.npz", allow_pickle=True)
    BIO = dict(zip(bio["sites"].astype(int), bio["bio1"]))
    rows = []
    for st in sorted(set(sm["sv"][GENS[0]][0])):
        rec = {"site": st, "bio1": BIO.get(st, np.nan)}
        sl_ = {}
        for kind, (_, p0) in sets.items():
            gg, D = [], []
            for g in GENS:
                ss, mat = sm[kind][g]
                j = np.where(ss == st)[0]
                if j.size:
                    gg.append(g); D.append(mat[j[0]] - p0)
            if not D:
                continue
            D = np.vstack(D); gf = np.r_[0.0, np.array(gg, float)]
            sl = np.full(D.shape[1], np.nan)
            for v in range(D.shape[1]):
                y = np.r_[0.0, D[:, v]]; ok = np.isfinite(y)
                if ok.sum() >= 3:
                    sl[v] = np.polyfit(gf[ok], y[ok], 1)[0]
            sl = sl[np.isfinite(sl)]; sl_[kind] = sl
            rec[f"{kind}_n"], rec[f"{kind}_slope"] = sl.size, sl.mean()
            rec[f"{kind}_p"] = stats.wilcoxon(sl).pvalue if sl.size > 10 else np.nan
            rec["n_gen"] = len(gg)
        if {"sv", "snp"} <= sl_.keys():
            rec["diff"] = rec["sv_slope"] - rec["snp_slope"]
            rec["mw_p"] = stats.mannwhitneyu(sl_["sv"], sl_["snp"]).pvalue
        rows.append(rec)
    R = pd.DataFrame(rows).sort_values("bio1")
    R.to_csv(f"{OUT}/cds_traj_pergarden.csv", index=False)
    print(f"\n[wrote] {OUT}/cds_traj_pergarden.csv")
    print(f"gardens {len(R)};  generations sampled: {R.n_gen.value_counts().to_dict()}")
    print(f"SV slope < 0 in {int((R.sv_slope < 0).sum())}/{len(R)};  "
          f"below matched SNPs in {int((R['diff'] < 0).sum())}/{len(R)}")
    for c, l in (("sv_slope", "CDS SV slope"), ("diff", "SV - matched SNP")):
        ok = np.isfinite(R.bio1) & np.isfinite(R[c])
        lr = stats.linregress(R.bio1[ok], R[c][ok])
        print(f"{l:20s} vs bio1: {lr.slope:+.3e}/degC  r={lr.rvalue:+.3f}  p={lr.pvalue:.3g}")
    print("\n" + R[["site", "bio1", "n_gen", "sv_n", "sv_slope", "sv_p",
                    "snp_slope", "diff", "mw_p"]].round(6).to_string(index=False))


if __name__ == "__main__":
    main()

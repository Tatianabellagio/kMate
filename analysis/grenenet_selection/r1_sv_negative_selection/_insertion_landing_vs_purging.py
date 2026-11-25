#!/usr/bin/env python
"""Which insertions are purged in warm gardens, and where do they land?

Joins the per-variant temporal result (the 4,946 MAC>=12 SVs behind
temporal_s_consolidated.ipynb, s per site + climate slope beta) to the insertion-site
context table (insertion_context.csv) and flags the insertions that decline MORE than a
SNP at the same initial frequency, under two definitions:

  A  beta(bio1) below the 10th percentile of p0-matched SNP betas
     (purging intensifies with heat; climate-differential by construction)
  B  per-site s below the p0-matched SNP 10th percentile in >= 3 of the hot gardens
     (bio1 >= HOT_BIO1); the same count at cold gardens (bio1 <= COLD_BIO1) is kept as a
     control, because s is correlated across sites through shared carriers.

The join is ORDINAL, not chrom:pos. The AF store's index_nonsnp is the panel's non-SNP
records in panel order, so each store row maps to its var_pa column (`rec`), which is the
key insertion_context.csv uses. 2.14% of positions carry several records and a position
join silently matches the wrong allele (verified: (chrom,pos,size) inflates 2,546 -> 8,047).

Env: kmate. Reads AF store + panel meta + results/sv_adaptive/{s_climate_slope.npz,
s_dist_by_stratum.npz,_sitemeta.csv,insertion_context.csv}.
Also computes beta for ALL 172k insertions with no MAC / p0 filter, from the per-site s the
no-filter arm stored (s_dist_by_stratum_nofilter.npz = every store SV in order), using the same
across-site OLS slope as production beta.

The same is done for all deletions (context from `_deletion_genomic_context.py`, run it first).

Writes results/sv_adaptive/insertion_landing_vs_purging.{csv,json}, insertion_beta_all.csv,
deletion_beta_all.csv.
"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

R = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
MIN_P0, SV_BP, MIN_MAC = 0.02, 50, 12          # same as _compute_s_dist_by_stratum.py
TAIL, NBIN = 0.10, 20
HOT_BIO1, COLD_BIO1, MIN_HOT = 17.0, 12.0, 3


def pct_in_matched_snp(val_sv, p0_sv, val_snp, p0_snp, nbin=NBIN):
    """Percentile of each SV value within the SNP distribution of its p0 bin."""
    e = np.quantile(np.concatenate([p0_sv, p0_snp]), np.linspace(0, 1, nbin + 1)); e[-1] += 1e-9
    bs = np.digitize(p0_sv, e[1:-1]); bn = np.digitize(p0_snp, e[1:-1])
    out = np.full(val_sv.size, np.nan)
    for b in range(nbin):
        ref = np.sort(val_snp[bn == b]); m = bs == b
        if ref.size and m.any():
            out[m] = np.searchsorted(ref, val_sv[m]) / ref.size
    return out


def main():
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz")
    ch = idx["chrom"].astype("U5"); pos = idx["pos"].astype(np.int64)
    rl = idx["ref_len"].astype(np.int64); al = idx["alt_len"].astype(np.int64)

    print("[1] store row -> panel record ordinal")
    rec = np.full(pos.size, -1, np.int64)
    for c in lib.CHROMS:
        meta = np.load(f"{lib.PROJ}/panel/arch3/{c.lower()}/var_pa_231_arch3_{c.lower()}.meta.npz",
                       allow_pickle=True)
        mp = meta["pos"].astype(np.int64); mr = meta["ref_len"].astype(np.int64); ma = meta["alt_len"].astype(np.int64)
        w = np.where(~((mr == 1) & (ma == 1)))[0]; m = ch == c
        assert w.size == m.sum() and np.array_equal(mp[w], pos[m]) and np.array_equal(mr[w], rl[m]) \
            and np.array_equal(ma[w], al[m]), f"{c}: AF store is not the panel's non-SNP records in order"
        rec[m] = w

    print("[2] reproduce the tested SV set")
    dlen = np.abs(al - rl)
    keep = lib.founder_panel_keep(ch, pos, min_mac=MIN_MAC)
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy").astype(np.float64)
    cols = np.where((dlen >= 1) & keep & (p0 >= MIN_P0) & (p0 <= 1 - MIN_P0))[0]
    sv = cols[dlen[cols] > SV_BP]
    Z = np.load(f"{R}/s_climate_slope.npz")
    assert sv.size == Z["p0_sv"].size and np.allclose(p0[sv], Z["p0_sv"], atol=1e-6)
    assert np.array_equal(rl[sv] > al[sv], Z["sv_isdel"].astype(bool))
    isins = al[sv] > rl[sv]
    print(f"    {sv.size} tested SVs: {isins.sum()} insertions, {(~isins).sum()} deletions")

    J = pd.DataFrame(dict(chrom=ch[sv], rec=rec[sv], pos=pos[sv], size=al[sv] - rl[sv], p0=p0[sv],
                          mean_s=Z["mean_sv"], beta_bio1=Z["beta_sv_bio1"], beta_bio18=Z["beta_sv_bio18"]))[isins]
    J = J.reset_index(drop=True)

    print("[3] ordinal join to insertion_context.csv")
    ctx = pd.read_csv(f"{R}/insertion_context.csv")
    assert not ctx.duplicated(["chrom", "rec"]).any()
    J = J.merge(ctx.drop(columns=["pos", "size"]), on=["chrom", "rec"], how="left", validate="one_to_one")
    assert J.genic.notna().all()
    r = np.corrcoef(J.p0, J.freq)[0, 1]
    print(f"    {len(J)} joined; corr(seed-mix p0, founder freq) = {r:.4f}")

    print("[4] flag A: climate slope vs p0-matched SNPs")
    J["pct_beta1"] = pct_in_matched_snp(J.beta_bio1.values, J.p0.values, Z["beta_snp_bio1"], Z["p0_snp"])
    J["pct_beta18"] = pct_in_matched_snp(J.beta_bio18.values, J.p0.values, Z["beta_snp_bio18"], Z["p0_snp"])
    J["flagA"] = J.pct_beta1 < TAIL
    print(f"    beta(bio1) < SNP {int(TAIL*100)}th pct: {J.flagA.sum()} ({100*J.flagA.mean():.1f}%); "
          f"beta(bio18) > {int((1-TAIL)*100)}th: {(J.pct_beta18 > 1 - TAIL).sum()}")

    print("[5] flag B: per-site s vs p0-matched SNPs")
    D = np.load(f"{R}/s_dist_by_stratum.npz")
    smeta = pd.read_csv(f"{R}/s_dist_by_stratum_sitemeta.csv").set_index("site")
    hot = [int(s) for s in smeta.index[smeta.bio1 >= HOT_BIO1]]
    cold = [int(s) for s in smeta.index[smeta.bio1 <= COLD_BIO1]]
    n_hot = np.zeros(len(J), int); n_cold = np.zeros(len(J), int)
    for s in smeta.index:
        s_sv = D[f"{s}_s_SV"][isins]
        J[f"s_{s}"] = s_sv
        pct = pct_in_matched_snp(s_sv, J.p0.values, D[f"{s}_s_SNP"], D[f"{s}_p0_SNP"])
        J[f"pct_s_{s}"] = pct
        if s in hot: n_hot += pct < TAIL
        if s in cold: n_cold += pct < TAIL
    J["n_hot_low"] = n_hot; J["n_cold_low"] = n_cold
    J["flagB"] = J.n_hot_low >= MIN_HOT
    print(f"    hot sites {hot}; cold sites {cold}")
    print(f"    s < SNP 10th pct in >= {MIN_HOT} hot gardens: {J.flagB.sum()} ({100*J.flagB.mean():.1f}%); "
          f"same rule at cold: {(J.n_cold_low >= MIN_HOT).sum()}; A&B: {(J.flagA & J.flagB).sum()}")

    print("[6] climate slope for ALL insertions, no MAC / p0 filter (per-site s from the no-filter arm)")
    D0 = np.load(f"{R}/s_dist_by_stratum_nofilter.npz")
    sv_all = np.where(dlen > SV_BP)[0]
    assert sv_all.size == D0["1_s_SV"].size and np.allclose(p0[sv_all], D0["1_p0_SV"], atol=1e-6), \
        "no-filter npz is not the store's SV set in order"
    ins_all = al[sv_all] > rl[sv_all]
    PM = f"{lib.GEA}/common/results/pool_matrices"
    cmeta = pd.concat([pd.read_csv(f"{PM}/pool_gen{g}_nonsnp.meta.csv") for g in (1, 2, 3)],
                      ignore_index=True).groupby("site")[["bio1", "bio18"]].mean()
    sites0 = [s for s in sorted({int(k.split("_")[0]) for k in D0.files if k.endswith("_s_SV")})
              if s in cmeta.index]
    S = np.vstack([D0[f"{s}_s_SV"][ins_all].astype(np.float64) for s in sites0])   # [n_sites x n_ins]
    Sc = S - S.mean(0)
    A = dict(chrom=ch[sv_all][ins_all], rec=rec[sv_all][ins_all], p0=p0[sv_all][ins_all],
             mean_s=S.mean(0), n_sites=len(sites0))
    from scipy import stats as st
    n0 = len(sites0)
    for k in ("bio1", "bio18"):
        c = cmeta.loc[sites0, k].to_numpy(); cc = c - c.mean(); sxx = (cc ** 2).sum()
        b = (cc[:, None] * Sc).sum(0) / sxx
        sse = ((Sc - b[None, :] * cc[:, None]) ** 2).sum(0)
        se = np.sqrt(sse / (n0 - 2) / sxx)
        A[f"beta_{k}"] = b; A[f"beta_{k}_se"] = se
        A[f"beta_{k}_p"] = 2 * st.t.sf(np.abs(b / se), df=n0 - 2)
    A = pd.DataFrame(A).merge(ctx, on=["chrom", "rec"], how="left", validate="one_to_one")
    assert A.genic.notna().all() and len(A) == len(ctx)
    neg = A.beta_bio1 < 0; sig = neg & (A.beta_bio1_p < 0.05)
    print(f"    {len(A):,} insertions over {n0} sites; beta(bio1) < 0: {100*neg.mean():.1f}%; "
          f"beta < 0 and p < 0.05: {sig.sum():,} ({100*sig.mean():.1f}%); beta > 0 and p < 0.05: {((~neg) & (A.beta_bio1_p < 0.05)).sum():,}")
    A.to_csv(f"{R}/insertion_beta_all.csv", index=False)

    print("[7] the same for ALL deletions (context from _deletion_genomic_context.py)")
    dctx = pd.read_csv(f"{R}/deletion_context.csv")
    del_all = ~ins_all
    Sd = np.vstack([D0[f"{s}_s_SV"][del_all].astype(np.float64) for s in sites0]); Sdc = Sd - Sd.mean(0)
    Bd = dict(chrom=ch[sv_all][del_all], rec=rec[sv_all][del_all], p0=p0[sv_all][del_all],
              mean_s=Sd.mean(0), n_sites=n0)
    for k in ("bio1", "bio18"):
        c = cmeta.loc[sites0, k].to_numpy(); cc = c - c.mean(); sxx = (cc ** 2).sum()
        b = (cc[:, None] * Sdc).sum(0) / sxx
        se = np.sqrt(((Sdc - b[None, :] * cc[:, None]) ** 2).sum(0) / (n0 - 2) / sxx)
        Bd[f"beta_{k}"] = b; Bd[f"beta_{k}_se"] = se; Bd[f"beta_{k}_p"] = 2 * st.t.sf(np.abs(b / se), df=n0 - 2)
    Bd = pd.DataFrame(Bd).merge(dctx, on=["chrom", "rec"], how="left", validate="one_to_one")
    assert Bd.genic.notna().all() and len(Bd) == len(dctx)
    # tested (MAC>=12) deletions: flag A against p0-matched SNP beta, exactly as for insertions
    tested_del = pd.DataFrame(dict(chrom=ch[sv][~isins], rec=rec[sv][~isins], tested=True))
    Bd = Bd.merge(tested_del, on=["chrom", "rec"], how="left"); Bd["tested"] = Bd.tested.fillna(False).astype(bool)
    t = Bd.tested.values
    pct = pct_in_matched_snp(Bd.beta_bio1.values[t], Bd.p0.values[t], Z["beta_snp_bio1"], Z["p0_snp"])
    Bd["pct_beta1"] = np.nan; Bd.loc[t, "pct_beta1"] = pct
    Bd["flagA"] = Bd.pct_beta1 < TAIL
    negd = Bd.beta_bio1 < 0; sigd = negd & (Bd.beta_bio1_p < 0.05)
    print(f"    {len(Bd):,} deletions; beta(bio1) < 0: {100*negd.mean():.1f}%; beta < 0 and p < 0.05: {sigd.sum():,} "
          f"({100*sigd.mean():.1f}%); beta > 0 and p < 0.05: {((~negd) & (Bd.beta_bio1_p < 0.05)).sum():,}; "
          f"tested {t.sum()}, flagged vs matched SNPs {Bd.flagA.sum()}")
    Bd.to_csv(f"{R}/deletion_beta_all.csv", index=False)

    J.to_csv(f"{R}/insertion_landing_vs_purging.csv", index=False)
    with open(f"{R}/insertion_landing_vs_purging.json", "w") as f:
        json.dump(dict(min_p0=MIN_P0, sv_bp=SV_BP, min_mac=MIN_MAC, tail=TAIL, nbin=NBIN,
                       hot_bio1=HOT_BIO1, cold_bio1=COLD_BIO1, min_hot=MIN_HOT,
                       hot_sites=hot, cold_sites=cold,
                       site_bio1={int(s): float(b) for s, b in smeta.bio1.items()},
                       n_tested_sv=int(sv.size), n_tested_ins=int(isins.sum()),
                       n_all_insertions=int(len(ctx)), corr_p0_freq=float(r)), f, indent=1)
    print(f"\n[wrote] {R}/insertion_landing_vs_purging.{{csv,json}}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Insertion landing vs purging -- RECORD-level (allele) filter. Corrected `_insertion_landing_vs_purging.py`.

WHAT WAS WRONG. `_insertion_landing_vs_purging.py` took its tested set, the SV climate slopes and the
SNP reference from `s_climate_slope.npz` / `s_dist_by_stratum.npz`, all built with
`lib.founder_panel_keep()`: a POSITION-level MAC floor, so at a multi-allelic position every allele
passed if any sibling cleared it (bug found 2026-09-18, see lib.founder_panel_keep docstring and
sv_purging_fixed.ipynb). Rare alleles were then analysed as alleles -- matched on their own p0,
labelled from their own carriers -- after being admitted on the bubble's frequency. The same filter
also defined the SNP reference, so both sides of the "purged more than a matched SNP" comparison
were affected.

THE RULE HERE. One unit throughout, the panel record: `lib.founder_record_keep()` (own carriers
>= MIN_MAC, called in >= 90% of 231 founders) for SVs AND SNPs, then the same p0 window
(0.02..0.98). Everything the old script loaded from the position-filtered npz is recomputed here
from the pool AF matrices with the same estimator (`plot_replicate_sz`, plot-replicate logit
slope) and the same across-site OLS climate slope:

  flag A  beta(bio1) below the 10th percentile of p0-matched SNP betas
  flag B  per-site s below the p0-matched SNP 10th pct in >= 3 hot gardens (control: cold gardens)

Unchanged and reused (no MAC filter, so the bug never touched them): insertion_context.csv,
deletion_context.csv (per-record carrier counts, binary var_pa), insertion_beta_all.csv /
deletion_beta_all.csv betas (no-filter arm) -- only the deletions' `tested`/`flagA` columns
depend on the filter and are rewritten here.

CONTROL. SV_FILTER=position SV_MIN_MAC=12 reruns the OLD filter through THIS code; it must
reproduce the old tested set (2,546 insertions / 4,946 SVs) and the old betas, else nothing
downstream means anything. The script asserts it.

FULL SPECTRUM. SV_FILTER=none keeps only the own-carrier floor (SV_MIN_MAC, default 2: private /
singleton records are dropped as possibly spurious) and drops the call-rate floor and the p0 window:
every record with MIN_MAC..N-MIN_MAC own carriers is scored, SVs and SNP reference alike, so the
frequency spectra of the flagged set are not truncated. p0 <= EPS (1e-3) records share one clipped starting logit;
they are ranked only against SNPs in the same p0 bin, which share the clip.

Env: kmate. SV_FILTER=record|position|none (default record), SV_MIN_MAC (default 2), SV_SNP_REF
(SNPs in the beta reference, default 300000).
Writes results/sv_adaptive/insertion_landing_vs_purging{TAG}.{csv,json},
deletion_landing_vs_purging{TAG}.csv, s_climate_slope{TAG}.npz, TAG = _rec{MAC}, _pos{MAC} or _all{MAC}.
"""
import os, sys, json, glob
import numpy as np, pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import importlib.util
_sp = importlib.util.spec_from_file_location(
    "pm", os.path.join(os.path.dirname(os.path.abspath(__file__)), "_temporal_s_plots_snp_vs_nonsnp.py"))
pm = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(pm)

os.chdir("/global/scratch/users/tbellg/kmate")
R = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
STORE = lib.AF_STORE; PM = f"{lib.GEA}/common/results/pool_matrices"
MIN_P0, SV_BP, TAIL, NBIN = 0.02, 50, 0.10, 20
HOT_BIO1, COLD_BIO1, MIN_HOT = 17.0, 12.0, 3
FILTER = os.environ.get("SV_FILTER", "record")
MIN_MAC = int(os.environ.get("SV_MIN_MAC", 2))
N_SNP_REF = int(os.environ.get("SV_SNP_REF", 300000))
N_SNP_SITE = 25000                                   # per-site SNP s kept for flag B (as before)
if FILTER == "none":
    MIN_P0 = 0.0
TAG = f"_all{MIN_MAC}" if FILTER == "none" else f"_{'rec' if FILTER == 'record' else 'pos'}{MIN_MAC}"


def keep_mask(kind, chrom, pos):
    if FILTER == "none":      # own-carrier floor only: no call-rate floor, no p0 window
        return lib.founder_record_keep(kind, np.arange(len(pos)), min_mac=MIN_MAC, called_min=0.0)
    if FILTER == "record":
        return lib.founder_record_keep(kind, np.arange(len(pos)), min_mac=MIN_MAC)
    return lib.founder_panel_keep(chrom, pos, min_mac=MIN_MAC)


def pct_in_matched_snp(val_sv, p0_sv, val_snp, p0_snp, nbin=NBIN):
    """Percentile of each SV value within the SNP distribution of its p0 bin (as in the old script)."""
    e = np.quantile(np.concatenate([p0_sv, p0_snp]), np.linspace(0, 1, nbin + 1)); e[-1] += 1e-9
    bs = np.digitize(p0_sv, e[1:-1]); bn = np.digitize(p0_snp, e[1:-1])
    out = np.full(val_sv.size, np.nan)
    for b in range(nbin):
        ref = np.sort(val_snp[bn == b]); m = bs == b
        if ref.size and m.any():
            out[m] = np.searchsorted(ref, val_sv[m]) / ref.size
    return out


def own_carriers(rows):
    """own carrier count and call count per af_store non-SNP row (binary var_pa)."""
    idx = np.load(f"{STORE}/index_nonsnp.npz"); ch = idx["chrom"].astype("U5")
    na = np.zeros(rows.size, np.int64); nc = np.zeros(rows.size, np.int64); off = 0
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        nonsnp = np.where(~((meta["ref_len"] == 1) & (meta["alt_len"] == 1)))[0]
        n = int((ch == C).sum()); sel = (rows >= off) & (rows < off + n)
        if sel.any():
            c = nonsnp[rows[sel] - off]
            na[sel] = np.asarray((sp.load_npz(f"{base}.var_pa.npz").tocsc()[:, c] > 0).sum(0)).ravel()
            nc[sel] = np.asarray((sp.load_npz(f"{base}.var_called.npz").tocsc()[:, c] > 0).sum(0)).ravel()
        off += n
    return na, nc


def main():
    print(f"filter={FILTER} min_mac={MIN_MAC} -> tag {TAG}")
    idx_non = np.load(f"{STORE}/index_nonsnp.npz"); idx_snp = np.load(f"{STORE}/index_snp.npz")
    ch = idx_non["chrom"].astype("U5"); pos = idx_non["pos"].astype(np.int64)
    rl = idx_non["ref_len"].astype(np.int64); al = idx_non["alt_len"].astype(np.int64)
    dlen = np.abs(al - rl)
    p0_non = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    p0_snp = np.load(f"{STORE}/p0_snp.npy").astype(np.float64)

    print("[1] store row -> panel record ordinal (the context-table key)")
    rec = np.full(pos.size, -1, np.int64)
    for c in lib.CHROMS:
        meta = np.load(f"{lib.PROJ}/panel/arch3/{c.lower()}/var_pa_231_arch3_{c.lower()}.meta.npz", allow_pickle=True)
        mr = meta["ref_len"].astype(np.int64); ma = meta["alt_len"].astype(np.int64)
        w = np.where(~((mr == 1) & (ma == 1)))[0]; m = ch == c
        assert w.size == m.sum() and np.array_equal(meta["pos"].astype(np.int64)[w], pos[m]) \
            and np.array_equal(mr[w], rl[m]) and np.array_equal(ma[w], al[m])
        rec[m] = w

    print("[2] tested sets")
    k_non = keep_mask("nonsnp", ch, pos) & (p0_non >= MIN_P0) & (p0_non <= 1 - MIN_P0)
    sv = np.where(k_non & (dlen > SV_BP))[0]
    k_snp = keep_mask("snp", idx_snp["chrom"].astype("U5"), idx_snp["pos"].astype(np.int64)) \
        & (p0_snp >= MIN_P0) & (p0_snp <= 1 - MIN_P0)
    snp_all = np.where(k_snp)[0]
    rng = np.random.default_rng(0)
    snp = np.sort(rng.choice(snp_all, min(N_SNP_REF, snp_all.size), replace=False))
    isins = al[sv] > rl[sv]
    na, nc = own_carriers(sv)
    print(f"    {sv.size:,} SVs ({isins.sum():,} insertions, {(~isins).sum():,} deletions); "
          f"{snp_all.size:,} SNPs pass, {snp.size:,} in the reference sample")
    print(f"    own carriers of tested SVs: min {na.min()}, below floor {(na < MIN_MAC).sum()}; "
          f"call rate min {nc.min() / 231:.3f}")
    if FILTER in ("record", "none"):
        assert (na >= MIN_MAC).all(), "carrier floor did not hold"
    if FILTER == "record":
        assert (nc >= 0.9 * 231).all(), "call-rate floor did not hold"

    print("[3] per-site s (plot-replicate logit slope) for tested SVs + SNP reference")
    cmeta = pd.concat([pd.read_csv(f"{PM}/pool_gen{g}_snp.meta.csv") for g in (1, 2, 3)],
                      ignore_index=True).groupby("site")[["bio1", "bio18"]].mean()
    sites_all = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                        for s in pd.read_csv(f).site.unique()})
    p0v, p0s = p0_non[sv], p0_snp[snp]
    i_site = np.sort(rng.choice(snp.size, min(N_SNP_SITE, snp.size), replace=False))
    Ssv, Ssn, sites = [], [], []
    for site in sites_all:
        s_sn, _, _, npl = pm.plot_replicate_sz(site, "snp", snp, p0s)
        if npl < 2 or site not in cmeta.index:
            continue
        s_sv, _, _, _ = pm.plot_replicate_sz(site, "nonsnp", sv, p0v)
        Ssv.append(s_sv.astype(np.float64)); Ssn.append(s_sn.astype(np.float64)); sites.append(site)
        print(f"    site {site:>2} bio1={cmeta.loc[site, 'bio1']:>5.1f} plots={npl}", flush=True)
    Ssv = np.vstack(Ssv); Ssn = np.vstack(Ssn)                     # [sites x variants]
    clim = cmeta.loc[sites]

    def ols(S, k):
        c = clim[k].to_numpy(); cc = c - c.mean()
        return (cc[:, None] * (S - S.mean(0))).sum(0) / (cc ** 2).sum()
    B = {f"beta_{w}_{k}": ols(S, k) for w, S in (("sv", Ssv), ("snp", Ssn)) for k in ("bio1", "bio18")}

    if FILTER == "position" and MIN_MAC == 12:
        Z = np.load(f"{R}/s_climate_slope.npz")
        assert sv.size == Z["p0_sv"].size and np.allclose(p0v, Z["p0_sv"], atol=1e-6), "old SV set not reproduced"
        r = np.corrcoef(B["beta_sv_bio1"], Z["beta_sv_bio1"])[0, 1]
        print(f"    CONTROL: old SV set reproduced ({sv.size}); corr(beta_bio1 new, old) = {r:.5f}; "
              f"max |diff| = {np.abs(B['beta_sv_bio1'] - Z['beta_sv_bio1']).max():.2e}")
        assert r > 0.999, "old betas not reproduced"

    if FILTER == "none":
        A0 = pd.read_csv(f"{R}/insertion_beta_all.csv", usecols=["chrom", "rec", "beta_bio1"])
        m = pd.DataFrame(dict(chrom=ch[sv], rec=rec[sv], b=B["beta_sv_bio1"]))[isins].merge(A0, on=["chrom", "rec"])
        r = np.corrcoef(m.b, m.beta_bio1)[0, 1]
        print(f"    CHECK vs no-filter arm (insertion_beta_all.csv): {len(m):,} insertions matched, corr(beta_bio1) = {r:.5f}")

    np.savez_compressed(f"{R}/s_climate_slope{TAG}.npz", p0_snp=p0s.astype(np.float32),
                        p0_sv=p0v.astype(np.float32), sv_isdel=~isins, sites=np.array(sites),
                        mean_sv=Ssv.mean(0).astype(np.float32), mean_snp=Ssn.mean(0).astype(np.float32),
                        **{k: v.astype(np.float32) for k, v in B.items()})

    print("[4] insertions: join context by record ordinal, flag A and B")
    base = pd.DataFrame(dict(chrom=ch[sv], rec=rec[sv], afrow=sv, pos=pos[sv], p0=p0v,
                             n_carriers_own=na, call_rate=nc / 231, mean_s=Ssv.mean(0),
                             beta_bio1=B["beta_sv_bio1"], beta_bio18=B["beta_sv_bio18"]))
    base["pct_beta1"] = pct_in_matched_snp(base.beta_bio1.values, p0v, B["beta_snp_bio1"], p0s)
    base["pct_beta18"] = pct_in_matched_snp(base.beta_bio18.values, p0v, B["beta_snp_bio18"], p0s)
    base["flagA"] = base.pct_beta1 < TAIL
    hot = [s for s in sites if clim.loc[s, "bio1"] >= HOT_BIO1]
    cold = [s for s in sites if clim.loc[s, "bio1"] <= COLD_BIO1]
    n_hot = np.zeros(sv.size, int); n_cold = np.zeros(sv.size, int)
    for j, s in enumerate(sites):
        base[f"s_{s}"] = Ssv[j]
        pct = pct_in_matched_snp(Ssv[j], p0v, Ssn[j, i_site], p0s[i_site])
        n_hot += (pct < TAIL) & (s in hot); n_cold += (pct < TAIL) & (s in cold)
    base["n_hot_low"] = n_hot; base["n_cold_low"] = n_cold; base["flagB"] = n_hot >= MIN_HOT

    ctx = pd.read_csv(f"{R}/insertion_context.csv")
    J = base[isins].merge(ctx.drop(columns=["pos"]), on=["chrom", "rec"], how="left", validate="one_to_one")
    assert J.genic.notna().all()
    assert (J.n_alt.astype(int) == J.n_carriers_own).all(), "context-table carriers != own carriers"
    r = np.corrcoef(J.p0, J.freq)[0, 1]
    print(f"    {len(J):,} insertions; corr(p0, founder freq) = {r:.4f}; flag A {J.flagA.sum()} "
          f"({100 * J.flagA.mean():.1f}%); flag B {J.flagB.sum()}; A&B {(J.flagA & J.flagB).sum()}")

    print("[5] deletions: same tested set logic, context by record ordinal")
    dctx = pd.read_csv(f"{R}/deletion_context.csv")
    Dd = base[~isins].merge(dctx.drop(columns=["pos"]), on=["chrom", "rec"], how="left", validate="one_to_one")
    assert Dd.genic.notna().all() and (Dd.n_alt.astype(int) == Dd.n_carriers_own).all()
    print(f"    {len(Dd):,} deletions; flag A {Dd.flagA.sum()} ({100 * Dd.flagA.mean():.1f}%)")

    J.to_csv(f"{R}/insertion_landing_vs_purging{TAG}.csv", index=False)
    Dd.to_csv(f"{R}/deletion_landing_vs_purging{TAG}.csv", index=False)
    with open(f"{R}/insertion_landing_vs_purging{TAG}.json", "w") as f:
        json.dump(dict(filter=FILTER, min_mac=MIN_MAC, min_p0=MIN_P0, sv_bp=SV_BP, tail=TAIL, nbin=NBIN,
                       hot_bio1=HOT_BIO1, cold_bio1=COLD_BIO1, min_hot=MIN_HOT, hot_sites=hot,
                       cold_sites=cold, sites=sites, n_snp_pass=int(snp_all.size), n_snp_ref=int(snp.size),
                       n_tested_sv=int(sv.size), n_tested_ins=int(isins.sum()),
                       n_tested_del=int((~isins).sum()), corr_p0_freq=float(r)), f, indent=1)
    print(f"[wrote] {R}/insertion_landing_vs_purging{TAG}.{{csv,json}}, deletion_landing_vs_purging{TAG}.csv, "
          f"s_climate_slope{TAG}.npz")


if __name__ == "__main__":
    main()

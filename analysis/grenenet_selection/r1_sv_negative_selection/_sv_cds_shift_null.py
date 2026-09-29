#!/usr/bin/env python
"""Random-CDS control for the shared-CDS insertion result (events).

Is the result about GENES, or would any similarly sized, similarly clustered set of genomic intervals
give it? The whole TAIR10 annotation (CDS, UTR, exon, gene, TE, promoter layers) is circularly shifted
along each chromosome by an independent random offset (>= 1 Mb from 0), keeping every interval's size
and the annotation's internal structure; contexts are recomputed (same rules, first match wins) and
the two headline statistics are recomputed on the shifted "pseudo-CDS":

  1. BURDEN  mean log-fitness ~ z(shared ins-event rate in CDS) + z(... intergenic) + call rates +
             cactus + kinship  (as _sv_burden_lmm.py) -> coefficient on the CDS rate
  2. LOAD    per garden, log load change of shared CDS insertion events minus shared CDS SNPs matched
             on carrier count x mean carrier bio1 (as sv_burden_lmm.ipynb 4c) -> mean over gardens,
             gardens below 0
Empirical p = share of shifts at least as extreme as the real annotation (one-sided, negative). Also a
SIZE-MATCHED null: insertions avoid real CDS, so a shifted pseudo-CDS holds ~6x more events than the
real one; there the pseudo-CDS set is subsampled to the real CDS count before refitting.

Units: insertion EVENTS (EVENT_MAP=truvari|strict), >= 2 carriers, called >= 90%; SNPs a random 10%
(same draw as _sv_grenenet_sfs.py). Env: kmate. Writes results/sv_adaptive/cds_shift_null{TAG}.npz.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import importlib.util
_H = os.path.dirname(os.path.abspath(__file__))
def _load(name, fn):
    s = importlib.util.spec_from_file_location(name, os.path.join(_H, fn)); m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m); return m
sc = _load("sfsctx", "_sv_sfs_context.py")
lm = _load("lmm", "_sv_burden_lmm.py")

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
EF = f"{lib.GEA}/r3_persite_gwas/results/ecotype_fitness/founder_site_dh.npz"
EVMAP = os.environ.get("EVENT_MAP", "truvari")
TAG = "_eventstrict" if EVMAP == "strict" else "_event"
NSHIFT, SV_BP, CALL_MIN, MIN_CAR, SNP_KEEP, EPS = int(os.environ.get("NSHIFT", 200)), 50, 0.9, 2, 0.10, 1e-5
CHLEN = {"Chr1": 30427671, "Chr2": 19698289, "Chr3": 23459830, "Chr4": 18585056, "Chr5": 26975502}
CB = np.array([2, 3, 4, 5, 6, 8, 11, 16, 26, 51, 1000])
iCDS, iINT = sc.CTX.index("CDS"), sc.CTX.index("intergenic")


def shift_layers(IV, offs):
    out = {}
    for k, d in IV.items():
        nd = {}
        for c, (s, e) in d.items():
            L = CHLEN[c]; o = offs[c]
            s2 = (s + o - 1) % L + 1; e2 = (e + o - 1) % L + 1
            w = e2 < s2
            S = np.r_[s2[~w], s2[w], np.ones(w.sum(), np.int64)]
            E = np.r_[e2[~w], np.full(w.sum(), L), e2[w]]
            o_ = np.argsort(S, kind="stable"); nd[c] = (S[o_], E[o_])
        out[k] = nd
    return out


def main():
    rng = np.random.default_rng(0)
    FO = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    GC = (FO.group == "cold").to_numpy(); GW = (FO.group == "warm").to_numpy(); B1 = FO.bio1.to_numpy(np.float64)
    N = len(FO)
    Z = np.load(EF, allow_pickle=True); assert np.array_equal(Z["founders"].astype(str), FO.founder.to_numpy())
    y = np.log((Z["H0"][None, :] + Z["DH_rel"] + EPS) / (Z["H0"][None, :] + EPS)).mean(0)
    H0 = Z["H0"] / Z["H0"].sum()
    HL = np.clip(Z["H0"][None, :] + Z["DH_rel"], 0, None); HL /= HL.sum(1, keepdims=True)
    Hall = np.vstack([H0[None, :], HL])                                  # [32 x founders]
    D = pd.read_csv(f"{OUT}/founder_burden{TAG}.csv", dtype={"founder": str})
    cov = np.column_stack([np.ones(N), lm.z(D.callrate_sv), lm.z(D.callrate_snp), D.is_cactus.astype(float)])
    L = lm.LMM(np.load(f"{OUT}/founder_kinship.npy"))
    EM = np.load(f"{OUT}/arch3_event_map{'_strict' if EVMAP == 'strict' else ''}.npz")

    ev = {k: [] for k in ("ch", "lo", "hi", "n", "b1", "p")}; EP, EQ = [], []
    sn = {k: [] for k in ("ch", "lo", "hi", "n", "b1", "p")}
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        # SNPs: same 10% draw as _sv_grenenet_sfs.py (its rng is consumed only here, in this order)
        m = (rl == 1) & (al == 1) & (rng.random(pos.size) < SNP_KEEP)
        cols = np.where(m)[0]
        P = (vp[:, cols] > 0).tocsc(); Q = (vc[:, cols] > 0).tocsc()
        c = np.asarray(P.sum(0)).ravel(); k = np.asarray(Q.sum(0)).ravel()
        nc = np.asarray(P[GC].sum(0)).ravel(); nw = np.asarray(P[GW].sum(0)).ravel()
        ok = np.where((k >= CALL_MIN * N) & (c >= MIN_CAR))[0]
        Pf = P[:, ok].astype(np.float64); Qf = Q[:, ok].astype(np.float64)
        b1 = np.asarray(Pf.T @ B1).ravel() / c[ok]
        shared = (nc[ok] > 0) & (nw[ok] > 0)
        p = (np.asarray((sp.csr_matrix(Hall) @ Pf).todense()) / np.asarray((sp.csr_matrix(Hall) @ Qf).todense())).T
        sn["ch"].append(np.full(ok.size, C)); sn["lo"].append(pos[cols][ok]); sn["hi"].append(pos[cols][ok])
        sn["n"].append(c[ok]); sn["b1"].append(b1); sn["p"].append(p)
        sn.setdefault("shared", []).append(shared)
        # insertion events (union over alleles), context from the first record, as _sv_burden_build.py
        rec, bub = EM[f"{C}_rec"], EM[f"{C}_bubble"]
        ins = (al[rec] - rl[rec]) > SV_BP; rec, bub = rec[ins], bub[ins]
        ub, li = np.unique(bub, return_inverse=True); li = li.ravel()
        M = sp.csr_matrix((np.ones(rec.size), (li, np.arange(rec.size))), shape=(ub.size, rec.size))
        UP = ((M @ (vp[:, rec] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64)   # events x founders
        UQ = ((M @ (vc[:, rec] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64)
        un = np.asarray(UP.sum(1)).ravel(); uk = np.asarray(UQ.sum(1)).ravel()
        unc = np.asarray(UP[:, GC].sum(1)).ravel(); unw = np.asarray(UP[:, GW].sum(1)).ravel()
        first = np.full(ub.size, -1, np.int64)
        for j in range(rec.size - 1, -1, -1):
            first[li[j]] = j
        keep = np.where((uk >= CALL_MIN * N) & (un >= MIN_CAR) & (unc > 0) & (unw > 0))[0]   # shared only
        UP, UQ = UP.tocsr()[keep], UQ.tocsr()[keep]
        fr = rec[first[keep]]
        ev["ch"].append(np.full(keep.size, C)); ev["lo"].append(pos[fr]); ev["hi"].append(pos[fr] + np.maximum(rl[fr] - 1, 0))
        ev["n"].append(un[keep]); ev["b1"].append(np.asarray(UP @ B1).ravel() / un[keep])
        ev["p"].append((np.asarray(UP @ Hall.T) / np.asarray(UQ @ Hall.T)))
        EP.append(UP); EQ.append(UQ)
        print(f"[{C}] shared insertion events {keep.size:,}; SNPs {ok.size:,}", flush=True)
    E = {k: np.concatenate(v) for k, v in ev.items()}; S = {k: np.concatenate(v) for k, v in sn.items()}
    EP = sp.vstack(EP).tocsr(); EQ = sp.vstack(EQ).tocsr()
    BQ = np.quantile(S["b1"], [0.2, 0.4, 0.6, 0.8])
    Ecell = np.digitize(E["n"], CB) * 10 + np.digitize(E["b1"], BQ)
    sh = S["shared"]
    Scell = (np.digitize(S["n"], CB) * 10 + np.digitize(S["b1"], BQ))[sh]
    Sp, Sch, Slo, Shi = S["p"][sh], S["ch"][sh], S["lo"][sh], S["hi"][sh]
    IV0 = sc.layers()

    def stats_for(IV, nsub=None, rs=None):
        ce = sc.context(IV, E["ch"], E["lo"], E["hi"]); cs = sc.context(IV, Sch, Slo, Shi)
        # 1. burden
        def rate(msk):
            num = np.asarray(EP[msk].sum(0)).ravel(); den = np.asarray(EQ[msk].sum(0)).ravel()
            return num / np.maximum(den, 1)
        mcds = ce == iCDS
        if nsub is not None and mcds.sum() > nsub:
            idx = np.where(mcds)[0]; mcds = np.zeros_like(mcds); mcds[rs.choice(idx, nsub, replace=False)] = True
        t, r = rate(mcds), rate(ce == iINT)
        b = L.fit(y, np.column_stack([cov, lm.z(t), lm.z(r)]))[0][-2] if t.std() > 0 else np.nan
        # 2. load, shared CDS events vs matched shared CDS SNPs
        ms, mn = ce == iCDS, cs == iCDS
        if nsub is not None and ms.sum() > nsub:                  # size-matched null: same number of
            idx = np.where(ms)[0]                                  # pseudo-CDS events as the real CDS set
            ms = np.zeros_like(ms); ms[rs.choice(idx, nsub, replace=False)] = True
        cse = pd.Series(Ecell[ms]).value_counts(); csn = pd.Series(Scell[mn]).value_counts()
        kp = np.isin(Ecell[ms], csn.index)
        Ls = E["p"][ms][kp].sum(0)
        w = (cse / csn).reindex(Scell[mn]).fillna(0).to_numpy()
        Ln = (Sp[mn] * w[:, None]).sum(0)
        d = np.log(Ls[1:] / Ls[0]) - np.log(Ln[1:] / Ln[0])
        return b, d.mean(), int((d < 0).sum()), int(mcds.sum())

    real = stats_for(IV0)
    print(f"\nREAL annotation: burden coef {real[0]:+.4f}; load diff {real[1]:+.4f}, below 0 in {real[2]}/31; "
          f"shared CDS insertion events {real[3]}", flush=True)
    rs = np.random.default_rng(1); NB, ND, NBL, NN, MB, MD = [], [], [], [], [], []
    nreal = real[3]
    for i in range(NSHIFT):
        offs = {c: int(rs.integers(1_000_000, CHLEN[c] - 1_000_000)) for c in lib.CHROMS}
        IVs = shift_layers(IV0, offs)
        b, dm, nb, nn = stats_for(IVs)
        mb, md, _, _ = stats_for(IVs, nsub=nreal, rs=rs)          # same n as the real CDS set
        NB.append(b); ND.append(dm); NBL.append(nb); NN.append(nn); MB.append(mb); MD.append(md)
        if (i + 1) % 25 == 0:
            print(f"  {i + 1} shifts", flush=True)
    NB, ND, NBL, NN, MB, MD = map(np.array, (NB, ND, NBL, NN, MB, MD))
    pb = (np.sum(NB <= real[0]) + 1) / (NSHIFT + 1); pd_ = (np.sum(ND <= real[1]) + 1) / (NSHIFT + 1)
    print(f"\nburden coef: real {real[0]:+.4f}; shifts median {np.nanmedian(NB):+.4f} [5-95%: {np.nanpercentile(NB,5):+.4f}, "
          f"{np.nanpercentile(NB,95):+.4f}]; empirical p = {pb:.3f}")
    print(f"load diff:   real {real[1]:+.4f}; shifts median {np.median(ND):+.4f} [5-95%: {np.percentile(ND,5):+.4f}, "
          f"{np.percentile(ND,95):+.4f}]; empirical p = {pd_:.3f}")
    pbm = (np.sum(MB <= real[0]) + 1) / (NSHIFT + 1); pdm = (np.sum(MD <= real[1]) + 1) / (NSHIFT + 1)
    print(f"pseudo-CDS size: real {real[3]} events; shifts median {int(np.median(NN))}")
    print(f"SIZE-MATCHED null (pseudo-CDS subsampled to {real[3]} events):")
    print(f"  burden coef: shifts median {np.nanmedian(MB):+.4f} [5-95%: {np.nanpercentile(MB,5):+.4f}, "
          f"{np.nanpercentile(MB,95):+.4f}]; empirical p = {pbm:.3f}")
    print(f"  load diff:   shifts median {np.median(MD):+.4f} [5-95%: {np.percentile(MD,5):+.4f}, "
          f"{np.percentile(MD,95):+.4f}]; empirical p = {pdm:.3f}")
    np.savez(f"{OUT}/cds_shift_null{TAG}.npz", real=np.array(real[:3], float), null_burden=NB, null_load=ND,
             null_below=NBL, null_n=NN, p_burden=pb, p_load=pd_,
             null_burden_matched=MB, null_load_matched=MD, p_burden_matched=pbm, p_load_matched=pdm)
    print(f"[wrote] {OUT}/cds_shift_null{TAG}.npz")


if __name__ == "__main__":
    main()

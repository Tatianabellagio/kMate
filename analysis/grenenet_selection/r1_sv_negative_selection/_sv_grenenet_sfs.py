#!/usr/bin/env python
"""Site-frequency spectra IN GrENE-Net: variant frequencies at the start (seed mix) and in the last
generation of each garden, from founder composition x founder genotypes.

  p_v(t) = sum_f h_f(t) g_fv / sum_f h_f(t) c_fv
  h(0)   = H0 (founding composition); h(last, garden) = H0 + DH_rel  (founder_site_dh.npz, relative
           plot weights; negatives clipped to 0, renormalised)
  g_fv   = founder f carries v; c_fv = founder f is called at v (missing calls drop out of the mean)

Variants as in `_sv_burden_build.py`: called in >= 90% of 231 founders, >= 2 carriers; SNPs and short
indels per record (SNPs a random 10%), insertions / deletions (> 50 bp) per bubble; origin
(shared / cold / warm) and context (CDS ... intergenic) labels as there.
Also the population LOAD of a class: L(t) = sum_v p_v(t) = expected number of the class's variants
carried per genome, at the start and per garden.
SV_UNIT=allele: insertions / deletions per record, filtered on own carriers; SV_UNIT=event: Truvari-collapsed
alleles (`_build_event_map.py`); default bubble. Also stores each variant's mean carrier bio1 (`b1`).
Env: kmate. Writes results/sv_adaptive/grenenet_sfs{,_allele}.npz.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import importlib.util
_H = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("sfsctx", os.path.join(_H, "_sv_sfs_context.py"))
sc = importlib.util.module_from_spec(_s); _s.loader.exec_module(sc)

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
EF = f"{lib.GEA}/r3_persite_gwas/results/ecotype_fitness/founder_site_dh.npz"
SV_BP, MIN_CAR, SNP_KEEP, SEED = 50, 2, 0.10, 0
CALL_MIN = float(os.environ.get("SV_CALL_MIN", 0.9))   # founder call-rate floor; != 0.9 -> outputs tagged _cr{pct}
CRTAG = "" if CALL_MIN == 0.9 else f"_cr{int(round(CALL_MIN * 100))}"
UNIT = os.environ.get("SV_UNIT", "bubble")          # bubble | allele | event
TAG = {"bubble": "", "allele": "_allele", "event": "_event"}[UNIT]
EVMAP = os.environ.get("EVENT_MAP", "truvari")    # truvari | strict (sensitivity: near-identical alleles only)
MAPF = (f"arch3_event_map{'_strict' if EVMAP == 'strict' else ''}.npz" if UNIT == "event" else "arch3_bubble_map.npz")
if UNIT == "event" and EVMAP == "strict":
    TAG = "_eventstrict"
TAG += CRTAG


def main():
    rng = np.random.default_rng(SEED)
    FO = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    GC = (FO.group == "cold").to_numpy(); GW = (FO.group == "warm").to_numpy(); N = len(FO)
    Z = np.load(EF, allow_pickle=True)
    assert np.array_equal(Z["founders"].astype(str), FO.founder.to_numpy())
    H0 = Z["H0"] / Z["H0"].sum()
    HL = np.clip(Z["H0"][None, :] + Z["DH_rel"], 0, None); HL /= HL.sum(1, keepdims=True)
    Hall = np.vstack([H0[None, :], HL])                       # [1 + 31 gardens] x founders
    IV = sc.layers(); BM = np.load(f"{OUT}/{MAPF}")
    B1 = FO.bio1.to_numpy()
    F = {k: [] for k in ("kind", "nc", "nw", "n", "ctx", "p", "b1")}

    def put(kind, P, Q, c, nc, nw, cx):
        # P, Q: founders x variants (sparse bool/float)
        num = np.asarray((sp.csr_matrix(Hall) @ P.astype(np.float64)).todense())
        den = np.asarray((sp.csr_matrix(Hall) @ Q.astype(np.float64)).todense())
        F["p"].append((num / den).T.astype(np.float32))
        F["kind"].append(np.full(c.size, kind)); F["n"].append(c.astype(np.int16))
        F["nc"].append(nc.astype(np.int16)); F["nw"].append(nw.astype(np.int16)); F["ctx"].append(cx.astype(np.int8))
        F["b1"].append((np.asarray(P.astype(np.float64).T @ B1).ravel() / np.maximum(c, 1)).astype(np.float32))   # mean carrier bio1

    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64)
        al = meta["alt_len"].astype(np.int64); dl = al - rl
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        ctx = sc.context(IV, np.full(pos.size, C, dtype=object), pos, pos + np.maximum(rl - 1, 0))
        for kind, m in (("snp", (rl == 1) & (al == 1) & (rng.random(pos.size) < SNP_KEEP)),
                        ("indel", (np.abs(dl) >= 1) & (np.abs(dl) <= SV_BP) & ~((rl == 1) & (al == 1)))):
            cols = np.where(m)[0]
            P = (vp[:, cols] > 0).tocsc(); Q = (vc[:, cols] > 0).tocsc()
            c = np.asarray(P.sum(0)).ravel(); k = np.asarray(Q.sum(0)).ravel()
            ok = np.where((k >= CALL_MIN * N) & (c >= MIN_CAR))[0]
            P, Q = P[:, ok], Q[:, ok]
            put(kind, P, Q, c[ok], np.asarray(P[GC].sum(0)).ravel(), np.asarray(P[GW].sum(0)).ravel(), ctx[cols][ok])
        if UNIT == "allele":
            for kind, m in (("ins", dl > SV_BP), ("del", -dl > SV_BP)):
                cols = np.where(m)[0]
                P = (vp[:, cols] > 0).tocsc(); Q = (vc[:, cols] > 0).tocsc()
                c = np.asarray(P.sum(0)).ravel(); k = np.asarray(Q.sum(0)).ravel()
                ok = np.where((k >= CALL_MIN * N) & (c >= MIN_CAR))[0]
                P, Q = P[:, ok], Q[:, ok]
                put(kind, P, Q, c[ok], np.asarray(P[GC].sum(0)).ravel(), np.asarray(P[GW].sum(0)).ravel(), ctx[cols][ok])
            print(f"[{C}] done", flush=True)
            continue
        rec2b = pd.Series(BM[f"{C}_bubble"], index=BM[f"{C}_rec"])
        for kind, m in (("ins", dl > SV_BP), ("del", -dl > SV_BP)):
            cols = np.where(m)[0]
            b = rec2b.reindex(cols).to_numpy().astype(np.int64)
            ub, li = np.unique(b, return_inverse=True); li = li.ravel()
            M = sp.csr_matrix((np.ones(cols.size), (li, np.arange(cols.size))), shape=(ub.size, cols.size))
            Pf = ((M @ (vp[:, cols] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64).T.tocsc()  # founders x bubbles
            Qf = ((M @ (vc[:, cols] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64).T.tocsc()
            c = np.asarray(Pf.sum(0)).ravel(); k = np.asarray(Qf.sum(0)).ravel()
            first = np.full(ub.size, -1, np.int64)
            for j in range(cols.size - 1, -1, -1):
                first[li[j]] = j
            bctx = ctx[cols[first]]
            ok = np.where((k >= CALL_MIN * N) & (c >= MIN_CAR))[0]
            P, Q = Pf[:, ok], Qf[:, ok]
            put(kind, P, Q, c[ok], np.asarray(P[GC].sum(0)).ravel(), np.asarray(P[GW].sum(0)).ravel(), bctx[ok])
        print(f"[{C}] done", flush=True)
    S = {k: np.concatenate(v) for k, v in F.items()}
    S["p"] = np.vstack(S["p"])                              # variants x (1 + 31)
    S["sites"] = Z["sites"]; S["bio1"] = Z["bio1"]; S["zone"] = Z["zone"].astype(str)
    np.savez_compressed(f"{OUT}/grenenet_sfs{TAG}.npz", **S)
    for k in ("snp", "indel", "ins", "del"):
        print(f"{k:<6} n={(S['kind'] == k).sum():,}")
    print(f"[wrote] {OUT}/grenenet_sfs{TAG}.npz")


if __name__ == "__main__":
    main()

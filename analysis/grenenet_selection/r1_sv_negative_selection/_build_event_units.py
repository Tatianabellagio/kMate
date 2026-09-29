#!/usr/bin/env python
"""Units for the temporal (trajectory) analysis on the af_store NON-SNP index, at the EVENT level.

  short indels (|dlen| <= 50) : one unit per record (as before)
  SVs (|dlen| > 50)           : one unit per Truvari EVENT (`arch3_event_map.npz`): all alleles of
                                a pangenome bubble collapsed with Truvari's default matching.
Per unit: member af_store rows (CSR), representative row (the member with most carriers -- gives
chrom / pos / size / context), union carriers (founder carries any member), union called (called at
any member), carriers in the cold / warm origin tertiles, mean carrier bio1, and kind.
The unit's pool frequency is the SUM of its members' frequencies (members are carried by disjoint
founders: only 0.61% of founder-event carrier pairs carry two alleles of one event), capped at 1.

Env: kmate. Writes results/sv_adaptive/nonsnp_event_units.npz.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SV_BP = 50


def main():
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    GC = (F.group == "cold").to_numpy(); GW = (F.group == "warm").to_numpy(); B1 = F.bio1.to_numpy(np.float64)
    EM = np.load(f"{OUT}/arch3_event_map.npz")
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz"); ch = idx["chrom"].astype("U5")
    rl = idx["ref_len"].astype(np.int64); al = idx["alt_len"].astype(np.int64)
    U = {k: [] for k in ("rep", "n", "ncal", "nc", "nw", "b1", "dlen", "chrom", "pos")}
    members = []; off = 0
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64); pal = meta["alt_len"].astype(np.int64); ppos = meta["pos"].astype(np.int64)
        panel = np.where(~((prl == 1) & (pal == 1)))[0]                     # non-SNP panel columns, af_store order
        n = int((ch == C).sum()); assert panel.size == n
        assert np.array_equal(prl[panel], rl[off:off + n]) and np.array_equal(pal[panel], al[off:off + n])
        ev = pd.Series(EM[f"{C}_bubble"], index=EM[f"{C}_rec"]).reindex(panel).to_numpy()   # NaN = not an SV
        dl = pal[panel] - prl[panel]
        sv = np.abs(dl) > SV_BP
        assert np.isnan(ev[~sv]).all() and not np.isnan(ev[sv]).any()
        key = np.where(sv, ev, -1 - np.arange(n)).astype(np.int64)          # indels: own unit
        uk, li = np.unique(key, return_inverse=True); li = li.ravel()
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        P = (vp[:, panel] > 0).astype(np.float64).T.tocsr(); Q = (vc[:, panel] > 0).astype(np.float64).T.tocsr()
        M = sp.csr_matrix((np.ones(n), (li, np.arange(n))), shape=(uk.size, n))
        UP = ((M @ P) > 0).astype(np.float64); UQ = ((M @ Q) > 0).astype(np.float64)   # units x founders
        un = np.asarray(UP.sum(1)).ravel()
        ncar = np.asarray(P.sum(1)).ravel()
        # representative member: most carriers (ties -> first)
        order = np.lexsort((np.arange(n), -ncar, li))
        first = np.r_[0, np.flatnonzero(np.diff(li[order])) + 1]
        rep_loc = order[first]
        U["rep"].append(off + rep_loc); U["n"].append(un)
        U["ncal"].append(np.asarray(UQ.sum(1)).ravel())
        U["nc"].append(np.asarray(UP[:, GC].sum(1)).ravel()); U["nw"].append(np.asarray(UP[:, GW].sum(1)).ravel())
        U["b1"].append(np.asarray(UP @ B1).ravel() / np.maximum(un, 1))
        U["dlen"].append(dl[rep_loc]); U["chrom"].append(np.full(uk.size, C)); U["pos"].append(ppos[panel][rep_loc])
        grp = np.argsort(li, kind="stable")                                 # members of unit u, in row order
        members += [off + g for g in np.split(grp, np.cumsum(np.bincount(li, minlength=uk.size))[:-1])]
        print(f"[{C}] {n:,} non-SNP records -> {uk.size:,} units ({int(sv.sum()):,} SV records -> "
              f"{np.unique(ev[sv]).size:,} events)", flush=True)
        off += n
    S = {k: np.concatenate(v) for k, v in U.items()}
    S["indptr"] = np.r_[0, np.cumsum([m.size for m in members])].astype(np.int64)
    S["indices"] = np.concatenate(members).astype(np.int64)
    S["n_founders"] = np.int64(len(F))
    np.savez_compressed(f"{OUT}/nonsnp_event_units.npz", **S)
    k = np.abs(S["dlen"]) > SV_BP
    print(f"units: {S['rep'].size:,} ({int(k.sum()):,} SV events: {int((S['dlen'] > SV_BP).sum()):,} ins, "
          f"{int((S['dlen'] < -SV_BP).sum()):,} del)")
    print(f"[wrote] {OUT}/nonsnp_event_units.npz")


if __name__ == "__main__":
    main()

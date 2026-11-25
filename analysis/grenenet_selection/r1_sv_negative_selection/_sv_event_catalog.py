#!/usr/bin/env python
"""Catalog of SV EVENTS (`arch3_event_map.npz`, truvari collapse defaults) for the key-results notebook.

One row per insertion / deletion event with >= 1 carrier: position, size, genomic context and gene
(from the event's first record, as `_sv_burden_build.py`), carriers (union over alleles), call rate,
origin label, mean carrier bio1, inserted-sequence class of the first record
(`insertion_tair10_class.csv`; insertions only), n alleles merged, and -- for events with >= 2
carriers called in >= 90% of founders -- the frequency in the GrENE-Net population at the start and in
the last generation of each garden (as `_sv_grenenet_sfs.py`).

Env: kmate. Writes results/sv_adaptive/sv_event_catalog.csv, sv_event_catalog_p.npy (rows aligned).
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import importlib.util
_H = os.path.dirname(os.path.abspath(__file__))
_sa = importlib.util.spec_from_file_location("ann", os.path.join(_H, "_sv_purging_annotate.py"))
ann = importlib.util.module_from_spec(_sa); _sa.loader.exec_module(ann)
_s = importlib.util.spec_from_file_location("sfsctx", os.path.join(_H, "_sv_sfs_context.py"))
sc = importlib.util.module_from_spec(_s); _s.loader.exec_module(sc)

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
EF = f"{lib.GEA}/r3_persite_gwas/results/ecotype_fitness/founder_site_dh.npz"
SV_BP, PROM = 50, 2000


def main():
    FO = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    GC = (FO.group == "cold").to_numpy(); GW = (FO.group == "warm").to_numpy(); B1 = FO.bio1.to_numpy(np.float64)
    N = len(FO)
    Z = np.load(EF, allow_pickle=True); assert np.array_equal(Z["founders"].astype(str), FO.founder.to_numpy())
    H0 = Z["H0"] / Z["H0"].sum()
    HL = np.clip(Z["H0"][None, :] + Z["DH_rel"], 0, None); HL /= HL.sum(1, keepdims=True)
    Hall = np.vstack([H0[None, :], HL])
    IV = sc.layers(); EM = np.load(f"{OUT}/arch3_event_map.npz")
    genes = lib.load_genes()
    rows, PP = [], []
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        _ch = np.full(pos.size, C, dtype=object); _hi = pos + np.maximum(rl - 1, 0)
        ctx = sc.context(IV, _ch, pos, _hi)
        # EVERY layer the record overlaps, not only the cascade winner. `context` keeps one label; these
        # keep the rest, so the arbitrariness of the precedence rule is measurable on the model's own
        # event set instead of having to be reconstructed from a different table.
        _H = {k: ann.hit(v, _ch, pos, _hi) for k, v in IV.items()}
        LYR = {"L_CDS": _H["CDS"], "L_UTR": _H["UTR"], "L_intron": _H["gene"] & ~_H["exon"],
               "L_TE": _H["TE"], "L_promoter": _H["promoter"]}
        rec, ev = EM[f"{C}_rec"], EM[f"{C}_bubble"]
        ue, li = np.unique(ev, return_inverse=True); li = li.ravel()
        M = sp.csr_matrix((np.ones(rec.size), (li, np.arange(rec.size))), shape=(ue.size, rec.size))
        UP = ((M @ (vp[:, rec] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64).tocsr()
        UQ = ((M @ (vc[:, rec] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64).tocsr()
        first = np.full(ue.size, -1, np.int64)
        for j in range(rec.size - 1, -1, -1):
            first[li[j]] = j
        fr = rec[first]
        n = np.asarray(UP.sum(1)).ravel(); k = np.asarray(UQ.sum(1)).ravel()
        nc = np.asarray(UP[:, GC].sum(1)).ravel(); nw = np.asarray(UP[:, GW].sum(1)).ravel()
        g = genes[genes.chrom == C].sort_values("start")
        gs, ge, gid, gst = g.start.to_numpy(), g.end.to_numpy(), g.gene.to_numpy(), g.strand.to_numpy()
        lo = pos[fr]; hi = lo + np.maximum(rl[fr] - 1, 0)
        gene = np.full(ue.size, "", dtype=object)
        for i in range(ue.size):                                   # overlapping gene, else promoter gene
            j = np.searchsorted(gs, hi[i], side="right")
            cand = np.where(ge[max(0, j - 50):j] >= lo[i])[0]
            if cand.size:
                gene[i] = gid[max(0, j - 50) + cand[-1]]; continue
            w = np.where(((gst == "+") & (gs - PROM <= hi[i]) & (gs - 1 >= lo[i])) |
                         ((gst == "-") & (ge + 1 <= hi[i]) & (ge + PROM >= lo[i])))[0]
            if w.size:
                gene[i] = gid[w[0]]
        p = np.asarray(UP @ Hall.T) / np.maximum(np.asarray(UQ @ Hall.T), 1e-12)
        dl = al[fr] - rl[fr]
        rows.append(pd.DataFrame(dict(
            chrom=C, pos=lo, kind=np.where(dl > 0, "ins", "del"), size_bp=np.abs(dl), n_alleles=np.bincount(li),
            n_carriers=n.astype(int), n_called=k.astype(int), n_cold=nc.astype(int), n_warm=nw.astype(int),
            carrier_bio1=np.asarray(UP @ B1).ravel() / np.maximum(n, 1),
            context=np.array(sc.CTX)[ctx[fr]], gene=gene, key=[f"{C}|{r}" for r in fr],
            **{k: v[fr] for k, v in LYR.items()})))
        PP.append(p.astype(np.float32))
        print(f"[{C}] {ue.size:,} events", flush=True)
    D = pd.concat(rows, ignore_index=True); P = np.vstack(PP)
    T = pd.read_csv(f"{OUT}/insertion_tair10_class.csv")[["key", "cls", "family"]].rename(
        columns={"cls": "seq_class", "family": "te_family"})
    D = D.merge(T, on="key", how="left")
    D["origin"] = np.where((D.n_cold > 0) & (D.n_warm > 0), "shared",
                   np.where((D.n_cold >= 2) & (D.n_warm == 0), "cold-specific",
                    np.where((D.n_warm >= 2) & (D.n_cold == 0), "warm-specific", "other")))
    D["passes"] = (D.n_carriers >= 2) & (D.n_called >= 0.9 * N)
    keep = (D.n_carriers >= 1).to_numpy()
    D = D[keep].reset_index(drop=True); P = P[keep]
    D.to_csv(f"{OUT}/sv_event_catalog.csv", index=False); np.save(f"{OUT}/sv_event_catalog_p.npy", P)
    print(pd.crosstab(D[D.passes].kind, D[D.passes].context).to_string())
    print(f"[wrote] {OUT}/sv_event_catalog.csv ({len(D):,} events; {int(D.passes.sum()):,} pass the filter)")


if __name__ == "__main__":
    main()

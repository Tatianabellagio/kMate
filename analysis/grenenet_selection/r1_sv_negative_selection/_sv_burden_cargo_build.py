#!/usr/bin/env python
"""Founder SV burdens by CONTENT x LOCATION, insertions and deletions pooled -- for the cargo burden model.

The reported burden model (`_sv_burden_final_models.R`, SV6) splits SVs by direction and by ONE context from
the cascade CDS > UTR > intron > TE > promoter > intergenic. That cascade treats "inside a reference TE" as a
LOCATION and ranks it above promoter, so a deletion removing a TE that sits in a promoter is called `TE`,
and what an SV brings or removes is not represented at all.

Here every SV event gets two independent labels:
  content   what the event adds (insertion: its ALT) or removes (deletion: its REF) --
            `TE` = at least half TE (te_full + te_most of the GraffiTE-protocol classification,
            graffite/sv_repeat_class_events.csv), `nonTE` = everything else
  location  where it is, by gene structure ONLY: CDS > UTR > intron > promoter (<= 1 kb upstream of a TSS,
            as the catalog's L_promoter) > intergenic, first match wins; the reference-TE layer is not used
and insertions and deletions are POOLED within a class. The burden counts NON-Col-0 alleles carried (an
insertion carrier has the extra sequence, a deletion carrier lacks Col-0's), as the reported model does, so
`TE_promoter` = promoter TE presence/absence sites where the founder differs from Col-0.

Everything else is the reported build (`_sv_burden_build.py`, SV_UNIT=event): Truvari events
(`arch3_event_map.npz`), carriers = union over the event's alleles, called = called at any allele, filter
called in >= 90% of founders and >= 2 carriers, origin `shared` (>= 1 cold- and >= 1 warm-origin carrier) and
`all`. An event's labels come from its FIRST record (lowest panel index), as the reported build takes its
context. Per-direction versions of the same classes (`ins_TE_promoter`, `del_TE_promoter`, ...) are written too,
for comparison.

Columns: <class>_<origin> (carried) and called_<class>_<origin>, class in {TE,nonTE}_{CDS,UTR,intron,promoter,
intergenic} and the ins_/del_ versions. Env: kmate. Compute node.
Writes results/sv_adaptive/founder_burden_cargo.csv and prints the event count per class.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SV_BP, MIN_CAR, CALL_MIN = 50, 2, 0.9
LOCS = ["CDS", "UTR", "intron", "promoter", "intergenic"]
CONT = ["TE", "nonTE"]


def main():
    FO = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    GC = (FO.group == "cold").to_numpy(); GW = (FO.group == "warm").to_numpy(); N = len(FO)
    EM = np.load(f"{OUT}/arch3_event_map.npz")
    CAT = pd.read_csv(f"{OUT}/sv_event_catalog.csv", usecols=["key", "kind", "L_CDS", "L_UTR", "L_intron",
                                                                "L_promoter", "passes"]).set_index("key")
    EVC = pd.read_csv(f"{OUT}/graffite/sv_repeat_class_events.csv", usecols=["rep_key", "cargo"]).set_index("rep_key")
    B = {}

    def add(name, v):
        B[name] = B.get(name, np.zeros(N)) + v

    n_ev = {}
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64); dl = al - rl
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        assert vp.shape[0] == N
        rec2e = pd.Series(EM[f"{C}_bubble"], index=EM[f"{C}_rec"])
        for kind, m in (("ins", dl > SV_BP), ("del", -dl > SV_BP)):
            cols = np.where(m)[0]
            e = rec2e.reindex(cols).to_numpy()
            assert not np.isnan(e).any(), f"{C} {kind}: SV records outside the event map"
            ue, li = np.unique(e.astype(np.int64), return_inverse=True); li = li.ravel()
            M = sp.csr_matrix((np.ones(cols.size), (li, np.arange(cols.size))), shape=(ue.size, cols.size))
            Pf = ((M @ (vp[:, cols] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64)   # events x founders
            Qf = ((M @ (vc[:, cols] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64)
            c = np.asarray(Pf.sum(1)).ravel(); k = np.asarray(Qf.sum(1)).ravel()
            ok = (k >= CALL_MIN * N) & (c >= MIN_CAR)
            first = np.full(ue.size, -1, np.int64)
            for j in range(cols.size - 1, -1, -1):
                first[li[j]] = j
            keys = [f"{C}|{r}" for r in cols[first]]
            cat = CAT.reindex(keys)
            # every event with >= 1 carrier is in the catalog, keyed by its first record
            assert cat.loc[c > 0].kind.notna().all(), f"{C} {kind}: carried events missing from the catalog"
            assert (cat.loc[c > 0, "passes"].astype(bool).to_numpy() == ok[c > 0]).all(), \
                f"{C} {kind}: pass filter disagrees with the catalog's `passes`"
            cargo = EVC.cargo.reindex(keys).to_numpy()
            assert pd.notna(cargo[ok]).all(), f"{C} {kind}: passing events without a repeat class"
            cont = np.where(np.isin(cargo, ["te_full", "te_most"]), "TE", "nonTE")
            loc = np.select([cat.L_CDS.fillna(False).to_numpy(bool), cat.L_UTR.fillna(False).to_numpy(bool),
                             cat.L_intron.fillna(False).to_numpy(bool), cat.L_promoter.fillna(False).to_numpy(bool)],
                            ["CDS", "UTR", "intron", "promoter"], "intergenic")
            nc = np.asarray(Pf[:, GC].sum(1)).ravel(); nw = np.asarray(Pf[:, GW].sum(1)).ravel()
            OM = {"all": np.ones(ue.size, bool), "shared": (nc > 0) & (nw > 0)}
            Pf = Pf.tocsr(); Qf = Qf.tocsr()
            for o, om in OM.items():
                for ct in CONT:
                    for lc in LOCS:
                        s = np.where(ok & om & (cont == ct) & (loc == lc))[0]
                        pc = np.asarray(Pf[s].sum(0)).ravel(); qc = np.asarray(Qf[s].sum(0)).ravel()
                        for name in (f"{ct}_{lc}", f"{kind}_{ct}_{lc}"):          # pooled, and per direction
                            add(f"{name}_{o}", pc); add(f"called_{name}_{o}", qc)
                            n_ev[(name, o)] = n_ev.get((name, o), 0) + s.size
        print(f"[{C}] done", flush=True)

    D = pd.DataFrame(B); D.insert(0, "founder", FO.founder.to_numpy())
    D.to_csv(f"{OUT}/founder_burden_cargo.csv", index=False)
    print(f"\nwrote {OUT}/founder_burden_cargo.csv ({D.shape[0]} founders, {D.shape[1] - 1} columns)\n")
    print("PASSING EVENTS PER CLASS (shared / all origins)")
    print(f"  {'class':22s} {'pooled':>14s} {'insertions':>14s} {'deletions':>14s}")
    for ct in CONT:
        for lc in LOCS:
            f = lambda nm: f"{n_ev[(nm, 'shared')]:6,d} /{n_ev[(nm, 'all')]:6,d}"
            print(f"  {ct + '_' + lc:22s} {f(f'{ct}_{lc}'):>14s} {f(f'ins_{ct}_{lc}'):>14s} {f(f'del_{ct}_{lc}'):>14s}")


if __name__ == "__main__":
    main()

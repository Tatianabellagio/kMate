#!/usr/bin/env python
"""Per-founder SV load under the two possible units, split into common and rare.

The panel is the biallelic (HPRC/PanGenie) representation of a multi-allelic graph VCF: one event
appears as one record per distinct inserted sequence, up to ~78 of them. So "how many insertions
does this founder carry" has two answers:

  ALLELE unit  -- count records. Frequency = the record's own carrier count.
  EVENT unit   -- count Truvari events: alleles of one bubble collapsed by sequence/size similarity
                  (`_build_event_map.py`, `arch3_event_map.npz`, truvari-collapse defaults).
                  Frequency = carriers of the event (union).

and each is split by frequency class, because the two units can only disagree where an event has
several alleles (the rare end):

  common: >= 12 carriers (and <= N-12)     rare: 2-11 carriers      (singletons dropped)

Per founder this script writes the number carried and the kb carried, for insertions and deletions,
for every unit x class combination.

Env: kmate. Writes results/sv_adaptive/sv_unit_compare_event.npz.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SV_BP, MIN_MAC = 50, 12
UNITS, CLASSES = ("allele", "event"), ("common", "rare")


def main():
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    n = len(F)
    BM = np.load(f"{OUT}/arch3_event_map.npz")
    acc = {f"{m}_{t}_{u}_{c}": np.zeros(n) for m in ("n", "kb") for t in ("ins", "del") for u in UNITS for c in CLASSES}
    nunit = {f"{t}_{u}_{c}": 0 for t in ("ins", "del") for u in UNITS for c in CLASSES}
    multi = {t: [0, 0] for t in ("ins", "del")}          # [founder-bubble pairs with >1 carried record, total pairs]

    for ch in lib.CHROMS:
        cl = ch.lower(); base = f"{lib.PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        assert np.array_equal(meta["founders"].astype("U6"), F.founder.to_numpy()), f"{ch}: founder order"
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
        a_all = np.asarray(vp.sum(0)).ravel()
        bmap = pd.Series(BM[f"{ch}_bubble"], index=BM[f"{ch}_rec"])
        for t, sel in (("ins", (al - rl) > SV_BP), ("del", (rl - al) > SV_BP)):
            cols = np.where(sel & (a_all >= 1))[0]
            size = np.abs(al - rl)[cols].astype(float); a = a_all[cols]
            V = vp[:, cols]                                                   # [231 x ncols]
            ub, li = np.unique(bmap.reindex(cols).to_numpy().astype(np.int64), return_inverse=True); li = li.ravel()
            M = sp.csr_matrix((np.ones(cols.size), (li, np.arange(cols.size))), shape=(ub.size, cols.size))
            Cb = np.asarray((M @ V.T).todense())                              # [nbub x 231] records carried
            P = Cb > 0                                                        # bubble carried by founder
            nb = P.sum(1)                                                     # carriers per bubble
            multi[t][0] += int((Cb > 1).sum()); multi[t][1] += int(P.sum())
            nb_rec = nb[li]                                                   # bubble carriers, per record
            for u, freq_unit, freq_rec in (("allele", a, a), ("event", nb, nb_rec)):
                for c in CLASSES:
                    keep_u = (freq_unit >= MIN_MAC) & (freq_unit <= n - MIN_MAC) if c == "common" else (freq_unit >= 2) & (freq_unit < MIN_MAC)
                    keep_r = (freq_rec >= MIN_MAC) & (freq_rec <= n - MIN_MAC) if c == "common" else (freq_rec >= 2) & (freq_rec < MIN_MAC)
                    acc[f"kb_{t}_{u}_{c}"] += V @ np.where(keep_r, size, 0) / 1e3
                    acc[f"n_{t}_{u}_{c}"] += (V @ keep_r.astype(float)) if u == "allele" else P[keep_u].sum(0)
                    nunit[f"{t}_{u}_{c}"] += int(keep_u.sum())
        print(f"[{ch}] done", flush=True)

    for t in ("ins", "del"):
        print(f"{t}: founder-event pairs where the founder carries >1 record of the same event: "
              f"{multi[t][0]:,} of {multi[t][1]:,} ({100*multi[t][0]/multi[t][1]:.2f}%)")
    for k, v in nunit.items():
        print(f"  units {k}: {v:,}")
    np.savez_compressed(f"{OUT}/sv_unit_compare_event.npz", founders=F.founder.to_numpy().astype(str),
                        bio1=F.bio1.to_numpy().astype(float), group=F.group.to_numpy().astype(str),
                        is_cactus=F.is_cactus.to_numpy().astype(bool),
                        **acc, **{f"nunit_{k}": v for k, v in nunit.items()})
    print(f"[wrote] {OUT}/sv_unit_compare_event.npz")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Record-level SV table carrying BOTH context rules, for the context-cascade decomposition.

The burden model and the GEA report different context compositions for the same panel. Three things
differ between the two tracks, and this table holds all three on one set of rows so they can be varied
one at a time:

  * unit    -- arch3 writes one SV as several allele records; `event` is the truvari collapse
               (`arch3_event_map.npz`), so record -> event is a groupby on this table.
  * filter  -- `n_carriers` / `n_called` over the 231 founders (the burden model's >=2 / >=90% floor)
               in the panel table, and `MAF` (the GEA's pool frequency -- a different denominator, not
               the founder frequency) in the GEA table.
  * rule    -- `context` is the burden rule (`_sv_sfs_context.context`: CDS > UTR > intron > TE >
               promoter <=2kb > intergenic); `tier_1kb` / `tier_2kb` are the GEA rule
               (`screen_sig_blocks.classify` via `annotate_gea_universe.py`: gene body > promoter > TE),
               carried by `gea_hit_landscape.csv.gz`.

SV records are `|alt_len - ref_len| > 50`, matching PANEL_STATS and `_sv_sfs_context.SV_BP`, so the row
count reconciles with the panel-stats variant-class figure (226,468).

NOT JOINED. There is no usable record key between the two tables: (chrom, pos, ref_len, alt_len) is
not unique in arch3 -- the GEA's 27,586 SV rows collapse to 24,057 such tuples, because two distinct
records can share a position AND both lengths and differ only in the inserted/deleted SEQUENCE (e.g.
Chr1:352683 carries three separate 301 bp insertions, MAF 0.054 / 0.128 / 0.061). That is the
chrom:pos multiallelic trap one level deeper, and any merge on this tuple fans rows out.

Both context rules are pure functions of (chrom, pos, ref_len), so no join is needed: the burden rule is
recomputed directly on the GEA's own rows. Two tables are written -- `..._panel.csv.gz` (every panel SV
record, burden rule, founder carriers/called, truvari event) and `..._gea.csv.gz` (the GEA's tested
records carrying BOTH rules). Within a duplicate tuple every row gets the same context under either
rule, which is correct: they sit at the same place.

Env: kmate. Writes results/sv_adaptive/sv_context_decomposition_{panel,gea}.csv.gz.
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
GEA_LANDSCAPE = f"{lib.GEA}/genes/convergence/results/gea_hit_landscape.csv.gz"
SV_BP = 50
KEY = ["chrom", "pos", "ref_len", "alt_len"]


def main():
    IV = sc.layers()
    EM = np.load(f"{OUT}/arch3_event_map.npz")
    rows = []
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64)
        rl = meta["ref_len"].astype(np.int64)
        al = meta["alt_len"].astype(np.int64)
        dl = al - rl
        sv = np.where(np.abs(dl) > SV_BP)[0]

        # burden-rule context, vectorized over the SV records of this chromosome
        ctx = sc.context(IV, np.full(sv.size, C, dtype=object), pos[sv],
                         pos[sv] + np.maximum(rl[sv] - 1, 0))

        # per-record carriers / called over the 231 founders
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()[:, sv]
        vc = sp.load_npz(f"{base}.var_called.npz").tocsc()[:, sv]
        n_car = np.asarray((vp > 0).sum(0)).ravel().astype(np.int32)
        n_cal = np.asarray((vc > 0).sum(0)).ravel().astype(np.int32)

        # truvari event id (the unit the burden model reports); every SV record is in the map
        rec, ev = EM[f"{C}_rec"], EM[f"{C}_bubble"]
        assert np.array_equal(np.sort(rec), sv), f"{C}: event map does not cover the SV records"
        e = pd.Series(ev, index=rec).reindex(sv).to_numpy()

        rows.append(pd.DataFrame(dict(
            chrom=C, pos=pos[sv], ref_len=rl[sv], alt_len=al[sv],
            kind=np.where(dl[sv] > 0, "ins", "del"), size_bp=np.abs(dl[sv]),
            context=np.array(sc.CTX)[ctx], event=[f"{C}|{x}" for x in e],
            n_carriers=n_car, n_called=n_cal)))
        print(f"[{C}] {sv.size:,} SV records, {np.unique(e).size:,} events", flush=True)

    D = pd.concat(rows, ignore_index=True)
    D.to_csv(f"{OUT}/sv_context_decomposition_panel.csv.gz", index=False)
    print(f"\n[panel] {len(D):,} SV records | {D.event.nunique():,} truvari events")
    print(f"[wrote] {OUT}/sv_context_decomposition_panel.csv.gz")

    # the GEA's tested SV records, with the BURDEN rule recomputed on the same rows (no join:
    # see the module docstring -- the 4-tuple is not a record key)
    G = pd.read_csv(GEA_LANDSCAPE, usecols=KEY + ["cls", "kind", "MAF", "tier_1kb", "tier_2kb"])
    G = G[G.cls == "sv"].drop(columns="cls").reset_index(drop=True)
    gctx = np.empty(len(G), dtype=object)
    for C, idx in G.groupby("chrom").groups.items():
        i = np.asarray(idx)
        p = G.pos.to_numpy(np.int64)[i]; r = G.ref_len.to_numpy(np.int64)[i]
        gctx[i] = np.array(sc.CTX)[sc.context(IV, np.full(i.size, C, dtype=object),
                                              p, p + np.maximum(r - 1, 0))]
    G["context"] = gctx
    G.to_csv(f"{OUT}/sv_context_decomposition_gea.csv.gz", index=False)
    print(f"\n[gea] {len(G):,} tested SV records ({G[KEY].drop_duplicates().shape[0]:,} distinct "
          f"{'/'.join(KEY)} tuples -- not a record key)")
    print(f"[wrote] {OUT}/sv_context_decomposition_gea.csv.gz")


if __name__ == "__main__":
    main()

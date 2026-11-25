#!/usr/bin/env python
"""Per-founder SV burdens split into strata -- ALT-FREQUENCY bins (STRATUM=freq) or SIZE bins
(STRATUM=size). Insertions and deletions are emitted separately and pooled downstream by the model
script, matching the main model's parameterisation.

Why a separate script. `_sv_burden_build.py` hard-codes MIN_CAR = 2, so private variants never enter its
table, and it emits only two frequency classes (`all` >= 2 carriers, `rare` 2-10). This computes the same
event-level burdens with the carrier floor lowered to 1 and an arbitrary set of disjoint frequency bins. It
writes its OWN table so that `founder_burden_event.csv` -- which every fitted model to date depends on --
is never touched.

FREQUENCY bins are on the FULL 231-founder panel (population frequency is a property of the variant, so it stays the
same whichever founders a model is later fitted on). Default edges:

    private  MAC 1        the variant is unique to one founder
    rare     MAC 2-12     <= 5.2% of the panel
    common   MAC 13-231   > 5.2%

*** READ THIS BEFORE USING THE `private` COLUMNS ***
A private variant can only be detected in a founder that has its own assembly, so the private burden is
overwhelmingly a statement about sequencing technology, not biology: cactus founders carry a median 52
private events against 3 for short-read (17x, p = 2e-33), 98 of 151 short-read founders carry ZERO private
CDS events, and assembly route alone explains 37% of the variance in the private CDS burden rate. The
private bin is therefore usable ONLY within the 80-founder cactus panel, where detection is homogeneous
(private burden vs that founder's SV call rate: r = -0.17, p = 0.15). Do not fit it on all 231 founders.

SIZE bins (STRATUM=size) use the event's span in bp, taken from the event's first record -- the same
record whose genomic context is used to label the event, so size and context are assigned consistently.
Default edges 50-99 / 100-499 / 500+, which is the most balanced 3-bin split available in CDS
(444 / 481 / 482 events). Size is NOT independent of insertion/deletion labelling: the largest CDS bin is
80% insertions and the smallest 38%, because CDS insertions are ~3x larger than CDS deletions in this
panel. Since ins/del is defined relative to Col-0 and is therefore arbitrary, the burdens are pooled --
but the bin composition is printed so the entanglement stays visible.

Also emitted in both modes: `bp_<kind>_<ctx>` and `called_bp_<kind>_<ctx>`, the total base pairs carried
and callable per context, so a dose-style bp burden can be fitted without another pass over the genotypes.

Unit is the truvari-collapsed SV EVENT, matching the reported models. Call-rate floor 90%, as everywhere.
Env: kmate (needs lib + the arch3 panel matrices).
Writes results/sv_adaptive/founder_burden_event_{freq,size}.csv. Runtime ~8 min.
Env vars: STRATUM (freq|size), FREQ_EDGES, SIZE_EDGES, FREQ_NAMES, SIZE_NAMES, SV_CALL_MIN (0.9).
"""
import os, sys, importlib.util
import numpy as np, pandas as pd, scipy.sparse as sp

_H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_H))
import lib
_s = importlib.util.spec_from_file_location("sfsctx", os.path.join(_H, "_sv_sfs_context.py"))
sc = importlib.util.module_from_spec(_s); _s.loader.exec_module(sc)

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SV_BP = 50
CALL_MIN = float(os.environ.get("SV_CALL_MIN", 0.9))
STRATUM = os.environ.get("STRATUM", "freq")
assert STRATUM in ("freq", "size", "callrate", "tierdens", "tiersize")
# Carrier floor. STRATUM=freq defaults to 1 so the `private` bin exists at all; STRATUM=size defaults to 2,
# matching the main model's `all` class -- otherwise the size burdens silently absorb the singleton class,
# which is an assembly-detection artifact (cactus founders carry 17x more private events) and would
# contaminate every size bin with it.
MIN_CAR = int(os.environ.get("MIN_CAR", 1 if STRATUM == "freq" else 2))
# "lo:hi,lo:hi,..." with a name per bin, in order; bounds are inclusive
if STRATUM == "freq":
    EDGES = os.environ.get("FREQ_EDGES", "1:1,2:12,13:231")
    NAMES = os.environ.get("FREQ_NAMES", "private,rare,common").split(",")
elif STRATUM == "size":
    EDGES = os.environ.get("SIZE_EDGES", "50:99,100:499,500:1000000000")
    NAMES = os.environ.get("SIZE_NAMES", "small,medium,large").split(",")
elif STRATUM == "tiersize":
    # 2x2: call-rate tier x event SIZE. Low-call coding events are ~5x larger than high-call ones, so
    # size is the other candidate axis. If size carried the sign flip, hi_big would look like the
    # low-call events -- but at the 0.9 floor the LARGEST coding events are the most deleterious.
    #   code 1 hi_small   2 hi_big   3 lo_small   4 lo_big
    EDGES = os.environ.get("TS_EDGES", "1:1,2:2,3:3,4:4")
    NAMES = os.environ.get("TS_NAMES", "hi_small,hi_big,lo_small,lo_big").split(",")
    SIZE_CUT = int(os.environ.get("SIZE_CUT", 500))     # bp
    HICALL = int(os.environ.get("HICALL", 208))
elif STRATUM == "tierdens":
    # 2x2: call-rate tier x local chromatin context. The low-call events differ from the high-call ones
    # in two independent ways (bigger, and in gene-poorer regions), so this asks which axis carries the
    # sign flip: if it is location, lo_rich should behave like the high-call events, not like lo_poor.
    #   code 1 hi_rich   2 hi_poor   3 lo_rich   4 lo_poor
    EDGES = os.environ.get("TD_EDGES", "1:1,2:2,3:3,4:4")
    NAMES = os.environ.get("TD_NAMES", "hi_rich,hi_poor,lo_rich,lo_poor").split(",")
    DENS_MIN = float(os.environ.get("DENS_MIN", 212))   # genes/Mb; bottom quartile of coding SV events
    HICALL = int(os.environ.get("HICALL", 208))         # >= 208/231 is the reported 0.9 set
else:
    # CALL-RATE tiers, cut on the number of founders the event is genotyped in. Run this with
    # SV_CALL_MIN=0.5 so the relaxed events exist at all: the `high` bin is then exactly the set the
    # reported 0.9 model uses, and `low` is what relaxing the floor newly admits. Putting both in one
    # model asks whether the reported effect survives when the newly admitted events are held beside it.
    EDGES = os.environ.get("CALL_EDGES", "116:207,208:231")
    NAMES = os.environ.get("CALL_NAMES", "lowcall,highcall").split(",")
BINS = [(int(a), int(b), n) for (a, b), n in
        zip((e.split(":") for e in EDGES.split(",")), NAMES)]


def main():
    FO = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    N = len(FO)
    IV = sc.layers()
    BM = np.load(f"{OUT}/arch3_event_map.npz")
    CTX = sc.CTX
    B = {}
    def add(k, v): B[k] = B.get(k, np.zeros(N)) + v

    # sanity: the bins must be disjoint and cover a contiguous range, or "common = all - rare" style
    # arithmetic downstream would silently double-count
    for (a1, b1, n1), (a2, b2, n2) in zip(BINS, BINS[1:]):
        assert b1 < a2, f"bins {n1} and {n2} overlap"
    unit = {"freq": "MAC", "size": "bp", "callrate": "founders called",
            "tierdens": "code", "tiersize": "code"}[STRATUM]
    print(f"STRATUM={STRATUM}  MIN_CAR={MIN_CAR}  bins: "
          + ", ".join(f"{n} {unit} {a}-{b}" for a, b, n in BINS))
    print(f"call-rate floor: {CALL_MIN:.0%} of {N} founders\n")

    if STRATUM == "tierdens":
        GN = lib.load_genes()
        GST = {c: np.sort(g.start.to_numpy()) for c, g in GN.groupby("chrom")}
        WIN = 500_000          # +-500 kb -> genes per Mb
    n_ev = {n: 0 for _, _, n in BINS}
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64)
        al = meta["alt_len"].astype(np.int64); dl = al - rl
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
        vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        assert vp.shape[0] == N
        ctx = sc.context(IV, np.full(pos.size, C, dtype=object), pos, pos + np.maximum(rl - 1, 0))
        rec2b = pd.Series(BM[f"{C}_bubble"], index=BM[f"{C}_rec"])   # event map, same keys as the bubble map
        for kind, m in (("ins", dl > SV_BP), ("del", -dl > SV_BP)):
            cols = np.where(m)[0]
            b = rec2b.reindex(cols).to_numpy().astype(np.int64)
            ub, li = np.unique(b, return_inverse=True); li = li.ravel()
            M = sp.csr_matrix((np.ones(cols.size), (li, np.arange(cols.size))), shape=(ub.size, cols.size))
            # event carried / called = union over the event's alleles, exactly as _sv_burden_build.py
            Pf = ((M @ (vp[:, cols] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64)
            Qf = ((M @ (vc[:, cols] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64)
            c = np.asarray(Pf.sum(1)).ravel(); k = np.asarray(Qf.sum(1)).ravel()
            okc = (k >= CALL_MIN * N) & (c >= MIN_CAR)
            first = np.full(ub.size, -1, np.int64)
            for j in range(cols.size - 1, -1, -1):
                first[li[j]] = j
            bctx = ctx[cols[first]]
            bsize = np.abs(dl[cols[first]])              # event span, from the same record that gives context
            Pf = Pf.tocsr(); Qf = Qf.tocsr()
            # total bp carried / callable per context -- the dose-style burden, free in this pass
            for i, x in enumerate(CTX):
                s_ = np.where(okc & (bctx == i))[0]
                if s_.size:
                    w = bsize[s_].astype(np.float64)
                    add(f"bp_{kind}_{x}", np.asarray(Pf[s_].T.dot(w)).ravel())
                    add(f"called_bp_{kind}_{x}", np.asarray(Qf[s_].T.dot(w)).ravel())
            if STRATUM == "tierdens":
                epos = pos[cols[first]]
                st = GST[C]
                dens = (np.searchsorted(st, epos + WIN) - np.searchsorted(st, epos - WIN)) / (2 * WIN / 1e6)
                BV = np.where(k >= HICALL, 1, 3) + np.where(dens < DENS_MIN, 1, 0)
            elif STRATUM == "tiersize":
                BV = np.where(k >= HICALL, 1, 3) + np.where(bsize >= SIZE_CUT, 1, 0)
            else:
                BV = {"freq": c, "size": bsize, "callrate": k}[STRATUM]   # the variable the bins are cut on
            for lo, hi, nm in BINS:
                fm = okc & (BV >= lo) & (BV <= hi)
                n_ev[nm] += int(fm.sum())
                for i, x in enumerate(CTX):
                    s = np.where(fm & (bctx == i))[0]
                    add(f"{kind}_{x}_{nm}", np.asarray(Pf[s].sum(0)).ravel())
                    add(f"called_{kind}_{x}_{nm}", np.asarray(Qf[s].sum(0)).ravel())
        print(f"[{C}] done", flush=True)

    D = pd.DataFrame(B)
    D.insert(0, "founder", FO.founder.to_numpy())
    D["is_cactus"] = FO.is_cactus.to_numpy()
    D["bio1"] = FO.bio1.to_numpy()
    D["group"] = FO.group.to_numpy()
    # founder genotyping quality, carried over so a model can control for it (the rare-bin burden tracks it)
    prev = pd.read_csv(f"{OUT}/founder_burden_event.csv", dtype={"founder": str})
    D = D.merge(prev[["founder", "callrate_sv", "callrate_snp"]], on="founder", how="left")
    assert len(D) == N and D.founder.is_unique
    f = f"{OUT}/founder_burden_event_{STRATUM}.csv"
    D.to_csv(f, index=False)
    print(f"\n[wrote] {f}   {len(D)} founders x {D.shape[1]} columns")
    print("\nevents per bin (all contexts, ins+del):")
    for _, _, nm in BINS:
        print(f"  {nm:<9} {n_ev[nm]:>7,}")
    if STRATUM == "size":
        print("\nins/del composition per bin (size is partly a kind proxy -- see the header):")
        for _, _, nm in BINS:
            ni = D[[c for c in D.columns if c.startswith("ins_") and c.endswith(f"_{nm}")]].sum().sum()
            nd = D[[c for c in D.columns if c.startswith("del_") and c.endswith(f"_{nm}")]].sum().sum()
            print(f"  {nm:<9} carried ins {ni:>9,.0f} / del {nd:>9,.0f}   ins share {ni/(ni+nd):.0%}")
    print("\ncarried events per founder, by bin and technology:")
    cac = D.is_cactus.astype(str) == "True"
    for _, _, nm in BINS:
        v = D[[c for c in D.columns if c.startswith(("ins_", "del_")) and c.endswith(f"_{nm}")]].sum(axis=1)
        print(f"  {nm:<9} all {v.median():>7.0f}   cactus {v[cac].median():>7.0f}   "
              f"short-read {v[~cac].median():>7.0f}")


if __name__ == "__main__":
    main()

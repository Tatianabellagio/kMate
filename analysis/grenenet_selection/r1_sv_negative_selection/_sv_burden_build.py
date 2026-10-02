#!/usr/bin/env python
"""Per-founder variant burdens + SNP kinship, for the burden LMM (`_sv_burden_lmm.py`).

Question: do founders carrying more insertions of a class under purifying selection in nature
(shared CDS/UTR insertions, see sv_purging_fixed.ipynb section 7b) lose frequency in GrENE-Net,
beyond relatedness? Per-variant selection is not identifiable (pool AF is a projection onto
founder genotypes); founder fitness is.

BURDEN_f = number of variants of a class the founder carries (ALT present), classes =
  kind   : snp / indel (1-50 bp) / ins (insertion BUBBLES, > 50 bp; carrier = any allele) / del (bubbles)
  origin : all / shared / cold-specific / warm-specific   (same labels as _sv_purging_classify.py)
  context: CDS / UTR / intron / TE / promoter / intergenic (as _sv_sfs_context.py)
  freq   : all (>= 2 carriers) / rare (2-10 carriers)
Filter: called in >= 90% of the 231 founders, >= 2 carriers (private / singleton variants excluded).
Rare 2-3-carrier variants are still in -- the end the trajectory analysis could not reach (p0 >= 0.02).
Also `called_<class>` = number of that class's variants the founder is genotyped at, so burden can be
read as a rate (carried / called) -- missing calls fall unevenly across contexts (intergenic worst).
Also per founder: fraction of filtered SV / SNP records it is called at (genotyping-quality covariate).

KINSHIP: standardized GRM from a random 60k SNPs with >= 5 carriers, called >= 90%, missing = mean.
SV_UNIT=allele: insertions / deletions per panel RECORD (allele), filtered on the allele's own carriers
(>= 2) and call rate -- a bubble made only of private alleles is excluded (10.2% of passing insertion
bubbles are). SV_UNIT=event: alleles of a bubble collapsed with Truvari's matching rules
(`_build_event_map.py`), carriers = union over the event's alleles. Default SV_UNIT=bubble.
Env: kmate. Writes results/sv_adaptive/founder_burden{,_allele}.csv, founder_kinship.npy.
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
SV_BP, MIN_CAR, RARE, NK, SEED = 50, 2, 10, 60000, 0
CALL_MIN = float(os.environ.get("SV_CALL_MIN", 0.9))   # founder call-rate floor; != 0.9 -> outputs tagged _cr{pct}
CRTAG = "" if CALL_MIN == 0.9 else f"_cr{int(round(CALL_MIN * 100))}"
# a precedence-flipped build must not overwrite the reported burdens
if os.environ.get("CTX_PRECEDENCE", "te_first") == "promoter_first": CRTAG += "_promfirst"
# CTX_MODE: "cascade" (default) assigns each event ONE context by first match, so the burdens are
# disjoint. "overcount" drops the precedence entirely -- an event contributes to EVERY layer it
# overlaps, so a deletion hitting CDS and UTR counts once in each. The burdens then overlap, which
# raises their correlation, but no arbitrary ordering decides anything and each burden answers
# "how many SVs touch this class" rather than "touch it and nothing else".
# `intergenic` stays well-defined under both: events overlapping NO annotation layer.
CTX_MODE = os.environ.get("CTX_MODE", "cascade")
assert CTX_MODE in ("cascade", "overcount")
if CTX_MODE == "overcount": CRTAG += "_overcount"
CTX = sc.CTX
UNIT = os.environ.get("SV_UNIT", "bubble")          # bubble | allele | event
TAG = {"bubble": "", "allele": "_allele", "event": "_event"}[UNIT]
EVMAP = os.environ.get("EVENT_MAP", "truvari")    # truvari | strict (sensitivity: near-identical alleles only)
MAPF = (f"arch3_event_map{'_strict' if EVMAP == 'strict' else ''}.npz" if UNIT == "event" else "arch3_bubble_map.npz")
if UNIT == "event" and EVMAP == "strict":
    TAG = "_eventstrict"   # event = Truvari-collapsed alleles
TAG += CRTAG
ORIG = ["all", "shared", "cold", "warm"]


def origin_masks(nc, nw):
    return {"all": np.ones(nc.size, bool), "shared": (nc > 0) & (nw > 0),
            "cold": (nc >= 2) & (nw == 0), "warm": (nw >= 2) & (nc == 0)}


def main():
    rng = np.random.default_rng(SEED)
    FO = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    GC = (FO.group == "cold").to_numpy(); GW = (FO.group == "warm").to_numpy()
    N = len(FO)
    IV = sc.layers()
    BM = np.load(f"{OUT}/{MAPF}")
    B = {}                                     # column name -> per-founder counts
    def add(name, v):
        B[name] = B.get(name, np.zeros(N)) + v
    Kg = []                                    # kinship SNP genotype blocks
    call_sv = np.zeros(N); tot_sv = 0; call_snp = np.zeros(N); tot_snp = 0
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64)
        al = meta["alt_len"].astype(np.int64); dl = al - rl
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        assert vp.shape[0] == N
        _ch = np.full(pos.size, C, dtype=object); _hi = pos + np.maximum(rl - 1, 0)
        ctx = sc.context(IV, _ch, pos, _hi)
        if CTX_MODE == "overcount":
            _H2 = {k: ann.hit(v, _ch, pos, _hi) for k, v in IV.items()}
            _L = {"CDS": _H2["CDS"], "UTR": _H2["UTR"], "intron": _H2["gene"] & ~_H2["exon"],
                  "TE": _H2["TE"], "promoter": _H2["promoter"]}
            _L["intergenic"] = ~(_L["CDS"] | _L["UTR"] | _L["intron"] | _L["TE"] | _L["promoter"])
            LMASK = np.vstack([_L[x] for x in CTX])          # len(CTX) x n_records
        # record-level: SNPs and indels
        for kind, m in (("snp", (rl == 1) & (al == 1)),
                        ("indel", (np.abs(dl) >= 1) & (np.abs(dl) <= SV_BP) & ~((rl == 1) & (al == 1)))):
            cols = np.where(m)[0]
            P = (vp[:, cols] > 0).tocsc(); Q = (vc[:, cols] > 0).tocsc()
            c = np.asarray(P.sum(0)).ravel(); k = np.asarray(Q.sum(0)).ravel()
            ok = (k >= CALL_MIN * N) & (c >= MIN_CAR)
            nc = np.asarray(P[GC].sum(0)).ravel(); nw = np.asarray(P[GW].sum(0)).ravel()
            OM = origin_masks(nc, nw)
            for o in ORIG:
                for fq, fm in (("all", ok), ("rare", ok & (c <= RARE))):
                    for i, x in enumerate(CTX):
                        sel = LMASK[i][cols] if CTX_MODE == "overcount" else (ctx[cols] == i)
                        s = np.where(fm & OM[o] & sel)[0]
                        add(f"{kind}_{o}_{x}_{fq}", np.asarray(P[:, s].sum(1)).ravel())
                        add(f"called_{kind}_{o}_{x}_{fq}", np.asarray(Q[:, s].sum(1)).ravel())
            if kind == "snp":
                call_snp += np.asarray(Q[:, np.where(ok)[0]].sum(1)).ravel(); tot_snp += ok.sum()
                kk = np.where((k >= 0.9 * N) & (c >= 5) & (c <= N - 5))[0]      # kinship SNPs: always the 0.9 floor
                kk = kk[rng.random(kk.size) < 0.02]
                G = P[:, kk].toarray().astype(np.float32); Cm = Q[:, kk].toarray()
                G[~Cm] = np.nan
                Kg.append(G)
        if UNIT == "allele":
            # SVs per record (allele): own carriers >= 2, own call rate >= 90%
            for kind, m in (("ins", dl > SV_BP), ("del", -dl > SV_BP)):
                cols = np.where(m)[0]
                P = (vp[:, cols] > 0).tocsc(); Q = (vc[:, cols] > 0).tocsc()
                c = np.asarray(P.sum(0)).ravel(); k = np.asarray(Q.sum(0)).ravel()
                ok = (k >= CALL_MIN * N) & (c >= MIN_CAR)
                OM = origin_masks(np.asarray(P[GC].sum(0)).ravel(), np.asarray(P[GW].sum(0)).ravel())
                for o in ORIG:
                    for fq, fm in (("all", ok), ("rare", ok & (c <= RARE))):
                        for i, x in enumerate(CTX):
                            sel = LMASK[i][cols] if CTX_MODE == "overcount" else (ctx[cols] == i)
                            s = np.where(fm & OM[o] & sel)[0]
                            add(f"{kind}_{o}_{x}_{fq}", np.asarray(P[:, s].sum(1)).ravel())
                            add(f"called_{kind}_{o}_{x}_{fq}", np.asarray(Q[:, s].sum(1)).ravel())
                call_sv += np.asarray(Q[:, np.where(ok)[0]].sum(1)).ravel(); tot_sv += ok.sum()
            print(f"[{C}] done", flush=True)
            continue
        # SVs: bubbles
        rec2b = pd.Series(BM[f"{C}_bubble"], index=BM[f"{C}_rec"])
        for kind, m in (("ins", dl > SV_BP), ("del", -dl > SV_BP)):
            cols = np.where(m)[0]
            b = rec2b.reindex(cols).to_numpy().astype(np.int64)
            ub, li = np.unique(b, return_inverse=True); li = li.ravel()
            M = sp.csr_matrix((np.ones(cols.size), (li, np.arange(cols.size))), shape=(ub.size, cols.size))
            Pf = ((M @ (vp[:, cols] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64)   # bubbles x founders
            Qf = ((M @ (vc[:, cols] > 0).astype(np.float64).T.tocsr()) > 0).astype(np.float64)
            c = np.asarray(Pf.sum(1)).ravel(); k = np.asarray(Qf.sum(1)).ravel()
            ok = (k >= CALL_MIN * N) & (c >= MIN_CAR)
            first = np.full(ub.size, -1, np.int64)
            for j in range(cols.size - 1, -1, -1):
                first[li[j]] = j
            bctx = ctx[cols[first]]
            bL = LMASK[:, cols[first]] if CTX_MODE == "overcount" else None
            nc = np.asarray(Pf[:, GC].sum(1)).ravel(); nw = np.asarray(Pf[:, GW].sum(1)).ravel()
            OM = origin_masks(nc, nw)
            Pf = Pf.tocsr(); Qf = Qf.tocsr()
            for o in ORIG:
                for fq, fm in (("all", ok), ("rare", ok & (c <= RARE))):
                    for i, x in enumerate(CTX):
                        sel = bL[i] if CTX_MODE == "overcount" else (bctx == i)
                        s = np.where(fm & OM[o] & sel)[0]
                        add(f"{kind}_{o}_{x}_{fq}", np.asarray(Pf[s].sum(0)).ravel())
                        add(f"called_{kind}_{o}_{x}_{fq}", np.asarray(Qf[s].sum(0)).ravel())
            call_sv += np.asarray(Qf[np.where(ok)[0]].sum(0)).ravel(); tot_sv += ok.sum()
        print(f"[{C}] done", flush=True)

    D = pd.DataFrame(B).copy(); D.insert(0, "founder", FO.founder.to_numpy())
    D["callrate_sv"] = call_sv / tot_sv; D["callrate_snp"] = call_snp / tot_snp
    D["is_cactus"] = FO.is_cactus.to_numpy(); D["bio1"] = FO.bio1.to_numpy(); D["group"] = FO.group.to_numpy()
    D.to_csv(f"{OUT}/founder_burden{TAG}.csv", index=False)
    G = np.hstack(Kg)
    mu = np.nanmean(G, 0); sd = np.nanstd(G, 0)
    G = np.where(np.isnan(G), mu, G); Z = (G - mu) / sd
    K = Z @ Z.T / Z.shape[1]
    if not CRTAG:                                # kinship is unit- and floor-independent: written by the default run only
        np.save(f"{OUT}/founder_kinship.npy", K)
    print(f"kinship from {Z.shape[1]:,} SNPs; burden columns {len(B)}")
    print(f"[wrote] {OUT}/founder_burden{TAG}.csv")


if __name__ == "__main__":
    main()

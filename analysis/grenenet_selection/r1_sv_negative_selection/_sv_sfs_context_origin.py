#!/usr/bin/env python
"""Site-frequency spectra split by GENOMIC CONTEXT: CDS / UTR / intron / TE / promoter / intergenic.

The classic purifying-selection test: a class of variants under selection has a spectrum shifted
towards rare variants relative to a neutral-ish class IN THE SAME CONTEXT. Here the question is
whether insertions show a larger rare-variant excess over SNPs in coding sequence than elsewhere.

Built on the full 231-founder panel:
  * context: same rules as `_sv_purging_annotate.py`, first match wins
             CDS > UTR > intron > TE > promoter (<= 2 kb upstream, strand-aware) > intergenic
  * classes: SNP and short indel (1-50 bp) -- random 5% of records as references; every insertion
             and deletion record (> 50 bp), as ALLELES and as BUBBLES (alleles of one graph bubble
             merged, carriers = union; `_build_bubble_map.py`).
  * filter:  called in >= 90% of ALL 231 founders (panel-level, not per group) and >= 1 carrier.
             Alleles are filtered on their OWN call rate and carriers, frequency = own carriers / own
             called, context = the record's own (bubbles: first record's). Allele is the primary unit.
             Singletons kept: this is a spectrum, and the rare end is the signal.
  * bubble context = context of the bubble's first record.
Origin version: also keeps each variant's carriers among cold- and warm-tertile founders
(`sv_origin_partition_founders.csv`), so spectra can be split cold-specific (>= 2 cold, 0 warm) /
warm-specific / shared -- the same labels as `_sv_purging_classify.py`. Bubbles: union carriers.
SV_UNIT=event: the `*_bubble` classes are built on Truvari-collapsed EVENTS instead (`_build_event_map.py`).
Env: kmate. Writes results/sv_adaptive/sv_sfs_context_origin{,_event}.npz.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import importlib.util
_s = importlib.util.spec_from_file_location("ann", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                                "_sv_purging_annotate.py"))
ann = importlib.util.module_from_spec(_s); _s.loader.exec_module(ann)

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SV_BP, KEEP, SEED = 50, 0.05, 0
CALL_MIN = float(os.environ.get("SV_CALL_MIN", 0.9))   # founder call-rate floor; != 0.9 -> outputs tagged _cr{pct}
CRTAG = "" if CALL_MIN == 0.9 else f"_cr{int(round(CALL_MIN * 100))}"
CTX = ["CDS", "UTR", "intron", "TE", "promoter", "intergenic"]
UNIT = os.environ.get("SV_UNIT", "bubble")          # the "_bubble" classes use this map: bubble | event
EVMAP = os.environ.get("EVENT_MAP", "truvari")    # truvari | strict
MAPF = (f"arch3_event_map{'_strict' if EVMAP == 'strict' else ''}.npz" if UNIT == "event" else "arch3_bubble_map.npz")
TAG = (("_eventstrict" if EVMAP == "strict" else "_event") if UNIT == "event" else "") + CRTAG
CLASSES = ["snp", "indel", "ins_allele", "ins_bubble", "del_allele", "del_bubble"]


def layers():
    g = pd.read_csv(lib.TAIR10_GENES_TE, sep="\t", header=None, comment="#",
                    names=["chrom", "src", "feat", "start", "end", "score", "strand", "frame", "attr"])
    g = g[g.chrom.isin(lib.CHROMS)]
    genes = g[g.feat == "gene"]
    prom = genes.copy()
    prom["start"] = np.where(genes.strand == "+", genes.start - 2000, genes.end + 1)
    prom["end"] = np.where(genes.strand == "+", genes.start - 1, genes.end + 2000)
    return {"CDS": ann.merged(g[g.feat == "CDS"]),
            "UTR": ann.merged(g[g.feat.isin(["five_prime_UTR", "three_prime_UTR"])]),
            "exon": ann.merged(g[g.feat == "exon"]), "gene": ann.merged(genes),
            "TE": ann.merged(g[g.feat == "transposable_element"]), "promoter": ann.merged(prom)}


def context(IV, ch, lo, hi):
    H = {k: ann.hit(v, ch, lo, hi) for k, v in IV.items()}
    cat = np.full(lo.size, 5, np.int8)
    cat[H["promoter"]] = 4; cat[H["TE"]] = 3
    cat[H["gene"] & ~H["exon"]] = 2; cat[H["UTR"]] = 1; cat[H["CDS"]] = 0
    return cat


def main():
    rng = np.random.default_rng(SEED)
    IV = layers()
    BM = np.load(f"{OUT}/{MAPF}")
    FO = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    GC = (FO.group == "cold").to_numpy(); GW = (FO.group == "warm").to_numpy()
    F = {f"{c}_{v}": [] for c in CLASSES for v in ("f", "n", "nc", "nw", "ctx")}
    def put(cname, f, n, nc, nw, cx):
        for v, a in (("f", f.astype(np.float32)), ("n", n.astype(np.int16)), ("nc", nc.astype(np.int16)), ("nw", nw.astype(np.int16)), ("ctx", cx.astype(np.int8))):
            F[f"{cname}_{v}"].append(a)
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64)
        al = meta["alt_len"].astype(np.int64); dl = al - rl
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        N = vp.shape[0]
        chv = np.full(pos.size, C, dtype=object)
        ctx = context(IV, chv, pos, pos + np.maximum(rl - 1, 0))
        sub = rng.random(pos.size) < KEEP
        # SNPs and indels: record-level
        for cname, m in (("snp", (rl == 1) & (al == 1) & sub),
                         ("indel", (np.abs(dl) >= 1) & (np.abs(dl) <= SV_BP) & ~((rl == 1) & (al == 1)) & sub)):
            cols = np.where(m)[0]
            c = np.asarray((vp[:, cols] > 0).sum(0)).ravel(); k = np.asarray((vc[:, cols] > 0).sum(0)).ravel()
            ok = (k >= CALL_MIN * N) & (c > 0)
            Pb = (vp[:, cols] > 0)
            nc = np.asarray(Pb[GC].sum(0)).ravel(); nw = np.asarray(Pb[GW].sum(0)).ravel()
            put(cname, c[ok] / k[ok], c[ok], nc[ok], nw[ok], ctx[cols][ok])
        # SVs: alleles and bubbles, paired filter
        rec2b = pd.Series(BM[f"{C}_bubble"], index=BM[f"{C}_rec"])
        for sv, m in (("ins", dl > SV_BP), ("del", -dl > SV_BP)):
            cols = np.where(m)[0]
            b = rec2b.reindex(cols).to_numpy().astype(np.int64)
            ub, li = np.unique(b, return_inverse=True); li = li.ravel()
            M = sp.csr_matrix((np.ones(cols.size), (li, np.arange(cols.size))), shape=(ub.size, cols.size))
            P = (vp[:, cols] > 0).astype(np.float64).T.tocsr(); Cc = (vc[:, cols] > 0).astype(np.float64).T.tocsr()
            cb = np.asarray(((M @ P) > 0).sum(1)).ravel(); mb = np.asarray(((M @ Cc) > 0).sum(1)).ravel()
            okb = (mb >= CALL_MIN * N) & (cb > 0)
            first = np.full(ub.size, -1, np.int64)
            for j in range(cols.size - 1, -1, -1):
                first[li[j]] = j
            bctx = ctx[cols[first]]
            ca = np.asarray(P.sum(1)).ravel(); ka = np.asarray(Cc.sum(1)).ravel()
            oka = (ka >= CALL_MIN * N) & (ca > 0)          # allele: its OWN call rate and carriers
            Pc = P[:, GC]; Pw = P[:, GW]
            bc = np.asarray(((M @ Pc) > 0).sum(1)).ravel(); bw = np.asarray(((M @ Pw) > 0).sum(1)).ravel()
            ac = np.asarray(Pc.sum(1)).ravel(); aw = np.asarray(Pw.sum(1)).ravel()
            put(f"{sv}_bubble", cb[okb] / mb[okb], cb[okb], bc[okb], bw[okb], bctx[okb])
            put(f"{sv}_allele", ca[oka] / ka[oka], ca[oka], ac[oka], aw[oka], ctx[cols][oka])
        print(f"[{C}] done", flush=True)
    S = {k: np.concatenate(v) for k, v in F.items()}
    np.savez_compressed(f"{OUT}/sv_sfs_context_origin{TAG}.npz", **S)
    for c in CLASSES:
        nc, nw = S[f"{c}_nc"], S[f"{c}_nw"]
        print(f"{c:<12} n={nc.size:>9,}  cold-spec={((nc>=2)&(nw==0)).sum():>8,}  warm-spec={((nw>=2)&(nc==0)).sum():>8,}  shared={((nc>0)&(nw>0)).sum():>8,}")
    print(f"[wrote] {OUT}/sv_sfs_context_origin{TAG}.npz")


if __name__ == "__main__":
    main()

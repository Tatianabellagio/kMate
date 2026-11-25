#!/usr/bin/env python
"""Site-frequency spectra split by GENOMIC CONTEXT: CDS / UTR / intron / TE / promoter / intergenic.

The classic purifying-selection test: a class of variants under selection has a spectrum shifted
towards rare variants relative to a neutral-ish class IN THE SAME CONTEXT. Here the question is
whether insertions show a larger rare-variant excess over SNPs in coding sequence than elsewhere.

Built on the full 231-founder panel:
  * context: same rules as `_sv_purging_annotate.py`, first match wins
             CDS > UTR > intron > TE > promoter (<= 1 kb upstream, strand-aware) > intergenic
  * classes: SNP and short indel (1-50 bp) -- random 5% of records as references; every insertion
             and deletion record (> 50 bp), as ALLELES and as BUBBLES (alleles of one graph bubble
             merged, carriers = union; `_build_bubble_map.py`). Bubbles are the unit to read --
             allele-splitting inflates the rare end, unevenly across contexts.
  * filter:  called in >= 90% of ALL 231 founders (panel-level, not per group) and >= 1 carrier.
             Alleles inherit their bubble's filter and use its called set as denominator.
             Singletons kept: this is a spectrum, and the rare end is the signal.
  * bubble context = context of the bubble's first record.
Env: kmate. Writes results/sv_adaptive/sv_sfs_context.npz.
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
SV_BP, KEEP, CALL_MIN, SEED = 50, 0.05, 0.9, 0
# PROJECT-WIDE promoter window, upstream of the TSS. 1 kb, not 2 kb: the median Arabidopsis gene sits
# only 863 bp from its nearest upstream neighbour, so a 2 kb window lies inside another gene's body for
# 73.4% of genes against 54.2% at 1 kb. 1 kb is also what TAIR distributes as its standard upstream set.
# The GEA classifier (genes/dissection/screen_sig_blocks.py, PROM_BP) uses the same value.
PROM_BP = 1000
PRECEDENCE = os.environ.get("CTX_PRECEDENCE", "te_first")
assert PRECEDENCE in ("te_first", "promoter_first")
CTX = ["CDS", "UTR", "intron", "TE", "promoter", "intergenic"]
CLASSES = ["snp", "indel", "ins_allele", "ins_bubble", "del_allele", "del_bubble"]


def layers():
    g = pd.read_csv(lib.TAIR10_GENES_TE, sep="\t", header=None, comment="#",
                    names=["chrom", "src", "feat", "start", "end", "score", "strand", "frame", "attr"])
    g = g[g.chrom.isin(lib.CHROMS)]
    genes = g[g.feat == "gene"]
    prom = genes.copy()
    prom["start"] = np.where(genes.strand == "+", genes.start - PROM_BP, genes.end + 1)
    prom["end"] = np.where(genes.strand == "+", genes.start - 1, genes.end + PROM_BP)
    return {"CDS": ann.merged(g[g.feat == "CDS"]),
            "UTR": ann.merged(g[g.feat.isin(["five_prime_UTR", "three_prime_UTR"])]),
            "exon": ann.merged(g[g.feat == "exon"]), "gene": ann.merged(genes),
            "TE": ann.merged(g[g.feat == "transposable_element"]), "promoter": ann.merged(prom)}


def context(IV, ch, lo, hi):
    H = {k: ann.hit(v, ch, lo, hi) for k, v in IV.items()}
    cat = np.full(lo.size, 5, np.int8)
    # Later assignment wins. CTX_PRECEDENCE flips the ONE ordering choice with no biological
    # justification: 26% of events called TE also sit in a promoter window, and nothing says TE should
    # outrank promoter. "te_first" (default, the project rule) assigns promoter then lets TE overwrite;
    # "promoter_first" does the reverse. Gene body still outranks both either way -- that ordering IS
    # justified, since disrupting coding sequence outranks sitting upstream of a neighbour.
    if PRECEDENCE == "promoter_first":
        cat[H["TE"]] = 3; cat[H["promoter"]] = 4
    else:
        cat[H["promoter"]] = 4; cat[H["TE"]] = 3
    cat[H["gene"] & ~H["exon"]] = 2; cat[H["UTR"]] = 1; cat[H["CDS"]] = 0
    return cat


def main():
    rng = np.random.default_rng(SEED)
    IV = layers()
    BM = np.load(f"{OUT}/arch3_bubble_map.npz")
    F = {f"{c}_{x}": [] for c in CLASSES for x in CTX}
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
            for i, x in enumerate(CTX):
                s = ok & (ctx[cols] == i)
                F[f"{cname}_{x}"].append((c[s] / k[s]).astype(np.float32))
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
            ca = np.asarray(P.sum(1)).ravel()
            oka = okb[li] & (ca > 0)
            for i, x in enumerate(CTX):
                sb = okb & (bctx == i)
                F[f"{sv}_bubble_{x}"].append((cb[sb] / mb[sb]).astype(np.float32))
                sa = oka & (bctx[li] == i)
                F[f"{sv}_allele_{x}"].append((ca[sa] / mb[li][sa]).astype(np.float32))
        print(f"[{C}] done", flush=True)
    S = {k: (np.concatenate(v) if v else np.zeros(0, np.float32)) for k, v in F.items()}
    np.savez_compressed(f"{OUT}/sv_sfs_context.npz", **S)
    print(f"\n{'class':<12}" + "".join(f"{x:>12}" for x in CTX))
    for c in CLASSES:
        print(f"{c:<12}" + "".join(f"{S[f'{c}_{x}'].size:>12,}" for x in CTX))
    print(f"\nshare <= 5%, and fold vs SNPs in the same context:")
    for c in CLASSES:
        row = f"{c:<12}"
        for x in CTX:
            f = S[f"{c}_{x}"]; r = (S[f"snp_{x}"] <= 0.05).mean()
            row += f"{100*(f<=0.05).mean():>6.0f}%{(f<=0.05).mean()/r:>5.2f}x"
        print(row)
    print(f"\n[wrote] {OUT}/sv_sfs_context.npz")


if __name__ == "__main__":
    main()

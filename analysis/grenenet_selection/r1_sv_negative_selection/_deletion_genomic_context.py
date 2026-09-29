#!/usr/bin/env python
"""What do SV deletions remove? Genomic context of the deleted interval, for every deletion.

Deletion counterpart of `_insertion_genomic_context.py`. A deletion has reference coordinates
[pos, pos+ref_len), so unlike an insertion it can be classified by what the removed interval
OVERLAPS (priority CDS > UTR > exon_noncoding > intron > intergenic) -- that is the biologically
relevant question ("does it remove coding sequence"). Because overlap probability grows with
deletion length (median 363 bp, p99 31 kb) while the genome-fraction baseline is point-based,
the MIDPOINT class is stored too (`genic_mid`) for a like-for-like comparison with insertions.

No MAC floor, same as the insertion table. Env: kmate.
Writes results/sv_adaptive/deletion_context.csv.
"""
import os, sys, importlib.util
import numpy as np, pandas as pd
import scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import lib
_sp = importlib.util.spec_from_file_location("insctx", os.path.join(HERE, "_insertion_genomic_context.py"))
insctx = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(insctx)

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
KEY = "/global/scratch/users/tbellg/gea_grene-net/key_files"
SV_BP = 50
CHROMS = [f"Chr{i}" for i in range(1, 6)]


def overlaps(a, b, iv, back=4000):
    """Boolean per query interval [a,b]: does it overlap any feature interval? iv=(starts,ends) sorted by start."""
    st, en = iv
    hit = np.zeros(a.size, bool)
    if st.size == 0:
        return hit
    emax = np.maximum.accumulate(en)
    idx = np.searchsorted(st, b, side="right") - 1          # last feature starting at or before b
    cand = (idx >= 0) & (emax[np.clip(idx, 0, None)] >= a)
    for i in np.where(cand)[0]:
        j = idx[i]; lo = max(0, j - back)
        hit[i] = np.any(en[lo:j + 1] >= a[i])
    return hit


def classify(in_cds, in_utr, in_exon, in_gene):
    cls = np.full(in_cds.size, "intergenic", dtype=object)
    cls[in_gene & ~in_exon] = "intron"
    cls[in_exon & ~in_cds & ~in_utr] = "exon_noncoding"
    cls[in_utr] = "UTR"
    cls[in_cds] = "CDS"
    return cls


def main():
    print("[1] loading TAIR10 features")
    feats = insctx.load_features()
    eco = pd.read_csv(f"{KEY}/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
    eco["ecotypeid"] = eco["ecotypeid"].astype(str); em = eco.set_index("ecotypeid")

    rows = []
    print("[2] scanning panel for SV deletions (no MAC floor)")
    for ch in CHROMS:
        cl = ch.lower()
        base = f"{lib.PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        founders = meta["founders"].astype("U6")
        pos = meta["pos"].astype(np.int64)
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        home1 = np.array([float(em.loc[f, "bio1"]) if f in em.index else np.nan for f in founders])
        home18 = np.array([float(em.loc[f, "bio18"]) if f in em.index else np.nan for f in founders])

        dele = (rl - al) > SV_BP
        n_alt = np.asarray(vp.sum(0)).ravel(); n_cal = np.asarray(vc.sum(0)).ravel()
        cols = np.where(dele & (n_alt >= 1) & (n_cal >= 1))[0]
        if cols.size == 0:
            continue
        V = vp[:, cols]
        na = np.asarray(V.sum(0)).ravel().astype(float); nc = np.asarray(vc[:, cols].sum(0)).ravel().astype(float)
        cb1 = np.asarray(np.nan_to_num(home1)[None, :] @ V).ravel() / np.maximum(na, 1)
        cb18 = np.asarray(np.nan_to_num(home18)[None, :] @ V).ravel() / np.maximum(na, 1)

        a = pos[cols] + 1; b = pos[cols] + rl[cols] - 1            # removed REF bases, excluding the shared anchor base (1-based, inclusive)
        mid = (a + b) // 2
        ov = {k: overlaps(a, b, feats[k][ch]) for k in ("CDS", "UTR", "exon", "gene", "TE")}
        pt = {k: insctx.covered(mid, feats[k][ch]) for k in ("CDS", "UTR", "exon", "gene", "TE")}
        rows.append(pd.DataFrame(dict(
            chrom=ch, rec=cols, pos=a, size=(rl - al)[cols], n_alt=na, n_called=nc,
            freq=na / np.maximum(nc, 1),
            genic=classify(ov["CDS"], ov["UTR"], ov["exon"], ov["gene"]), te_overlap=ov["TE"],
            genic_mid=classify(pt["CDS"], pt["UTR"], pt["exon"], pt["gene"]), te_mid=pt["TE"],
            carrier_bio1=cb1, carrier_bio18=cb18)))
        print(f"  {ch}: {cols.size:,} SV deletions  (overlap CDS {int(ov['CDS'].sum())}, "
              f"midpoint CDS {int(pt['CDS'].sum())})", flush=True)

    D = pd.concat(rows, ignore_index=True)
    D.to_csv(f"{OUT}/deletion_context.csv", index=False)
    print(f"\n[3] {len(D):,} SV deletions; median size {D['size'].median():,.0f} bp, p99 {D['size'].quantile(.99):,.0f}")
    print("overlap class:\n" + D.genic.value_counts().to_string())
    print("midpoint class:\n" + D.genic_mid.value_counts().to_string())
    print(f"\n[wrote] {OUT}/deletion_context.csv")


if __name__ == "__main__":
    main()

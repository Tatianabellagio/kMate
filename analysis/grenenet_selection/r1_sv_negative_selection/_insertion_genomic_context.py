#!/usr/bin/env python
"""Where do SV insertions land, and does the frequency spectrum differ by genomic context?

The question: of the inserted sequence that is enriched in cold-origin ecotypes and declines
in warm gardens, what part of the genome does it sit in -- coding, intronic, intergenic?

WHAT THIS CAN AND CANNOT ANSWER. An "insertion" here is sequence present in a founder and
ABSENT from TAIR10, so the inserted sequence itself has no reference coordinates. The GFF can
therefore only classify **where the insertion landed** (the insertion site), never what the
inserted sequence contains. Characterising the inserted sequence needs the ALT allele plus a
repeat/homology analysis and is deliberately out of scope for this pass.

NO ANCESTRAL POLARIZATION. This project has no outgroup (see SV_TEMPORAL_PURGING_SUMMARY.md),
so frequency here is the ALT (insertion-present) frequency among called founders, NOT derived
allele frequency. An insertion absent from Col-0 may be derived, or may be ancestral with a
Col-0 deletion. The absolute spectrum is therefore not a DAF spectrum. The CROSS-CLASS
comparison is still valid, because whatever polarization error exists is shared across
CDS / intron / intergenic.

NO MAC FLOOR. Every other section of this arm filters MAC>=12. That is wrong here: purifying
selection shows up in the RARE tail, which that filter removes. This script uses all
segregating SV records. Consequence: the per-variant climate slope beta, which only exists for
the MAC>=12 subset, is available for a minority of them -- handled as a separate split.

Classification is two independent axes, because a TE insertion inside an intron is both:
  genic class  CDS > UTR > exon(non-coding) > intron > intergenic   (priority order)
  te_overlap   does the site fall in transposable_element / transposon_fragment /
               transposable_element_gene

Env: kmate. Reads panel/arch3 + the TAIR10 genes_transposons GFF + founder home climate.
Writes results/sv_adaptive/insertion_context.{npz,csv}.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PROJ = lib.PROJ
OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
KEY = "/global/scratch/users/tbellg/gea_grene-net/key_files"
SV_BP = 50
CHROMS = [f"Chr{i}" for i in range(1, 6)]
GENIC = ["CDS", "UTR", "exon_noncoding", "intron", "intergenic"]


def load_features():
    """Per chrom, sorted (start,end) arrays for each feature group."""
    g = pd.read_csv(lib.TAIR10_GENES_TE, sep="\t", header=None, comment="#",
                    names=["chrom", "src", "feat", "start", "end", "score",
                           "strand", "frame", "attr"])
    g = g[g.chrom.isin(CHROMS)]
    groups = {
        "CDS": ["CDS"],
        "UTR": ["five_prime_UTR", "three_prime_UTR"],
        "exon": ["exon", "pseudogenic_exon"],
        "gene": ["gene", "pseudogene"],
        "TE": ["transposable_element", "transposon_fragment", "transposable_element_gene"],
    }
    out = {}
    for name, feats in groups.items():
        sub = g[g.feat.isin(feats)]
        out[name] = {}
        for ch in CHROMS:
            s = sub[sub.chrom == ch]
            a = np.sort(np.column_stack([s.start.to_numpy(), s.end.to_numpy()]), axis=0) \
                if len(s) else np.zeros((0, 2), int)
            # sort by start, keep pairs intact
            st = s.start.to_numpy(); en = s.end.to_numpy()
            o = np.argsort(st)
            out[name][ch] = (st[o], en[o])
        n = sum(len(out[name][c][0]) for c in CHROMS)
        print(f"  [gff] {name:<6} {n:,} intervals")
    return out


def covered(pos, iv):
    """Boolean: does each pos fall inside any interval? iv=(starts,ends) sorted by start."""
    st, en = iv
    if st.size == 0:
        return np.zeros(pos.size, bool)
    # running max of ends lets a single searchsorted answer containment
    emax = np.maximum.accumulate(en)
    idx = np.searchsorted(st, pos, side="right") - 1
    ok = idx >= 0
    hit = np.zeros(pos.size, bool)
    # candidate check against running max, then exact scan on the ambiguous ones
    cand = ok & (emax[np.clip(idx, 0, None)] >= pos)
    for i in np.where(cand)[0]:
        j = idx[i]
        lo = max(0, j - 200)
        hit[i] = np.any((st[lo:j + 1] <= pos[i]) & (en[lo:j + 1] >= pos[i]))
    return hit


def genome_fraction(feats):
    """bp fraction of the genome in each class, for the enrichment null."""
    tot = 0; frac = {}
    lens = {}
    for ch in CHROMS:
        st, en = feats["gene"][ch]
        lens[ch] = max(en.max() if en.size else 0, 1)
    tot = sum(lens.values())
    for name in ("CDS", "UTR", "exon", "gene", "TE"):
        bp = 0
        for ch in CHROMS:
            st, en = feats[name][ch]
            if st.size == 0:
                continue
            o = np.argsort(st); s2, e2 = st[o], en[o]
            # merge overlapping intervals
            m_s, m_e = [s2[0]], [e2[0]]
            for a, b in zip(s2[1:], e2[1:]):
                if a <= m_e[-1] + 1:
                    m_e[-1] = max(m_e[-1], b)
                else:
                    m_s.append(a); m_e.append(b)
            bp += int(np.sum(np.array(m_e) - np.array(m_s) + 1))
        frac[name] = bp / tot
    return frac, tot


def main():
    print("[1] loading TAIR10 features")
    feats = load_features()
    frac, glen = genome_fraction(feats)
    print(f"  genome span {glen/1e6:.1f} Mb | genome fraction: "
          + ", ".join(f"{k} {v:.3f}" for k, v in frac.items()))

    # founder home climate, aligned to panel order
    eco = pd.read_csv(f"{KEY}/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
    eco["ecotypeid"] = eco["ecotypeid"].astype(str)
    em = eco.set_index("ecotypeid")

    rows = []
    print("[2] scanning panel for SV insertions (no MAC floor)")
    for ch in CHROMS:
        cl = ch.lower()
        base = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        founders = meta["founders"].astype("U6")
        pos = meta["pos"].astype(np.int64)
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
        vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        home1 = np.array([float(em.loc[f, "bio1"]) if f in em.index else np.nan
                          for f in founders])
        home18 = np.array([float(em.loc[f, "bio18"]) if f in em.index else np.nan
                           for f in founders])

        ins = (al - rl) > SV_BP                       # insertions only
        n_alt = np.asarray(vp.sum(0)).ravel()
        n_cal = np.asarray(vc.sum(0)).ravel()
        keep = ins & (n_alt >= 1) & (n_cal >= 1)
        cols = np.where(keep)[0]
        if cols.size == 0:
            continue
        V = vp[:, cols]
        na = np.asarray(V.sum(0)).ravel().astype(float)
        nc = np.asarray(vc[:, cols].sum(0)).ravel().astype(float)
        # carrier-mean origin climate (no experimental data involved)
        cb1 = np.asarray(np.nan_to_num(home1)[None, :] @ V).ravel() / np.maximum(na, 1)
        cb18 = np.asarray(np.nan_to_num(home18)[None, :] @ V).ravel() / np.maximum(na, 1)

        p = pos[cols]
        in_cds = covered(p, feats["CDS"][ch])
        in_utr = covered(p, feats["UTR"][ch])
        in_exon = covered(p, feats["exon"][ch])
        in_gene = covered(p, feats["gene"][ch])
        in_te = covered(p, feats["TE"][ch])
        cls = np.full(p.size, "intergenic", dtype=object)
        cls[in_gene & ~in_exon] = "intron"
        cls[in_exon & ~in_cds & ~in_utr] = "exon_noncoding"
        cls[in_utr] = "UTR"
        cls[in_cds] = "CDS"

        rows.append(pd.DataFrame(dict(
            chrom=ch, rec=cols, pos=p, size=(al - rl)[cols], n_alt=na, n_called=nc,
            freq=na / np.maximum(nc, 1), genic=cls, te_overlap=in_te,
            carrier_bio1=cb1, carrier_bio18=cb18)))
        print(f"  {ch}: {cols.size:,} SV insertions  "
              f"(CDS {int(in_cds.sum())}, UTR {int(in_utr.sum())}, "
              f"intron {int((in_gene & ~in_exon).sum())}, TE {int(in_te.sum())})", flush=True)

    D = pd.concat(rows, ignore_index=True)
    D.to_csv(f"{OUT}/insertion_context.csv", index=False)
    np.savez_compressed(f"{OUT}/insertion_context.npz",
                        genome_frac=np.array([frac.get(k, np.nan)
                                              for k in ("CDS", "UTR", "exon", "gene", "TE")]),
                        genome_frac_keys=np.array(["CDS", "UTR", "exon", "gene", "TE"]))
    print(f"\n[3] {len(D):,} SV insertions total")
    print(D.genic.value_counts().to_string())
    print(f"  TE-overlapping: {int(D.te_overlap.sum()):,} ({100*D.te_overlap.mean():.1f}%)")
    print(f"  frequency: median {D.freq.median():.4f}, "
          f"{100*(D.freq <= 1/231).mean():.1f}% are singletons")
    print(f"  size: median {D['size'].median():,.0f} bp, p99 {D['size'].quantile(0.99):,.0f} bp")
    print(f"\n[wrote] {OUT}/insertion_context.csv + .npz")


if __name__ == "__main__":
    main()

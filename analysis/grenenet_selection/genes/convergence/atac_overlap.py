#!/usr/bin/env python
"""Do our candidate non-SNP variants sit in open chromatin more than chance?

The ATAC track in the MOI-LAB zoom-Manhattan figures is decoration: a static multi-tissue
peak union drawn as pale spans behind the association panel. The useful version is a
number -- does a candidate indel/SV overlap an ATAC peak, in how many tissues, and is that
more than a matched non-candidate indel/SV would?

Peaks: weiwei/synthetic_evolution/ATAC-seq/ATAC-seq_multitissue.bed -- the source BED the
Drive's ATAC-seq_multitissue.csv was derived from. Replicates intersected within each of
flower / leaf / root / shoot, then the four tissues unioned. 99,332 disjoint peaks,
17.19 Mb = 14.4% of the nuclear genome; `num` = in how many tissues the peak was called.

Candidates: the pooled Bonferroni non-SNP hits (results/variants_classified.csv, sv +
smallindel; MNPs excluded as substitutions -- they change no length and the ATAC question
is about the same variant class the rest of this layer uses).

Background: results/repeat_variants_background_v2.csv.gz, built by build_repeat_variants.py
-- every panel indel/SV record with seed-mix p0 MAF > 0.05, i.e. approximately the set the
GEA could have found. Reused rather than rebuilt so "testable" means the same thing here as
it does there.

Three things the test has to survive, in increasing order of how much they hurt:

1. COMPOSITION. Candidates are 15% SV against 5% in the background, and are shifted in MAF.
   Both correlate with overlap. Matched on vclass x MAF bin x |size| bin.
2. LD. The 2,469 candidates are not 2,469 independent draws -- they pile into haplotype
   blocks, so a per-variant Fisher test is anti-conservative in exactly the way this
   project keeps getting burned by. The unit of the headline test is the LOCUS (200 kb
   collapse, as in build_convergence.add_loci), not the variant.
3. CIRCULARITY. ATAC peaks concentrate at promoters, and the candidate set is promoter-rich
   (743/2,469). "Candidates are in open chromatin" could just be "candidates are in
   promoters". Reported stratified by region as well, which is the non-circular version:
   among promoter variants only, are candidates still more often in a peak?

Env: kmate. Compute node, ~1 min. Writes results/atac_overlap{,_candidates}.csv.
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, GEA); sys.path.insert(0, HERE)

OUT = f"{HERE}/results"
ATAC = ("/global/scratch/projects/fc_moilab/weiwei/synthetic_evolution/"
        "ATAC-seq/ATAC-seq_multitissue.bed")
TISSUES = ["flower", "leaf", "root", "shoot"]
LOCUS_KB = 200
N_PERM = 10000
SEED = 0


# ------------------------------------------------------------------------------------------------
def load_peaks() -> dict:
    """chrom -> (start_1b, end_1b, num, tissue_matrix). BED is 0-based half-open."""
    A = pd.read_csv(ATAC, sep="\t")
    A = A.rename(columns={f"{t}/{t}_intersect.bed": t for t in TISSUES})
    A["chrom"] = "Chr" + A.chrom.astype(str)
    out = {}
    for ch, g in A.groupby("chrom"):
        g = g.sort_values("start")
        assert (g.start.values[1:] >= g.end.values[:-1]).all(), f"{ch}: peaks overlap"
        out[ch] = (g.start.values + 1, g.end.values, g.num.values,
                   g[TISSUES].values.astype(np.int8))
    return out


def peak_overlap(peaks: dict, chrom, pos, ref_len) -> pd.DataFrame:
    """n_tissues (0 = no peak) and the per-tissue flags for each variant.

    Variant span is the REF footprint [pos, pos + ref_len - 1] -- a point for an insertion,
    the deleted bases for a deletion. A variant spanning several peaks takes the max `num`
    and the union of tissues.
    """
    chrom = np.asarray(chrom); pos = np.asarray(pos, np.int64)
    ref_len = np.asarray(ref_len, np.int64)
    lo, hi = pos, pos + np.maximum(ref_len, 1) - 1
    n_tis = np.zeros(len(pos), np.int8)
    n_pk = np.zeros(len(pos), np.int32)
    tis = np.zeros((len(pos), len(TISSUES)), np.int8)
    for ch, (ps, pe, num, tm) in peaks.items():
        m = np.where(chrom == ch)[0]
        if not len(m):
            continue
        i0 = np.searchsorted(pe, lo[m], side="left")       # first peak ending at/after lo
        i1 = np.searchsorted(ps, hi[m], side="right")      # first peak starting after hi
        k = i1 - i0
        n_pk[m] = np.maximum(k, 0)
        one = np.where(k == 1)[0]                          # the overwhelming majority
        n_tis[m[one]] = num[i0[one]]
        tis[m[one]] = tm[i0[one]]
        for j in np.where(k > 1)[0]:                       # a long SV across several peaks
            sl = slice(i0[j], i1[j])
            n_tis[m[j]] = num[sl].max()
            tis[m[j]] = tm[sl].max(axis=0)
    return pd.DataFrame({"atac_n_tissues": n_tis, "atac_n_peaks": n_pk,
                         "in_atac": n_tis > 0,
                         **{f"atac_{t}": tis[:, i] for i, t in enumerate(TISSUES)}})


# ------------------------------------------------------------------------------------------------
def add_locus(D: pd.DataFrame, kb=LOCUS_KB) -> pd.Series:
    """Collapse variants within `kb` on a chromosome into one locus -- the LD-ish unit."""
    lab = np.empty(len(D), object)
    n = 0
    for ch, g in D.groupby("chrom"):
        o = g.pos.values.argsort()
        idx = g.index.values[o]
        gap = np.diff(g.pos.values[o]) > kb * 1000
        grp = np.concatenate([[0], np.cumsum(gap)])
        lab[D.index.get_indexer(idx)] = [f"{ch}:{n + k}" for k in grp]
        n += grp.max() + 1
    return pd.Series(lab, index=D.index)


def strata(vclass, maf, size) -> np.ndarray:
    """vclass x MAF bin x |size| bin -- the three covariates that predict peak overlap."""
    mafb = np.digitize(np.asarray(maf, float), [0.075, 0.10, 0.15, 0.20, 0.30])
    szb = np.digitize(np.asarray(size, float), [2, 5, 10, 25, 50, 200, 1000])
    return np.array([f"{c}|{m}|{s}" for c, m, s in zip(vclass, mafb, szb)])


def main():
    rng = np.random.default_rng(SEED)
    peaks = load_peaks()
    tot = sum(int((pe - ps + 1).sum()) for ps, pe, _, _ in peaks.values())
    print(f"ATAC: {sum(len(p[0]) for p in peaks.values()):,} peaks, {tot/1e6:.2f} Mb "
          f"({100*tot/119.7e6:.1f}% of the nuclear genome)\n")

    # ---- candidates ----------------------------------------------------------------------
    V = pd.read_csv(f"{OUT}/variants_classified.csv")
    C = V[V.vclass.isin(["sv", "smallindel"])].copy().reset_index(drop=True)
    C["maf"] = C.gea_MAF.fillna(C.gwas_maf)
    C["abs_size"] = C["size"].abs()
    C = pd.concat([C, peak_overlap(peaks, C.chrom, C.pos, C.ref_len)], axis=1)
    C["locus"] = add_locus(C)
    print(f"candidates: {len(C):,} non-SNP hits "
          f"({(C.vclass=='sv').sum()} SV, {(C.vclass=='smallindel').sum()} small indel) "
          f"in {C.locus.nunique():,} loci")

    # ---- background ----------------------------------------------------------------------
    B = pd.read_csv(f"{OUT}/repeat_variants_background_v2.csv.gz",
                    usecols=["chrom", "pos", "ref_len", "alt_len", "maf", "size", "vclass",
                             "in_cds"])
    B = pd.concat([B, peak_overlap(peaks, B.chrom, B.pos, B.ref_len)], axis=1)
    B["abs_size"] = B["size"].abs()
    print(f"background: {len(B):,} testable indel/SV records "
          f"({(B.vclass=='sv').sum():,} SV)\n")

    # ---- raw rates -----------------------------------------------------------------------
    print("== raw overlap with an ATAC peak ==")
    print(f"{'set':<28}{'in peak':>16}{'rate':>8}{'mean tissues|in':>17}")
    for lab, D in [("candidates (all non-SNP)", C), ("  SV", C[C.vclass == "sv"]),
                   ("  small indel", C[C.vclass == "smallindel"]),
                   ("background (testable)", B), ("  SV", B[B.vclass == "sv"]),
                   ("  small indel", B[B.vclass == "smallindel"])]:
        mt = D.loc[D.in_atac, "atac_n_tissues"].mean() if D.in_atac.any() else float("nan")
        print(f"{lab:<28}{int(D.in_atac.sum()):>9,}/{len(D):<6,}{100*D.in_atac.mean():>7.1f}%"
              f"{mt:>17.2f}")

    # ---- matched, locus-level permutation ------------------------------------------------
    # Only candidates that exist in the background's stratum space can be matched; the rest
    # (mostly GWAS hits below the MAF>0.05 floor) are reported but not tested.
    C["stratum"] = strata(C.vclass, C.maf, C.abs_size)
    B["stratum"] = strata(B.vclass, B.maf, B.abs_size)
    pool = {s: g.in_atac.values for s, g in B.groupby("stratum")}
    ok = C.stratum.isin(pool).values & C.maf.notna().values
    T = C[ok]
    print(f"\n== matched permutation ==\n{len(T):,}/{len(C):,} candidates have a matched "
          f"background stratum ({(~ok).sum()} without: MAF below the GEA floor or an "
          f"unrepresented size)")

    # statistic: mean over LOCI of the within-locus overlap rate (each locus weighs 1)
    loci = T.groupby("locus")
    obs = loci.in_atac.mean().mean()
    sizes = loci.size()
    strat_by_locus = [g.stratum.values for _, g in loci]
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        null[i] = np.mean([np.mean([pool[s][rng.integers(len(pool[s]))] for s in ss])
                           for ss in strat_by_locus])
    p = (1 + (null >= obs).sum()) / (1 + N_PERM)
    print(f"locus-level overlap rate: observed {100*obs:.1f}%  "
          f"matched null {100*null.mean():.1f}% (sd {100*null.std():.1f})  "
          f"fold {obs/null.mean():.2f}x  p = {p:.4f}   "
          f"[{len(sizes):,} loci, median {int(sizes.median())} variants each]")

    # per-variant Fisher, for reference only -- anti-conservative under LD
    from scipy.stats import fisher_exact
    exp = np.mean([pool[s].mean() for s in T.stratum])
    print(f"per-variant (matched expectation): observed {100*T.in_atac.mean():.1f}%  "
          f"expected {100*exp:.1f}%  fold {T.in_atac.mean()/exp:.2f}x")
    print("   (no p-value quoted per variant -- LD makes it anti-conservative)")

    # ---- the non-circular version: stratified by region ----------------------------------
    print("\n== by region (candidates only; the background has no region annotation) ==")
    print(f"{'region':<18}{'n':>7}{'in peak':>10}{'rate':>8}{'mean tis|in':>13}")
    for reg, g in C.groupby("region"):
        if len(g) < 20:
            continue
        mt = g.loc[g.in_atac, "atac_n_tissues"].mean() if g.in_atac.any() else float("nan")
        print(f"{reg:<18}{len(g):>7,}{int(g.in_atac.sum()):>10,}"
              f"{100*g.in_atac.mean():>7.1f}%{mt:>13.2f}")

    C.drop(columns=["stratum"]).to_csv(f"{OUT}/atac_overlap.csv", index=False)
    print(f"\nwrote {OUT}/atac_overlap.csv")


if __name__ == "__main__":
    main()

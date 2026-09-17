#!/usr/bin/env python
"""Is the ATAC overlap of our candidates more than a matched non-candidate indel/SV? -- part 2.

atac_overlap.py established the raw rates: candidates 21.6% vs testable background 21.5%,
i.e. nothing, with one apparent exception -- a locus-weighted rate of 25.4% against a
matched null of 21.1% (1.20x, p=0.0099). This script exists to decide whether that
exception survives, because two things about it are suspect:

1. THE NULL'S VARIANCE IS TOO SMALL. Drawing k independent background variants for a
   k-variant locus gets the null MEAN right but not its SPREAD: real variants inside one
   locus are metres apart and share peak status (all inside one promoter, all inside one
   TE), so real per-locus rates are far more extreme than independent draws produce. An
   understated null sd makes the p anti-conservative -- the same failure mode as
   `persite-gwas-low-mac-tail-inflation` and `wza-block-calibration`. Null B below draws
   CONTIGUOUS runs of background variants instead, so the null carries the clustering too.

2. IT MAY BE REGION COMPOSITION, NOT ACCESSIBILITY. Peaks concentrate at promoters
   (candidates: promoter 34.5%, 5'UTR 54.5% in a peak, against TE 4.7%, intron 9.4%).
   "Candidates are in open chromatin" and "candidates are in promoters" are different
   claims and only the first is interesting. The background is therefore region-annotated
   with the SAME classifier the candidates went through
   (dissection/screen_sig_blocks.classify) and compared region by region.

Env: kmate. Compute node. The background classification is an iterrows loop over the GFF,
so it runs on a subsample (SUB records, still ~30x the candidate set) and is cached.
Writes results/atac_enrich_background.csv.gz and prints the two tests.
"""
from __future__ import annotations
import os, sys, importlib.util as ilu
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, GEA); sys.path.insert(0, HERE)
import atac_overlap as AO                                        # noqa: E402

OUT = f"{HERE}/results"
SUB = 80_000           # background records to region-annotate
N_PERM = 10000
SEED = 0


def screen_mod():
    p = os.path.join(os.path.dirname(HERE), "dissection", "screen_sig_blocks.py")
    spec = ilu.spec_from_file_location("screen_sig_blocks", p)
    mod = ilu.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def main():
    rng = np.random.default_rng(SEED)
    peaks = AO.load_peaks()

    C = pd.read_csv(f"{OUT}/atac_overlap.csv")
    C["abs_size"] = C["size"].abs()
    print(f"candidates: {len(C):,} in {C.locus.nunique()} loci\n")

    B = pd.read_csv(f"{OUT}/repeat_variants_background_v2.csv.gz",
                    usecols=["chrom", "pos", "ref_len", "alt_len", "maf", "size", "vclass"])
    B = pd.concat([B, AO.peak_overlap(peaks, B.chrom, B.pos, B.ref_len)], axis=1)
    B["abs_size"] = B["size"].abs()
    B = B.sort_values(["chrom", "pos"]).reset_index(drop=True)

    # ---- NULL B: contiguous runs of background variants (keeps the clustering) -----------
    print("== the locus-level result, under an LD-honest null ==")
    loci = C.groupby("locus")
    obs = loci.in_atac.mean().mean()
    ks = loci.size().values
    starts = {ch: (B.chrom.values == ch).argmax() for ch in B.chrom.unique()}
    counts = B.chrom.value_counts()
    chs = np.array(list(starts)); w = np.array([counts[c] for c in chs], float); w /= w.sum()
    hit = B.in_atac.values
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        tot = 0.0
        for k in ks:
            c = rng.choice(len(chs), p=w)
            lo = starts[chs[c]]
            j = lo + rng.integers(0, max(counts[chs[c]] - k, 1))
            tot += hit[j:j + k].mean()
        null[i] = tot / len(ks)
    p = (1 + (null >= obs).sum()) / (1 + N_PERM)
    print(f"  observed {100*obs:.1f}%   contiguous-block null {100*null.mean():.1f}% "
          f"(sd {100*null.std():.1f})   fold {obs/null.mean():.2f}x   p = {p:.4f}")
    print(f"  [null sd was {1.7:.1f} under independent draws -- the clustering roughly "
          f"{null.std()*100/1.7:.1f}x it]")

    # how concentrated is the observed locus statistic?
    lr = loci.in_atac.mean().sort_values(ascending=False)
    print(f"  loci fully in peaks: {(lr == 1).sum()}/{len(lr)}; fully outside: {(lr == 0).sum()}; "
          f"top-10 loci contribute {100*lr.head(10).sum()/lr.sum():.0f}% of the statistic")

    # ---- region-matched comparison -------------------------------------------------------
    cache = f"{OUT}/atac_enrich_background.csv.gz"
    if os.path.exists(cache):
        BR = pd.read_csv(cache)
        print(f"\nreusing region-annotated background ({len(BR):,})")
    else:
        idx = rng.choice(len(B), size=min(SUB, len(B)), replace=False)
        S = B.iloc[np.sort(idx)].reset_index(drop=True)
        print(f"\nregion-annotating {len(S):,} background records with the candidates' own "
              f"classifier (iterrows over the TAIR10 GFF, several minutes) ...", flush=True)
        sm = screen_mod()
        V, _ = sm.classify(S)
        BR = pd.concat([S, V[["tier", "region"]].reset_index(drop=True)], axis=1) \
            if "region" in V.columns else pd.concat([S, V.reset_index(drop=True)], axis=1)
        BR.to_csv(cache, index=False)
        print(f"  wrote {cache}")

    print("\n== ATAC overlap by region: candidates vs background in the SAME region ==")
    print(f"{'region':<18}{'cand n':>8}{'cand':>8}{'bg n':>9}{'bg':>8}{'fold':>8}{'p':>10}")
    from scipy.stats import fisher_exact
    for reg in ["promoter", "5'UTR", "3'UTR", "CDS", "intron", "TE", "intergenic"]:
        c = C[C.region == reg]; b = BR[BR.region == reg]
        if len(c) < 20 or len(b) < 20:
            continue
        tab = [[int(c.in_atac.sum()), len(c) - int(c.in_atac.sum())],
               [int(b.in_atac.sum()), len(b) - int(b.in_atac.sum())]]
        _, pv = fisher_exact(tab)
        fold = (c.in_atac.mean() / b.in_atac.mean()) if b.in_atac.mean() else float("nan")
        print(f"{reg:<18}{len(c):>8,}{100*c.in_atac.mean():>7.1f}%{len(b):>9,}"
              f"{100*b.in_atac.mean():>7.1f}%{fold:>8.2f}{pv:>10.3g}")
    print("  (Fisher per region, variant-level: still LD-optimistic, read the fold not the p)")

    print("\n== does region composition alone explain the candidate rate? ==")
    rate = BR.groupby("region").in_atac.mean()
    exp = C.region.map(rate).dropna()
    print(f"  candidates observed {100*C.in_atac.mean():.1f}%   "
          f"expected from the candidates' region mix {100*exp.mean():.1f}%   "
          f"fold {C.in_atac.mean()/exp.mean():.2f}x")


if __name__ == "__main__":
    main()

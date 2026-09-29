#!/usr/bin/env python
"""End-to-end audit of the GPX6 candidate -- every link in the chain, re-derived.

GPX6 (AT4G11600) is the strongest candidate in the set, which is exactly why it needs
this. Each check below either confirms a claim we are making or replaces it with a
weaker one. Run it whenever the upstream tables change.

The claim being audited: a 1,164 bp deletion 330 bp upstream of GPX6's TSS is favoured in
cold gardens and purged in warm ones, removes open chromatin and destroys TF binding
sites, and is therefore a causal candidate for climate adaptation.

  1  allele identity     the position is multiallelic (11 panel records, alt_len 50-63)
                         and two af_store rows share the (chrom,pos,ref_len,alt_len) key
  2  inflation           lambda reported as context; genomic control is deliberately NOT
                         applied (see the note in the check itself)
  3  frequency           p0 vs gen9 mean; is the reported MAF the right quantity
  4  gene attribution    is GPX6 really the nearest gene, and is it a promoter
  5  climate gradient    garden-level, which is the pseudo-replication-free version
  6  empirical calibration  where that gradient sits against matched background variants
  7  locus attribution   lead or passenger, against non-SNP AND SNP neighbours
  8  panel support       call rate, which everything downstream inherits

env: kmate. Compute node, ~3 min (reads the pool matrices).
"""
from __future__ import annotations
import os, sys, subprocess
import numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
import lib                                                       # noqa: E402
import axis_clusters as ac                                       # noqa: E402

OUT = f"{HERE}/results"
PROJ = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
BCF = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/bcftools"
ROW, CHROM, POS, REF_LEN, ALT_LEN = 1517046, "Chr4", 7011705, 1225, 61
BONF = {"sv": 5.41, "nonsnp": 7.02, "smallindel": 7.02}
WIN, NBG = 25_000, 4000


def hdr(n, t):
    print(f"\n{'='*78}\n{n}. {t}\n{'='*78}")


def garden_r(M, p0, cols, site, clim):
    """Per-garden mean change from p0, correlated with the climate axis.

    The GEA's unit is the pool (352) but pools are nested in 31 gardens, so a pool-level
    test is pseudo-replicated -- that is most of what drives lambda. Aggregating to the
    garden first is the same question asked honestly.
    """
    X = np.asarray(M[:, cols], dtype=np.float32) - p0[cols][None, :]
    S = pd.DataFrame(X); S["site"] = site
    G = S.groupby("site").mean()
    c = clim.reindex(G.index).values
    ok = np.isfinite(c); cc = c[ok] - c[ok].mean()
    Z = G.values[ok] - np.nanmean(G.values[ok], 0)
    return np.nansum(Z * cc[:, None], 0) / np.sqrt(np.nansum(Z ** 2, 0) * np.sum(cc ** 2))


def main():
    rng = np.random.default_rng(0)
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())
    pc1 = A["pc1"]

    # -- 1 allele identity ----------------------------------------------------------------
    hdr(1, "ALLELE IDENTITY at a multiallelic position")
    M9 = pd.read_csv(f"{OUT}/screen_record_gea_match.csv")
    m = M9[(M9.chrom == CHROM) & (M9.pos == POS)]
    print(m[["ref_len", "alt_len", "store_row", "key_ambiguous", "maf_g9",
             "gea_sig_record"]].to_string(index=False))
    print(f"\n  -> the MAF rule keeps store_row {ROW}; the collider (maf 0.034) is rejected.")
    print("     plot_variant_garden_grid refuses this key without an explicit store_row,")
    print("     and results/figures/INDEX.csv records 1517046 -- the figure used it.")

    # -- 2 inflation ----------------------------------------------------------------------
    hdr(2, "INFLATION: reported as context, NOT corrected for")
    P = pd.read_csv(f"{OUT}/gea_pool_records.csv")
    p = P[(P.chrom == CHROM) & (P.pos == POS)].copy()
    p["bonf"] = p.cls.map(BONF)
    print(p[["cls", "axis", "nlp", "lam", "bonf"]].round(3).to_string(index=False))
    print(f"\n  {int((p.nlp > p.bonf).sum())}/{len(p)} records clear their per-class "
          f"Bonferroni on RAW p. lambda spans {p.lam.min():.2f}-{p.lam.max():.2f}.")
    print("\n  DO NOT apply genomic control here. Standing project decision "
          "(memory `raw-lfmm-over-wza-decision`,")
    print("  README 'Why'): climate adaptation IS correlated with climate, so a large part "
          "of lambda is")
    print("  real polygenic signal rather than confounding. Dividing it out removes the "
          "signal along with")
    print("  the artifact and would empty the candidate list. These p-values are a "
          "CANDIDATE NET, not a")
    print("  calibrated test -- which is exactly why this layer ranks on independent lines "
          "of evidence")
    print("  (checks 5-8 below, plus ATAC and TFBS) instead of on p. The two class scans "
          "at this key are")
    print("  the same variant tested in the `sv` and `nonsnp` runs, not two alleles; "
          "n_axes deduplicates.")

    # -- 3 frequency ----------------------------------------------------------------------
    hdr(3, "FREQUENCY: p0 vs gen9")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy",
                mmap_mode="r")
    a = np.asarray(M[:, ROW]); mu = np.nanmean(a)
    print(f"  p0 = {p0[ROW]:.4f}   gen9 pool mean = {mu:.4f}   "
          f"MAF = {min(mu, 1-mu):.4f}   net {mu - p0[ROW]:+.4f}")
    print(f"  pools with data: {int(np.isfinite(a).sum())}/{len(a)} -- the POOL layer is "
          f"complete; the missingness in check 8 is in the founder VCF, a different thing.")

    # -- 4 gene attribution ---------------------------------------------------------------
    hdr(4, "GENE ATTRIBUTION")
    G = lib.load_genes(); g = G[G.gene == "AT4G11600"].iloc[0]
    vs, ve = POS, POS + REF_LEN - 1
    nb = G[(G.chrom == CHROM) & (G.end > vs - 15000) & (G.start < ve + 15000)].copy()
    nb["dist"] = np.where((nb.start <= ve) & (nb.end >= vs), 0,
                          np.minimum((nb.start - ve).abs(), (vs - nb.end).abs()))
    print(nb.sort_values("dist")[["gene", "start", "end", "strand", "dist"]]
          .head(4).to_string(index=False))
    print(f"\n  -> GPX6 is on the minus strand, so its TSS is its END ({g.end:,}); the "
          f"deletion starts {vs - g.end} bp upstream of it and does not touch the gene body.")

    # -- 5 climate gradient ---------------------------------------------------------------
    hdr(5, "CLIMATE GRADIENT, per garden")
    r_gpx = float(garden_r(M, p0, np.array([ROW]), site, pc1)[0])
    X = np.asarray(M[:, ROW]) - p0[ROW]
    D = pd.DataFrame({"site": site, "d": X}).groupby("site").d.mean()
    c = pc1.reindex(D.index)
    pr, pp = pearsonr(c, D); sr, sp = spearmanr(c, D)
    print(f"  gardens = {len(D)};  Pearson r = {pr:+.3f} (p={pp:.2e});  "
          f"Spearman rho = {sr:+.3f} (p={sp:.2e})")
    print(f"  gardens where the deletion rose: {int((D > 0).sum())}/{len(D)}")

    # -- 6 empirical calibration ----------------------------------------------------------
    hdr(6, "IS THAT GRADIENT UNUSUAL? empirical null from matched background")
    B = pd.read_csv(f"{OUT}/repeat_variants_background_v2.csv.gz",
                    usecols=["chrom", "pos", "ref_len", "alt_len"])
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    key = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                        "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    key["row"] = np.arange(len(key))
    rows = key.merge(B, on=["chrom", "pos", "ref_len", "alt_len"])["row"].values
    samp = np.unique(np.append(rng.choice(rows, NBG, replace=False), ROW))
    rb = garden_r(M, p0, np.sort(samp), site, pc1)
    pct = 100 * np.nanmean(np.abs(rb) < abs(r_gpx))
    emp = float(np.nanmean(np.abs(rb) >= abs(r_gpx)))
    print(f"  background: {len(rows):,} testable indel/SV records, {len(samp):,} sampled")
    print(f"  |r| median {np.nanmedian(np.abs(rb)):.3f}  p99 {np.nanpercentile(np.abs(rb),99):.3f}"
          f"  max {np.nanmax(np.abs(rb)):.3f}")
    print(f"  GPX6 r = {r_gpx:+.3f}  ->  {pct:.2f}nd percentile, empirical p = {emp:.4f}")
    print(f"\n  -> top 0.1%, but 0.1% of {len(rows):,} testable variants is ~{emp*len(rows):.0f}")
    print("     variants genome-wide. Strong, NOT genome-wide significant on its own.")

    # -- 7 locus attribution --------------------------------------------------------------
    hdr(7, "LEAD OR PASSENGER? neighbours within +-25 kb")
    ch, pos_ = idx["chrom"].astype(str), idx["pos"]
    sel = np.where((ch == CHROM) & (pos_ >= POS - WIN) & (pos_ <= POS + WIN))[0]
    sel = sel[np.minimum(p0[sel], 1 - p0[sel]) > 0.05]
    rn = garden_r(M, p0, sel, site, pc1)
    print(f"  non-SNP: rank {int((rn < r_gpx).sum())+1} of {len(sel)}")
    Ms = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_snp_af.npy", mmap_mode="r")
    p0s = np.load(f"{lib.AF_STORE}/p0_snp.npy")
    i2 = np.load(f"{lib.AF_STORE}/index_snp.npz", allow_pickle=True)
    c2, q2 = i2["chrom"].astype(str), i2["pos"]
    s2 = np.where((c2 == CHROM) & (q2 >= POS - WIN) & (q2 <= POS + WIN))[0]
    s2 = s2[np.minimum(p0s[s2], 1 - p0s[s2]) > 0.05]
    rs = garden_r(Ms, p0s, s2, site, pc1)
    nbeat = int((rs <= r_gpx).sum())
    print(f"  SNP    : {nbeat} of {len(s2)} SNPs match or beat it; best SNP r = {np.nanmin(rs):+.3f}"
          f" at {int(q2[s2][np.nanargmin(rs)]):,}")
    print("\n  -> lead among non-SNP variants, but NOT the best marker in the window:")
    print("     SNPs track the same haplotype at least as well. Consistent with the")
    print("     standing 'no kMate gain' result; the value is causal candidacy, not detection.")

    # -- 8 panel support ------------------------------------------------------------------
    hdr(8, "PANEL SUPPORT (call rate) -- everything downstream inherits this")
    vcf = f"{PROJ}/panel/arch3/chr4/merged_231_chr4_final.vcf.gz"
    q = subprocess.run([BCF, "query", "-r", f"{CHROM}:{POS-1000000}-{POS+1000000}",
                        "-f", "%POS\t%REF\t%ALT[\t%GT]\n", vcf],
                       capture_output=True, text=True).stdout
    recs = []
    for line in q.splitlines():
        f = line.split("\t")
        rl, al = len(f[1]), len(f[2])
        gts = f[3:]
        cr = 100 * sum(x != "." for x in gts) / max(len(gts), 1)
        recs.append((int(f[0]), abs(al - rl), cr))
    R = pd.DataFrame(recs, columns=["pos", "size", "cr"])
    SV = R[R["size"] > 50]
    our = R[(R.pos == POS) & (R["size"].between(1160, 1168))].cr.iloc[0]
    print(f"  GPX6 deletion call rate: {our:.1f}% of 231 founders")
    print(f"  SV records within 1 Mb : median {SV.cr.median():.1f}%  "
          f"(n={len(SV):,});  GPX6 is at the {100*(SV.cr < our).mean():.1f}th percentile")
    print("\n  -> not merely 'typical for an SV' -- it is among the worst-called SVs in the")
    print("     region. r2 to SNPs is untestable, so the eQTL track cannot be read, and the")
    print("     pool AF rests on a panel where two thirds of founders have no call.")


if __name__ == "__main__":
    main()

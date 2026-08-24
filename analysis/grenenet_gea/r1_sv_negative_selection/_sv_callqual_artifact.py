#!/usr/bin/env python
"""Calling-confidence proxy for the insertion-polarity artifact-vs-biology question
(2026-07-03 session). No outgroup genome/alignment exists in this repo, so true
ancestral polarization isn't available this session -- instead we ask: does the
SV-insertion climate-purging signal (SV_TEMPORAL_PURGING_SUMMARY.md) depend on
call quality, using the founder panel's OWN Minigraph-Cactus per-record quality
fields (F_MISSING = fraction of the 231 founders missing a genotype; MA = alleles
missing in panel haplotypes) as a confidence proxy. CONFLICT was also pulled but
is never populated in this VCF (0/2.25M records) -- dropped.

  (A) Genome-wide: do insertions have worse call quality than deletions among common
      SVs? (independent of the temporal test -- if insertions are just harder to call,
      that alone is a flag.)
  (B) Does the climate-slope beta = d s / d {bio1,bio18} purging excess (insertions vs
      matched SNP, from _compute_s_climate_slope.py's saved arrays) survive when
      restricted to WELL-supported insertions, or is it concentrated in poorly-
      supported ones (call-quality artifact signature)?

Note: the raw VCF has biallelic-split records, so (chrom,pos,ref_len,alt_len) is not
always a unique key (~5.6% of chr1 rows share a length-key with a same-position,
different-sequence allele) -- duplicates are mean-aggregated before joining, matching
how kMate's TSV already collapses to one row per (locus,ref_len,alt_len) (lib.rec_key).

Env: kmate. Reads analysis/grenenet_gea/sv_adaptive/{s_climate_slope.npz,vcf_callqual_chr*.tsv}.
Writes analysis/grenenet_gea/sv_adaptive/sv_callqual_artifact.npz + prints a text summary.
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

os.chdir("/global/scratch/users/tbellg/kmate")
STORE = lib.AF_STORE
MIN_P0, SV_BP, MIN_MAC, NBIN, NDRAW = 0.02, 50, 12, 25, 20000
rng = np.random.default_rng(0)


def match_to(tp0, ss, sp0):
    """Same frequency-matching scheme as _compute_s_climate_slope.py."""
    e = np.quantile(np.concatenate([tp0, sp0]), np.linspace(0, 1, NBIN + 1)); e[-1] += 1e-9
    tb = np.digitize(tp0, e[1:-1]); sb = np.digitize(sp0, e[1:-1]); out = []
    for b in range(NBIN):
        k = int(round((tb == b).mean() * NDRAW)); pool = ss[sb == b]
        if k > 0 and pool.size: out.append(rng.choice(pool, k, replace=True))
    return np.concatenate(out)


def load_callqual() -> pd.DataFrame:
    """Concat the per-chrom bcftools-extracted call-quality TSVs (non-SNP records),
    mean-aggregated over duplicate (chrom,pos,ref_len,alt_len) length-keys."""
    cols = ["chrom", "pos", "ref_len", "alt_len", "f_missing", "conflict", "ma"]
    fs = sorted(glob.glob(f"{lib.GEA}/sv_adaptive/vcf_callqual_chr*.tsv"))
    if not fs:
        raise FileNotFoundError("run _extract_vcf_callqual.sh first")
    d = pd.concat([pd.read_csv(f, sep="\t", header=None, names=cols) for f in fs],
                  ignore_index=True)
    d["f_missing"] = pd.to_numeric(d["f_missing"], errors="coerce")
    d["ma"] = pd.to_numeric(d["ma"], errors="coerce").fillna(0.0)
    key = (d.chrom.astype(str) + ":" + d.pos.astype(str) + ":"
           + d.ref_len.astype(str) + ":" + d.alt_len.astype(str))
    n_before = len(d)
    d = d.groupby(key)[["f_missing", "conflict", "ma"]].mean()
    print(f"[callqual] {n_before:,} raw non-SNP VCF rows -> {len(d):,} unique length-keys "
          f"({n_before - len(d):,} biallelic-split length-collisions mean-aggregated)")
    return d


def main():
    cq = load_callqual()
    print(f"[callqual] F_MISSING nan={cq.f_missing.isna().sum():,}, "
          f"CONFLICT=1 frac={cq.conflict.mean():.4f} (uninformative, dropped below), "
          f"MA==0 frac={(cq.ma == 0).mean():.4f}")

    # ---- (A) genome-wide: common-SV insertions vs deletions, call quality ----
    idx_non = np.load(f"{STORE}/index_nonsnp.npz")
    ch_non = idx_non["chrom"].astype("U5"); pos_non = idx_non["pos"].astype(np.int64)
    rl = idx_non["ref_len"].astype(np.int64); al = idx_non["alt_len"].astype(np.int64)
    dlen = np.abs(al - rl)
    common_non = lib.founder_panel_keep(ch_non, pos_non, min_mac=MIN_MAC)
    sv_all = common_non & (dlen > SV_BP)
    isdel_all = rl > al
    key_all = pd.Series(ch_non).astype(str) + ":" + pd.Series(pos_non).astype(str) + \
              ":" + pd.Series(rl).astype(str) + ":" + pd.Series(al).astype(str)
    qA = cq.reindex(key_all[sv_all].to_numpy())
    insA = ~isdel_all[sv_all]
    print("\n" + "=" * 70)
    print(f"(A) genome-wide common SVs (MAC>={MIN_MAC}, n={int(sv_all.sum()):,}, "
          f"join hit rate={qA.f_missing.notna().mean():.3f}): "
          f"ins={int(insA.sum()):,} del={int((~insA).sum()):,}")
    fm_ins, fm_del = qA.f_missing[insA], qA.f_missing[~insA]
    ma_ins, ma_del = qA.ma[insA], qA.ma[~insA]
    print(f"  F_MISSING median: ins={np.nanmedian(fm_ins):.4f} del={np.nanmedian(fm_del):.4f}  "
          f"MWU p={stats.mannwhitneyu(fm_ins.dropna(), fm_del.dropna()).pvalue:.2e}")
    print(f"  MA median:        ins={np.nanmedian(ma_ins):.1f} del={np.nanmedian(ma_del):.1f}  "
          f"(frac MA>0: ins={np.nanmean(ma_ins > 0):.3f} del={np.nanmean(ma_del > 0):.3f})  "
          f"MWU p={stats.mannwhitneyu(ma_ins.dropna(), ma_del.dropna()).pvalue:.2e}")

    # ---- (B) does the beta purging signal survive within quality tiers? ----
    print("\n" + "=" * 70)
    print("(B) climate-slope beta, insertions only, split by call quality")
    z = np.load(f"{lib.GEA}/sv_adaptive/s_climate_slope.npz")
    # recover chrom/pos for the SAME cols_non/sv_k ordering used to build z (deterministic,
    # matches _compute_s_climate_slope.py / _audit_s_classes.py exactly)
    p0_non = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    non_m = (dlen >= 1) & common_non & (p0_non >= MIN_P0) & (p0_non <= 1 - MIN_P0)
    cols_non = np.where(non_m)[0]
    dl = dlen[cols_non]; sv_k = dl > SV_BP
    ch_sv = ch_non[cols_non][sv_k]; pos_sv = pos_non[cols_non][sv_k]
    rl_sv = rl[cols_non][sv_k]; al_sv = al[cols_non][sv_k]
    key_sv = pd.Series(ch_sv).astype(str) + ":" + pd.Series(pos_sv).astype(str) + \
             ":" + pd.Series(rl_sv).astype(str) + ":" + pd.Series(al_sv).astype(str)
    assert key_sv.size == z["p0_sv"].size, "chrom/pos recompute count != saved beta array size"
    qB = cq.reindex(key_sv.to_numpy())
    sv_isdel = z["sv_isdel"]
    ins_mask = ~sv_isdel
    print(f"  n SV total={sv_isdel.size}  insertions={int(ins_mask.sum())}  "
          f"(join hit rate: {qB.f_missing.notna().mean():.3f})")

    well = (qB.f_missing.to_numpy() <= np.nanmedian(qB.f_missing)) & (qB.ma.to_numpy() == 0)
    poor = ~well
    print(f"  well-supported (F_MISSING<=median & MA==0): n={int(well.sum())}/{well.size}")
    for k in ("bio1", "bio18"):
        bsv = z[f"beta_sv_{k}"]; bsnp = z[f"beta_snp_{k}"]
        bmatch = match_to(z["p0_sv"], bsnp, z["p0_snp"])
        print(f"\n  -- beta = d s / d {k} --")
        print(f"     ALL insertions (n={int(ins_mask.sum())}): median={np.median(bsv[ins_mask]):+.5f}  "
              f"KS vs matched-SNP p={stats.ks_2samp(bsv[ins_mask], bmatch).pvalue:.2e}")
        for tag, m in (("well-supported", well), ("poorly-supported", poor)):
            sel = ins_mask & m
            if sel.sum() < 20:
                print(f"     {tag:>17} insertions: n={int(sel.sum())} (too few, skip)")
                continue
            bmatch_tier = match_to(z["p0_sv"][sel], bsnp, z["p0_snp"])
            print(f"     {tag:>17} insertions (n={int(sel.sum())}): "
                  f"median={np.median(bsv[sel]):+.5f}  "
                  f"KS vs matched-SNP p={stats.ks_2samp(bsv[sel], bmatch_tier).pvalue:.2e}")
        # is call quality itself confounded with p0 (already freq-matched, but check)?
        r = stats.spearmanr(qB.f_missing.to_numpy()[ins_mask], z["p0_sv"][ins_mask],
                            nan_policy="omit").statistic
        print(f"     Spearman(F_MISSING, p0) among insertions: {r:+.2f}")
        # continuous check: if the purging signal were a call-quality artifact, worse
        # quality (higher F_MISSING / MA) should predict a MORE NEGATIVE beta.
        rf, pf = stats.spearmanr(qB.f_missing.to_numpy()[ins_mask], bsv[ins_mask], nan_policy="omit")
        rm, pm = stats.spearmanr(qB.ma.to_numpy()[ins_mask], bsv[ins_mask], nan_policy="omit")
        print(f"     Spearman(F_MISSING, beta) among insertions: {rf:+.3f} (p={pf:.2g})  "
              f"Spearman(MA, beta): {rm:+.3f} (p={pm:.2g})  [artifact predicts negative r]")

    np.savez_compressed(f"{lib.GEA}/sv_adaptive/sv_callqual_artifact.npz",
                        key_sv=key_sv.to_numpy(), f_missing=qB.f_missing.to_numpy(),
                        ma=qB.ma.to_numpy(), well=well, ins_mask=ins_mask,
                        beta_sv_bio1=z["beta_sv_bio1"], beta_sv_bio18=z["beta_sv_bio18"],
                        p0_sv=z["p0_sv"])
    print("\n[wrote] sv_callqual_artifact.npz")


if __name__ == "__main__":
    main()

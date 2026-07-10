#!/usr/bin/env python3
"""Score per-SV alt-AF for kMate and vg against the independent SV truth (§3).

Companion to score_competitors.py (SNPs). The SV panel (BENCHMARK_DESIGN.md §3)
compares:
  - kMate  : direct AF for SV records, subset from the per-record TSV (|indel|>=SVLEN).
  - vg     : graph-native AF = AD_alt/(AD_ref+AD_alt) from a `vg call` VCF.
             Preferred input is `vg call -v <panel_sv.vcf>` (genotype-given-variants),
             which emits one record per PANEL SV -> exact key join, the apples-to-apples
             analog of how the SNP benchmark scores vg at the shared panel sites.
             A de-novo `vg call` VCF can be matched fuzzily with --vg-denovo (sparse;
             vg discovers few SVs at low coverage -- reported as a sanity artifact).
  - hapFIRE: BLANK by construction -- HARP's likelihood is SNP-base only and has no
             SV term (§3.1), so it has no native SV output. Always emitted as an
             explicit NA row so the table states this rather than omitting it.

All tools are joined to the SAME SV truth on the panel key (chrom,pos,ref_len,alt_len),
so metrics are directly comparable. Reports MAE + RMSE + R2 (kMate convention) AND
Pearson r (freqk-paper convention), against the physical truth (missing->REF, fair to
read-based tools) and on the convention-free fully-called subset (n_called==F).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # benchmarks/
from _accuracy_panel import grid_panel  # noqa: E402


def norm_chrom(s):
    return s.astype(str).str.replace("Chr", "", regex=False).str.replace("chr", "", regex=False)


def load_truth(p, svlen):
    """SV truth, in panel-meta order. (pos,ref_len,alt_len) is NOT a unique key for
    SVs (two ALT sequences of equal length can share a position), so the panel row
    INDEX over the SV subset is the join key. truth_sv is already that ordered subset."""
    t = pd.read_csv(p, sep="\t")
    t = t[(t.alt_len - t.ref_len).abs() >= svlen].reset_index(drop=True)
    t["svidx"] = np.arange(len(t))
    return t


def load_kmate(p, svlen):
    """kMate per-record output is the full panel meta in the SAME order as truth.
    Apply the SAME SV mask and index positionally -> svidx aligns to truth's svidx."""
    d = pd.read_csv(p, sep="\t")
    d = d[(d.alt_len - d.ref_len).abs() >= svlen].reset_index(drop=True)
    d["svidx"] = np.arange(len(d))
    return d.rename(columns={"alt_freq": "est"})[["svidx", "pos", "ref_len", "alt_len", "est"]]


def _vcf_records(p):
    """Yield (chrom, pos, ref, alt, ad_ref, ad_alt) per ALT allele of a vg-call VCF."""
    with open(p) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            chrom, pos, ref, alts = f[0], int(f[1]), f[3], f[4].split(",")
            fmt = f[8].split(":")
            if "AD" not in fmt:
                continue
            ad = f[9].split(":")[fmt.index("AD")].split(",")
            try:
                ad = [int(x) for x in ad]
            except ValueError:
                continue
            ref_ad = ad[0] if ad else 0
            for i, alt in enumerate(alts):
                if alt in (".", "*", "<NON_REF>"):
                    continue
                alt_ad = ad[i + 1] if i + 1 < len(ad) else 0
                yield chrom, pos, ref, alt, ref_ad, alt_ad


def seq_key(chrom, pos, ref, alt):
    return f"{chrom.replace('Chr', '').replace('chr', '')}:{pos}:{ref}:{alt}"


def build_panel_index(panel_sv_vcf, svlen):
    """Map each panel SV's exact (chrom,pos,REF,ALT) sequence -> svidx (its row in
    the SV subset, in VCF order == truth/kMate order). Lets vg's own records be
    placed into panel-index space by sequence identity, not by fragile pos+length."""
    key2idx, idx = {}, 0
    op = (lambda f: __import__("gzip").open(f, "rt")) if str(panel_sv_vcf).endswith(".gz") else open
    with op(panel_sv_vcf) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            chrom, pos, ref, alts = f[0], int(f[1]), f[3], f[4].split(",")
            for alt in alts:
                if abs(len(alt) - len(ref)) < svlen:
                    continue
                key2idx[seq_key(chrom, pos, ref, alt)] = idx
                idx += 1
    return key2idx


def load_vg_genotyped(p, svlen, key2idx):
    """vg call -v <panel_sv.vcf>: AF = AD_alt/(AD_ref+AD_alt), placed at the panel
    svidx by exact (chrom,pos,REF,ALT) sequence match against the panel SV VCF."""
    rows = []
    for chrom, pos, ref, alt, ra, aa in _vcf_records(p):
        if abs(len(alt) - len(ref)) < svlen:
            continue
        sv = key2idx.get(seq_key(chrom, pos, ref, alt))
        if sv is None:
            continue
        dp = ra + aa
        rows.append((sv, aa / dp if dp > 0 else np.nan))
    d = pd.DataFrame(rows, columns=["svidx", "est"])
    return d.groupby("svidx", as_index=False).agg(est=("est", "first"))


def load_vg_denovo(truth, p, svlen, pos_tol=25, len_frac=0.25):
    """De-novo vg call: no panel key. Match each panel SV to the nearest vg SV
    (same svlen sign, |dpos|<=pos_tol, |dlen|/len<=len_frac); best = smallest dpos.
    Returns est aligned to truth['k'] (NaN where no vg SV matched)."""
    vg = [(pos, len(alt) - len(ref), (aa / (ra + aa) if (ra + aa) > 0 else np.nan))
          for _, pos, ref, alt, ra, aa in _vcf_records(p)
          if abs(len(alt) - len(ref)) >= svlen]
    vg_arr = np.array([(p_, l_) for p_, l_, _ in vg], dtype=float)
    vg_af = np.array([a_ for _, _, a_ in vg])
    est = np.full(len(truth), np.nan)
    if len(vg) == 0:
        return pd.DataFrame({"svidx": truth["svidx"].values, "est": est})
    tpos = truth["pos"].values.astype(float)
    tlen = (truth["alt_len"] - truth["ref_len"]).values.astype(float)
    for i in range(len(truth)):
        dpos = np.abs(vg_arr[:, 0] - tpos[i])
        m = (dpos <= pos_tol) & (np.sign(vg_arr[:, 1]) == np.sign(tlen[i]))
        denom = max(abs(tlen[i]), 1.0)
        m &= (np.abs(vg_arr[:, 1] - tlen[i]) / denom <= len_frac)
        if m.any():
            j = np.where(m)[0][np.argmin(dpos[m])]
            est[i] = vg_af[j]
    return pd.DataFrame({"svidx": truth["svidx"].values, "est": est})


def metrics(truth, est):
    m = np.isfinite(truth) & np.isfinite(est)
    truth, est = truth[m], est[m]
    if len(truth) == 0:
        return dict(n=0, MAE=np.nan, RMSE=np.nan, R2=np.nan, pearson_r=np.nan)
    err = est - truth
    mae = float(np.abs(err).mean())
    rmse = float(np.sqrt((err ** 2).mean()))
    ss_res = float((err ** 2).sum())
    ss_tot = float(((truth - truth.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    pear = float(np.corrcoef(truth, est)[0, 1]) if len(truth) > 1 else np.nan
    return dict(n=int(m.sum()), MAE=mae, RMSE=rmse, R2=r2, pearson_r=pear)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, help="truth_sv_*.tsv")
    ap.add_argument("--kmate")
    ap.add_argument("--vg", help="vg call -v <panel_sv.vcf> output (genotype-given-variants)")
    ap.add_argument("--panel-sv-vcf", help="panel SV VCF (seq->svidx map; required with --vg)")
    ap.add_argument("--vg-denovo", help="de-novo vg call VCF (fuzzy match; sanity only)")
    ap.add_argument("--vg-cov", help="Option-2 pack-coverage AF table (svidx<TAB>est) from vg_cov_sv.py")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--F", type=int, default=80, help="#founders (for fully-called subset)")
    ap.add_argument("--svlen", type=int, default=50, help="min |alt_len-ref_len| to call SV")
    ap.add_argument("--fullcalled", action="store_true",
                    help="restrict to fully-called SVs (n_called==F): convention-free")
    ap.add_argument("--truth-col", default="truth_af_phys",
                    choices=["truth_af_phys", "truth_af"],
                    help="physical (missing->REF, fair to read-based tools) vs MAR")
    a = ap.parse_args()

    tr = load_truth(a.truth, a.svlen)
    tr["truth_af"] = tr[a.truth_col]
    print(f"[score_sv] {len(tr):,} SV truth records (|indel|>={a.svlen}bp), "
          f"truth column: {a.truth_col}", file=sys.stderr)
    if a.fullcalled:
        keep = tr["n_called"] == a.F
        print(f"[score_sv] fully-called (n_called=={a.F}): {keep.sum():,}/{len(tr):,} kept",
              file=sys.stderr)
        tr = tr[keep]

    rows, cells = [], []

    def record(name, truth_v, est_v):
        mt = metrics(truth_v, est_v)
        mt["tool"] = name
        rows.append(mt)
        cells.append((name, truth_v, est_v))
        print(f"{name}: n={mt['n']:,} MAE={mt['MAE']:.4f} RMSE={mt['RMSE']:.4f} "
              f"R2={mt['R2']:.4f} r={mt['pearson_r']:.4f}", file=sys.stderr)

    if a.kmate:
        km = load_kmate(a.kmate, a.svlen)
        j = tr.merge(km, on="svidx", how="inner", suffixes=("", "_km"))
        bad = ((j.pos != j.pos_km) | (j.ref_len != j.ref_len_km)
               | (j.alt_len != j.alt_len_km))
        if bad.any():
            sys.exit(f"[score_sv] FATAL: {bad.sum()} positional mismatches kMate vs truth "
                     "-- panel order differs, cannot align SVs by index.")
        record("kMate", j["truth_af"].values, j["est"].values)
    if a.vg:
        if not a.panel_sv_vcf:
            sys.exit("[score_sv] --vg requires --panel-sv-vcf (seq->svidx map)")
        key2idx = build_panel_index(a.panel_sv_vcf, a.svlen)
        vg = load_vg_genotyped(a.vg, a.svlen, key2idx)
        j = tr.merge(vg, on="svidx", how="inner")
        record("vg-giraffe", j["truth_af"].values, j["est"].values)
    if a.vg_denovo:
        vg = load_vg_denovo(tr, a.vg_denovo, a.svlen)  # already aligned to tr rows
        record("vg-denovo", tr["truth_af"].values, vg["est"].values)
    if a.vg_cov:
        vc = pd.read_csv(a.vg_cov, sep="\t")
        j = tr.merge(vc, on="svidx", how="inner")
        record("vg-cov", j["truth_af"].values, j["est"].values)

    # hapFIRE: explicit blank (no native SV estimation, §3.1)
    rows.append(dict(tool="hapFIRE", n=0, MAE=np.nan, RMSE=np.nan, R2=np.nan, pearson_r=np.nan))

    res = pd.DataFrame(rows)[["tool", "n", "MAE", "RMSE", "R2", "pearson_r"]]
    res.to_csv(a.out_prefix + "_metrics.tsv", sep="\t", index=False)
    if cells:
        grid_panel(cells, "", a.out_prefix + "_scatter.png", ncols=len(cells), cell=4.0)
    print(res.to_string(index=False))

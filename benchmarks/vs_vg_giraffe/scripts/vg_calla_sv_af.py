#!/usr/bin/env python3
"""Extract vg-SV allele frequency from a `vg call -a` VCF (freqk recipe) and map to
the panel svidx for SV scoring.

Per record (freqk build_per_sample_p_parquet_vg.py): AF = AD_alt / sum(AD), split
multi-allelic into one row per ALT. Match each vg ALT to a panel SV by exact
(chrom,pos,REF,ALT) sequence, or SWAP (chrom,pos,ALT,REF) -> use 1-AF. Emits a
(svidx, est) table for score_sv.py --vg-cov.

Usage:
  vg_calla_sv_af.py --vg-vcf <vgcalla.vcf> --panel-sv-vcf <panel_sv.vcf.gz> \
      --out <svidx_est.tsv> [--svlen 50] [--keep-filter PASS,lowdepth,lowad]
"""
import argparse, gzip, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_sv import build_panel_index, seq_key  # noqa: E402


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vg-vcf", required=True)
    ap.add_argument("--panel-sv-vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--svlen", type=int, default=50)
    ap.add_argument("--keep-filter", default="PASS,lowdepth,lowad,.",
                    help="FILTER values to keep (default: keep all common vg labels)")
    a = ap.parse_args()

    key2idx = build_panel_index(a.panel_sv_vcf, a.svlen)   # exact seq_key -> svidx
    keep = set(a.keep_filter.split(","))

    best = {}                                  # svidx -> est (prefer exact over swap)
    n_rec = n_sv = n_exact = n_swap = 0
    with opener(a.vg_vcf) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, ref, altf, filt, fmt, samp = f[0], int(f[1]), f[3], f[4], f[6], f[8], f[9]
            if filt not in keep:
                continue
            fd = dict(zip(fmt.split(":"), samp.split(":")))
            ad = fd.get("AD", "")
            if not ad or ad == ".":
                continue
            try:
                adv = [int(x) for x in ad.split(",")]
            except ValueError:
                continue
            tot = sum(adv)
            if tot <= 0 or len(adv) < 2:
                continue
            n_rec += 1
            for i, alt in enumerate(altf.split(",")):
                if alt in (".", "*") or abs(len(alt) - len(ref)) < a.svlen:
                    continue
                n_sv += 1
                af = adv[i + 1] / tot if i + 1 < len(adv) else np.nan
                if not np.isfinite(af):
                    continue
                sv = key2idx.get(seq_key(chrom, pos, ref, alt))
                if sv is not None:                       # exact match
                    if sv not in best or best[sv][1] != "exact":
                        best[sv] = (af, "exact"); n_exact += 1
                    continue
                svs = key2idx.get(seq_key(chrom, pos, alt, ref))  # swap (REF<->ALT)
                if svs is not None and svs not in best:
                    best[svs] = (1.0 - af, "swap"); n_swap += 1

    rows = [(sv, est) for sv, (est, _) in best.items()]
    out = pd.DataFrame(rows, columns=["svidx", "est"]).sort_values("svidx")
    out.to_csv(a.out, sep="\t", index=False)
    print(f"[vg call -a] records w/ AD={n_rec:,} | SV-alt alleles={n_sv:,} | "
          f"matched panel SVs: exact={n_exact:,} swap={n_swap:,} total={len(out):,} "
          f"of {len(key2idx):,} panel SVs -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()

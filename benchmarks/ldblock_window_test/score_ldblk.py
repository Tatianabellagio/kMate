#!/usr/bin/env python
"""Score window-kMate output vs sim truth: % windows locally fit (status==0) + AF accuracy
(MAE/R2/r on SNPs and all records). Two modes:
  --validate <out.tsv> <truth.gz> <hblocks.npz>   # one run, print metrics
  --sweep                                          # score all sweep_s*_*.tsv, aggregate per floor
"""
import os, sys, glob, argparse
import numpy as np, pandas as pd

SIMS = "benchmarks/p231/sims"
DIR = "benchmarks/ldblock_window_test"


def metrics(obs, tru):
    """Agent convention (plot_ldblock/_accuracy_panel): drop non-finite, then drop
    records that are 0 in BOTH truth and est (true-negative mass)."""
    m = np.isfinite(obs) & np.isfinite(tru); obs, tru = obs[m], tru[m]
    keep = (tru > 0) | (obs > 0); obs, tru = obs[keep], tru[keep]
    e = obs - tru
    sstot = np.sum((tru - tru.mean()) ** 2)
    r2 = 1 - np.sum(e ** 2) / sstot if sstot > 0 else np.nan
    r = np.corrcoef(obs, tru)[0, 1] if len(obs) > 2 else np.nan
    return float(np.mean(np.abs(e))), float(r2), float(r), len(obs)


def score(out_tsv, truth_gz, hblocks_npz):
    o = pd.read_csv(out_tsv, sep="\t")
    t = pd.read_csv(truth_gz, sep="\t")                # row-aligned to est (NO merge)
    assert len(o) == len(t), (len(o), len(t))
    # defence-in-depth: confirm genuine key-alignment, not just same length
    for k in ("pos", "ref_len", "alt_len"):
        assert (o[k].values == t[k].values).all(), f"{k} mismatch: est and truth not row-aligned"
    e = o.alt_freq.values.astype(np.float64); tr = t.truth_af.values.astype(np.float64)
    snp = ((o.ref_len == 1) & (o.alt_len == 1)).values
    mae_s, r2_s, r_s, n_s = metrics(e[snp], tr[snp])
    mae_a, r2_a, r_a, n_a = metrics(e, tr)
    st = np.load(hblocks_npz, allow_pickle=True)["Chr1_status"]
    return dict(n_snp=n_s, n_blocks=len(st),
                pct_local=round(100 * (st == 0).mean(), 1),
                pct_fallback=round(100 * (st == 1).mean(), 1),
                pct_empty=round(100 * (st == 2).mean(), 1),
                mae_snp=round(mae_s, 4), r2_snp=round(r2_s, 4), r_snp=round(r_s, 4),
                mae_all=round(mae_a, 4), r2_all=round(r2_a, 4))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", nargs=3, metavar=("OUT", "TRUTH", "HBLOCKS"))
    ap.add_argument("--sweep", action="store_true")
    a = ap.parse_args()
    if a.validate:
        print(score(*a.validate))
        return
    if a.sweep:
        rows = []
        for f in sorted(glob.glob(f"{DIR}/sweep_s*_*.tsv")):
            base = os.path.basename(f)[:-4]                  # sweep_s42_8
            _, sseed, level = base.split("_")
            seed = int(sseed[1:])
            pool = f"cov10_n50_g3_s{seed}_self97_hotspots_dom500_p231_chr1"
            truth = f"{SIMS}/{pool}/recomb_truth_raw.tsv.gz"   # panel-aligned, row-by-row
            hb = f"{DIR}/{base}.h_blocks_per_chrom.npz"
            if not (os.path.exists(truth) and os.path.exists(hb)):
                continue
            r = score(f, truth, hb); r["seed"] = seed; r["level"] = level
            rows.append(r)
        if not rows:
            print("no sweep outputs scored yet"); return
        D = pd.DataFrame(rows)
        D.to_csv(f"{DIR}/coarse_sweep_scores.csv", index=False)
        order = ["base", "8", "15", "25", "40"]
        print(f"{'level':>6}{'nSeed':>6}{'%local':>9}{'%empty':>8}{'R2_snp':>9}{'MAE_snp':>9}{'R2_all':>9}")
        for lv in order:
            d = D[D.level == lv]
            if len(d):
                print(f"{lv:>6}{len(d):>6}{d.pct_local.mean():>8.1f}%{d.pct_empty.mean():>7.1f}%"
                      f"{d.r2_snp.mean():>9.4f}{d.mae_snp.mean():>9.4f}{d.r2_all.mean():>9.4f}"
                      f"   (R2_snp sd {d.r2_snp.std():.4f})")
        print(f"\nwrote {DIR}/coarse_sweep_scores.csv ({len(D)} runs)")


if __name__ == "__main__":
    main()

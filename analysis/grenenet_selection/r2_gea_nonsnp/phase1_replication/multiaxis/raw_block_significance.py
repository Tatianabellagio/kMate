#!/usr/bin/env python
"""How much of the genome is significant in the RAW (pre-WZA) per-record GEA?

Counterpart to class_peaks_overlap_table.py, which reports before/after-WZA counts.
This one never touches WZA: it takes the per-record model p-values that WZA *would*
aggregate, calls a clq0.9 tiling block a HIT when >=1 of its records of that class
clears the bar, and converts hit blocks to base pairs -- so the output answers
"what fraction of the genome does this scan call" rather than only "how many hits".

Reported per (class x axis): significant records, hit blocks, Mb, % of genome, and
% of that class's own testable footprint. Bonferroni and BH-FDR side by side, each
on the raw p and on its GIF-corrected twin, so the four thresholding regimes are
comparable in one table without rerunning anything.

CAVEATS THE OUTPUT CANNOT CARRY ITSELF
--------------------------------------
1. Raw p is UNCALIBRATED. Genomic inflation on this data is lambda ~1.7 median and
   up to 2.7 (bio15/bio10/bio1); 19/20 axes have lambda>1. The `lam` column is in
   the output for exactly this reason -- report it beside any raw count. BH on raw
   p is not a candidate list: it calls 40-60% of the genome on the strong axes.
2. Hit blocks are MUCH larger than average, so the Mb/% columns are size-driven.
   The tiling block size distribution is severely right-skewed: median 756 bp,
   mean 2.04 kb, p90 4.5 kb, p99 18.1 kb, max 1.22 Mb. Bonferroni-hit blocks run
   8-28x the mean (snp bio1: 188 blocks = 4.60 Mb = 24.5 kb/block, i.e. sitting
   around the 99th percentile of block size). Bigger blocks hold more records and
   so get more chances to clear a per-record bar -- and because a handful of blocks
   are hundreds of kb, hitting one or two of them moves "% of genome" by whole
   percentage points on its own. The block-count column is therefore the honest
   headline; treat Mb/% as an upper bound. The `*_kb_per_blk` columns exist to make
   this visible per cell. This size confound is exactly what WZA's SNP-number
   correction absorbs -- dropping WZA does not remove it, it relocates it to the
   reader.
3. GIF is applied only when lambda>1 (matching gif_manhattan_*_tile.ipynb):
   dividing chi2 by lambda<1 would manufacture signal.

SPAN RULE
---------
Each block id owns exactly its own tiling interval (end_{i-1}, end_i]; block 0 owns
[0, end_0] and the last owns (end_{n-2}, chrom_len]. Absent block ids are NOT
absorbed by their neighbours: reblock_blockdef.py applies merge_small_blocks per
class, so an id missing from a class mixes "merged into the previous block" with
"this class has no record here", and the two are indistinguishable afterwards.
Not absorbing keeps spans additive and bounded by the genome; the cost is that a
merged block's bp is under-attributed. For sparse classes this is why `testable_mb`
is far below the genome (sv occupies ~4.5k of 58,376 blocks) -- for a class-contrast
claim use the *_pcttest columns, not *_pctgen.

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY phase1_replication/multiaxis/raw_block_significance.py --model lfmm
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import lib

MA = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis"
IND = f"{MA}/wza_in_clq09_tile"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
CLASSES = ["snp", "sv", "smallindel", "nonsnp"]
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
CHROM_LEN = {"Chr1": 30427671, "Chr2": 19698289, "Chr3": 23459830,
             "Chr4": 18585056, "Chr5": 26975502}          # TAIR10
GENOME = sum(CHROM_LEN.values())
EXP_MED = stats.chi2.isf(0.5, 1)


def bh(p):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    q[o] = p[o] * n / (np.arange(n) + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)


def gif_correct(p):
    """(p_used, lambda, applied). Deflate chi2 by lambda only when lambda>1."""
    p = np.clip(np.asarray(p, float), 1e-300, 1.0)
    c2 = stats.chi2.isf(p, 1)
    lam = float(np.median(c2) / EXP_MED)
    if lam > 1.0:
        return np.clip(stats.chi2.sf(c2 / lam, 1), 1e-300, 1.0), lam, True
    return p, lam, False


def block_spans(r2=0.9):
    """bp owned by each clq{r2} tiling block id; sums to exactly the genome."""
    span = {}
    for ci, ch in enumerate(CHROMS, start=1):
        f = f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_clq{r2}_blocks_clq{r2}.tsv"
        e = pd.read_csv(f, sep="\t").sort_values("start_pos")["end_pos"].to_numpy(np.int64)
        lo = np.concatenate([[0], e[:-1]])
        hi = np.concatenate([e[:-1], [CHROM_LEN[ch]]])
        for k in range(len(e)):
            span[f"{ch}_{k}"] = int(hi[k] - lo[k])
    assert sum(span.values()) == GENOME, (sum(span.values()), GENOME)
    return span


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="lfmm", choices=["lfmm", "kendall", "binomial"])
    ap.add_argument("--maf", type=float, default=0.05)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = args.out or f"{MA}/raw_block_significance_{args.model}.csv"

    SPAN = block_spans()
    print(f"clq0.9 partition: {len(SPAN):,} blocks tiling {GENOME:,} bp", flush=True)

    rows, union = [], {}
    for cls in CLASSES:
        # cols 1-6 are byte-identical across axes for a class, so read the record
        # metadata once and pull only `pval` per axis.
        d0 = pd.read_csv(f"{IND}/{args.model}_{cls}_gen9_bio1.csv", usecols=["MAF", "block"])
        keep = (d0.MAF > args.maf).to_numpy()
        blk = d0.block.astype(str).to_numpy()[keep]
        n = int(keep.sum())
        tested = sorted(set(blk) - {"", "nan"})
        test_bp = sum(SPAN.get(b, 0) for b in tested)
        print(f"\n{cls}: {n:,} records (MAF>{args.maf}) in {len(tested):,}/{len(SPAN):,} blocks, "
              f"{test_bp/1e6:.1f} Mb testable ({100*test_bp/GENOME:.1f}% of genome)", flush=True)

        for axis in AXES:
            f = f"{IND}/{args.model}_{cls}_gen9_{axis}.csv"
            if not os.path.exists(f):
                print(f"  {axis:5s} MISSING"); continue
            p_raw = pd.read_csv(f, usecols=["pval"]).pval.to_numpy(float)[keep]
            p_raw = np.clip(np.nan_to_num(p_raw, nan=1.0), 1e-300, 1.0)
            p_gif, lam, applied = gif_correct(p_raw)
            r = dict(model=args.model, cls=cls, axis=axis, n_records=n,
                     n_blocks_tested=len(tested), testable_mb=round(test_bp / 1e6, 2),
                     lam=round(lam, 3), gif_applied=applied)
            for tag, p in (("raw", p_raw), ("gif", p_gif)):
                for bar, sig in (("bonf", p < args.alpha / n), ("bh", bh(p) < args.alpha)):
                    hb = sorted(set(blk[sig]) - {"", "nan"})
                    bp = sum(SPAN.get(b, 0) for b in hb)
                    r[f"{tag}_{bar}_rec"] = int(sig.sum())
                    r[f"{tag}_{bar}_blk"] = len(hb)
                    r[f"{tag}_{bar}_mb"] = round(bp / 1e6, 2)
                    r[f"{tag}_{bar}_pctgen"] = round(100 * bp / GENOME, 3)
                    r[f"{tag}_{bar}_pcttest"] = round(100 * bp / test_bp, 3) if test_bp else 0.0
                    # mean size of a hit block, the size-confound diagnostic
                    r[f"{tag}_{bar}_kb_per_blk"] = round(bp / len(hb) / 1e3, 2) if hb else 0.0
                    union.setdefault((cls, tag, bar), set()).update(hb)
            rows.append(r)
            print(f"  {axis:5s} lam={lam:4.2f} | RAW bonf {r['raw_bonf_blk']:>6,}blk "
                  f"{r['raw_bonf_pctgen']:>6.2f}%gen {r['raw_bonf_kb_per_blk']:>6.1f}kb/blk "
                  f"| RAW bh {r['raw_bh_blk']:>6,}blk {r['raw_bh_pctgen']:>6.2f}%gen "
                  f"|| GIF bonf {r['gif_bonf_blk']:>4,}blk  GIF bh {r['gif_bh_blk']:>5,}blk",
                  flush=True)

    R = pd.DataFrame(rows)
    R.to_csv(out, index=False)

    mean_kb = GENOME / len(SPAN) / 1e3
    for lab, col in [("RAW Bonferroni : blocks hit", "raw_bonf_blk"),
                     ("RAW Bonferroni : % of genome", "raw_bonf_pctgen"),
                     ("RAW Bonferroni : kb per hit block "
                      f"(genome mean block = {mean_kb:.2f} kb)", "raw_bonf_kb_per_blk"),
                     ("RAW BH q<0.05 : blocks hit", "raw_bh_blk"),
                     ("RAW BH q<0.05 : % of genome", "raw_bh_pctgen"),
                     ("GIF Bonferroni : blocks hit", "gif_bonf_blk"),
                     ("GIF BH q<0.05 : blocks hit", "gif_bh_blk")]:
        print(f"\n===== {lab} =====")
        print(R.pivot(index="axis", columns="cls", values=col)
               .reindex(index=AXES, columns=CLASSES).to_string())

    # union across axes: a block counts once if ANY of the 20 axes hits it. This is
    # the "how much of the genome does the whole scan implicate" number -- the
    # per-axis rows cannot be summed, since the same block recurs across the
    # correlated bioclim axes.
    u = []
    for cls in CLASSES:
        tb = int(R.loc[R.cls == cls, "n_blocks_tested"].iloc[0])
        row = dict(cls=cls, blocks_tested=tb)
        for tag in ("raw", "gif"):
            for bar in ("bonf", "bh"):
                hb = union.get((cls, tag, bar), set())
                bp = sum(SPAN.get(b, 0) for b in hb)
                row[f"{tag}_{bar}_blk"] = len(hb)
                row[f"{tag}_{bar}_pct_of_tested"] = round(100 * len(hb) / tb, 2) if tb else 0.0
                row[f"{tag}_{bar}_mb"] = round(bp / 1e6, 2)
                row[f"{tag}_{bar}_pctgen"] = round(100 * bp / GENOME, 2)
        u.append(row)
    U = pd.DataFrame(u)
    U.to_csv(out.replace(".csv", "_union.csv"), index=False)
    print("\n===== UNION over the 20 axes (block hit by >=1 axis; NOT the sum of rows) =====")
    print(U[["cls", "blocks_tested", "raw_bonf_blk", "raw_bonf_pct_of_tested",
             "raw_bonf_mb", "raw_bonf_pctgen", "raw_bh_blk", "raw_bh_pctgen",
             "gif_bonf_blk", "gif_bh_blk"]].to_string(index=False))

    print("\n===== totals summed over 20 axes (record level) =====")
    print(R.groupby("cls")[["raw_bonf_rec", "raw_bh_rec", "gif_bonf_rec", "gif_bh_rec"]]
           .sum().reindex(CLASSES).to_string())
    print(f"\nwrote {out}")
    print(f"wrote {out.replace('.csv', '_union.csv')}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Consolidate the 2026-09-30 GEA re-run into one per-record/per-unit production table.

The re-run changed two things at once, so its outputs live in two places:

  snp, smallindel   RECORD level, results/multiaxis/lfmm_cr50/  (ma_lfmm_cr50.sbatch)
                    A SNP is one record; a short indel is one record per unit
                    (build_event_matrix.py: "the small-indel scan is unchanged at the unit
                    level and is not repeated"), so for these two classes record == unit and
                    the record-level fit IS the event-level fit. Verified: the record-level
                    smallindel scan and the non-SV part of the event nonsnp scan both hold
                    603,630 rows.
  sv, nonsnp        EVENT level, genes/convergence/event_gea/cr50/ (run_event_lfmm_cr50.sbatch)
                    arch3 splits one SV into one record per assembly path; these two scans
                    collapse those into truvari events, which is the unit the SV literature
                    uses and the unit the burden/purging work already used.

Both arms carry MAF >= 0.05 and the founder call-rate floor >= 0.50. Both come out of LFMM
ridge K=16 with RAW p (no GIF, project rule). What they do NOT share is the block column:
both write hapFIRE LD blocks, while every downstream notebook keys on the clq r2=0.9 TILING
partition. This script reassigns that partition and writes one directory with one schema.

BLOCK RULE -- and why this differs from the old wza_in_clq09_tile
-----------------------------------------------------------------
The clq r2=0.9 partition is defined ONCE, by BigLD, on a dense common-variant set (MAF > 0.05,
call rate >= 90%, one variant per position -- 291k variants on Chr1 alone, overwhelmingly
SNPs), and `blocks_tiling.assign_tiling` then drops any variant of any class into the interval
that contains it. That is the property the cross-class comparisons depend on: a block id names
the same stretch of LD for a SNP and for an SV, so "this SV is in the same block as a
significant SNP" is a statement about linkage.

The old production tables destroyed that property. They applied `merge_small_blocks` on top,
which folds any block holding < 2 records OF THAT CLASS into its neighbour -- so the merge runs
per class and bites hardest where the class is sparsest. Measured on this re-run:

    class        stored ids   true blocks occupied   ids spanning >1 true block   records mis-keyed
    snp              56,832                 58,132                         2.2%                0.1%
    smallindel       41,496                 49,388                        15.7%                1.3%
    nonsnp           41,697                 49,505                        15.5%                1.3%
    sv                3,583                  8,955                        55.8%               26.7%

One merged SV id covers up to 24 real blocks. A quarter of SV records therefore carried a
block id naming a different interval than the same position's SNPs did, which both invents
"shared with SNP" overlaps and hides real ones -- and it is the mechanism behind the SV
new-peak gene mis-attribution already on record (CARK8/9).

That merge existed for HapFM/hapFIRE's block-level aggregation, where a block needs >= 2
markers to be worth aggregating. This project retired WZA: p-values are per record/unit and
Bonferroni is per class, and blocks are only ever used afterwards to ask which records sit
together. Nothing here aggregates within a block, so the merge buys nothing and costs
comparability. It is therefore OFF by default; `--merge` restores the old behaviour for
reproducing the archived tables.

Even unmerged, blocks are a coarse instrument for gene attribution (SV-occupied blocks are
sparse). Prefer loci-with-a-window (_build_peak_alignment_nb.py) where attribution matters.

Out: results/multiaxis/tile_cr50/lfmm_{cls}_gen9_{axis}.csv
     columns chrom,pos,ref_len,alt_len,MAF,block,pval  (the old schema, so existing readers
     work unchanged) + unit,dlen,n_alleles,is_sv,n_called,B,z where the source has them.
env: kmate. Compute node.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GS = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, GS)
sys.path.insert(0, os.path.join(GS, "blocks"))
import lib                                                        # noqa: E402
import blocks_tiling as bt                                        # noqa: E402

MA = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis"
REC = f"{MA}/lfmm_cr50{lib.VSUF}"                                   # every path is the panel version's (KMATE_PANEL)
EVT = f"{lib.GEA}/genes/convergence/event_gea/cr50{lib.VSUF}"
OUTDIR = f"{MA}/tile_cr50{lib.VSUF}"
CM = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices_cr50{lib.VSUF}"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1", "pc2", "pc3"]
CLASSES = ["snp", "smallindel", "sv", "nonsnp"]
HEAD = ["chrom", "pos", "ref_len", "alt_len", "MAF", "block", "pval"]
EXTRA = ["unit", "dlen", "n_alleles", "is_sv", "n_called", "B", "z"]


def smallindel_unit_ids() -> np.ndarray:
    """`unit` id per row of the record-level small-indel tables, matched on INTEGER ids.

    Short indels are one unit per record, so the small-indel scan and the non-SV part of the
    event scan cover the same records -- but in different row orders (event tables are ordered
    by unit id, record tables by af_store order). The class matrix's `col` is the record's
    index into af_store/index_nonsnp, and nonsnp_event_units.npz lists each unit's members as
    those same indices, so the two match exactly on integers. Carrying `unit` here lets
    build_gea_pool.py collapse on it instead of on (chrom,pos,ref_len,alt_len) -- 4.9% of
    small-indel units share such a tuple, and collapsing on it merges distinct variants.
    """
    E = np.load(f"{lib.SV_ADAPTIVE}/nonsnp_event_units.npz")
    U = pd.read_csv(f"{EVT}/nonsnp_gen9_units.csv", usecols=["unit", "is_sv", "n_alleles"])
    N = U[~U.is_sv]
    assert (N.n_alleles == 1).all(), "a non-SV unit holds >1 record; the 1:1 mapping is invalid"
    member = E["indices"][E["indptr"][N.unit.to_numpy()]]
    unit_of = pd.Series(N.unit.to_numpy(), index=member)
    rec = pd.read_csv(f"{CM}/smallindel_gen9.records.csv", usecols=["pos", "col"])
    out = unit_of.reindex(rec.col.to_numpy())
    assert out.notna().all(), f"{int(out.isna().sum()):,} small-indel records have no unit"
    return out.to_numpy().astype(np.int64), rec.pos.to_numpy()


def source(cls: str, axis: str) -> tuple[str, str]:
    if cls in ("snp", "smallindel"):
        return f"{REC}/lfmm_{cls}_gen9_{axis}.csv", "record"
    return f"{EVT}/lfmm_event_{cls}_gen9_{axis}.csv.gz", "event"


def one(cls: str, axis: str, outdir: str, merge: bool) -> dict:
    src, level = source(cls, axis)
    d = pd.read_csv(src).drop(columns=["block"], errors="ignore")
    if cls == "smallindel":
        u, upos = smallindel_unit_ids()
        assert len(u) == len(d) and np.array_equal(upos, d.pos.to_numpy()), \
            "smallindel LFMM table is not row-aligned with class_matrices_cr50 records.csv"
        d["unit"] = u
    b = bt.assign_tiling(d.chrom.to_numpy(str), d.pos.to_numpy(np.int64), r2=0.9)
    n_raw = pd.unique(b).size
    if merge:
        b = bt.merge_small_blocks(b)
    d["block"] = b
    assert (d.block.to_numpy() != "").all(), f"{cls} {axis}: tiling left a record unassigned"
    cols = HEAD + [c for c in EXTRA if c in d.columns]
    os.makedirs(outdir, exist_ok=True)
    out = f"{outdir}/lfmm_{cls}_gen9_{axis}.csv"
    d[cols].to_csv(out, index=False)
    span = d.groupby("block").pos.agg(lambda s: s.max() - s.min())
    thr = 0.05 / len(d)
    r = dict(cls=cls, axis=axis, level=level, n=len(d), blocks_raw=n_raw,
             blocks=int(d.block.nunique()), median_span=float(span.median()),
             p99_span=float(span.quantile(0.99)), bonf=thr, hits=int((d.pval < thr).sum()),
             sig_blocks=int(d.loc[d.pval < thr, "block"].nunique()))
    print(f"  {cls:11s} {axis:5s} [{level}] {len(d):>9,} rows | blocks {n_raw:>7,}"
          f"{' -> ' + format(r['blocks'], ',') if merge else ''}"
          f" | median span {r['median_span']:>7,.0f} bp | Bonf {thr:.2e}: "
          f"{r['hits']:>5,} hits in {r['sig_blocks']:>4,} blocks", flush=True)
    return r


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", nargs="+", default=CLASSES, choices=CLASSES)
    ap.add_argument("--axes", nargs="+", default=AXES)
    ap.add_argument("--outdir", default=OUTDIR)
    ap.add_argument("--merge", action="store_true",
                    help="apply merge_small_blocks per class, as the archived "
                         "wza_in_clq09_tile tables did. OFF by default -- it makes block ids "
                         "mean different intervals in different classes (see module docstring)")
    a = ap.parse_args()
    lib.check_version_dir(a.outdir, "--outdir")         # never mix panel versions (KMATE_PANEL)
    rows = [one(c, ax, a.outdir, a.merge) for c in a.classes for ax in a.axes]
    S = pd.DataFrame(rows)
    f = f"{a.outdir}/_reblock_summary.csv"
    if os.path.exists(f):                                   # upsert: one row per (class, axis), newest wins
        S = pd.concat([pd.read_csv(f), S], ignore_index=True).drop_duplicates(["cls", "axis"], keep="last")
    S.to_csv(f, index=False)
    print(f"\nwrote {len(rows)} tables to {a.outdir}  (summary updated in {f}: {len(S)} rows)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""The production GEA tables: both arms in one schema, with tiling LD blocks (step 5 of _submit_gea.sh).

  snp, smallindel   record arm, results/multiaxis/lfmm_cr50_arch3v2/ (ma_lfmm_cr50.sbatch). A SNP or short indel
                    is one record, so record == unit; short indels get their integer event-unit id (`unit`).
  sv, nonsnp        event arm, genes/convergence/event_gea/cr50_arch3v2/ (run_event_lfmm_cr50.sbatch).

Block = the clq r2 = 0.9 tiling partition (blocks/blocks_tiling.assign_tiling): defined once on dense common
variants, then every variant of every class falls into the interval that contains it, so a block id names the same
stretch of LD in every class. Blocks are not merged per class (merge_small_blocks would make the same id mean
different intervals in different classes).

Writes results/multiaxis/tile_cr50_arch3v2/lfmm_{cls}_gen9_{axis}.csv: chrom,pos,ref_len,alt_len,MAF,block,pval,
plus unit,dlen,n_alleles,is_sv,n_called,B,z where the arm has them; and _reblock_summary.csv (rows, blocks,
Bonferroni hits per class x axis).
Environment: kmate.
"""
from __future__ import annotations
import argparse, functools, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GS = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, GS)
sys.path.insert(0, os.path.join(GS, "blocks"))
import lib                                                        # noqa: E402
import blocks_tiling as bt                                        # noqa: E402

MA = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis"
REC = f"{MA}/lfmm_cr50{lib.VSUF}"
EVT = f"{lib.GEA}/genes/convergence/event_gea/cr50{lib.VSUF}"
CM = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices_cr50{lib.VSUF}"
OUTDIR = f"{MA}/tile_cr50{lib.VSUF}"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1", "pc2", "pc3"]
CLASSES = ["snp", "smallindel", "sv", "nonsnp"]
HEAD = ["chrom", "pos", "ref_len", "alt_len", "MAF", "block", "pval"]
EXTRA = ["unit", "dlen", "n_alleles", "is_sv", "n_called", "B", "z"]


@functools.lru_cache(maxsize=1)
def smallindel_unit_ids() -> tuple[np.ndarray, np.ndarray]:
    """(unit id, pos) per row of the record-arm small-indel tables.

    Matched on integers: records.csv `col` is the record's index in af_store/index_nonsnp, and nonsnp_event_units.npz
    lists each unit's members as those same indices. Lets downstream code collapse on `unit` rather than on
    (chrom, pos, ref_len, alt_len), which is not unique.
    """
    E = np.load(f"{lib.SV_ADAPTIVE}/nonsnp_event_units.npz")
    U = pd.read_csv(f"{EVT}/nonsnp_gen9_units.csv", usecols=["unit", "is_sv", "n_alleles"])
    N = U[~U.is_sv]
    assert (N.n_alleles == 1).all(), "a non-SV unit holds >1 record; the 1:1 mapping is invalid"
    unit_of = pd.Series(N.unit.to_numpy(), index=E["indices"][E["indptr"][N.unit.to_numpy()]])
    rec = pd.read_csv(f"{CM}/smallindel_gen9.records.csv", usecols=["pos", "col"])
    out = unit_of.reindex(rec.col.to_numpy())
    assert out.notna().all(), f"{int(out.isna().sum()):,} small-indel records have no unit"
    return out.to_numpy().astype(np.int64), rec.pos.to_numpy()


def one(cls: str, axis: str, outdir: str) -> dict:
    if cls in ("snp", "smallindel"):
        level, src = "record", f"{REC}/lfmm_{cls}_gen9_{axis}.csv"
    else:
        level, src = "event", f"{EVT}/lfmm_event_{cls}_gen9_{axis}.csv.gz"
    d = pd.read_csv(src).drop(columns=["block"], errors="ignore")
    if cls == "smallindel":
        u, upos = smallindel_unit_ids()
        assert len(u) == len(d) and np.array_equal(upos, d.pos.to_numpy()), \
            "smallindel LFMM table is not row-aligned with the class matrix records.csv"
        d["unit"] = u
    d["block"] = bt.assign_tiling(d.chrom.to_numpy(str), d.pos.to_numpy(np.int64), r2=0.9)
    assert (d.block.to_numpy() != "").all(), f"{cls} {axis}: tiling left a record unassigned"
    d[HEAD + [c for c in EXTRA if c in d.columns]].to_csv(f"{outdir}/lfmm_{cls}_gen9_{axis}.csv", index=False)

    span = d.groupby("block").pos.agg(lambda s: s.max() - s.min())
    thr = 0.05 / len(d)
    r = dict(cls=cls, axis=axis, level=level, n=len(d), blocks=int(d.block.nunique()),
             median_span=float(span.median()), p99_span=float(span.quantile(0.99)), bonf=thr,
             hits=int((d.pval < thr).sum()), sig_blocks=int(d.loc[d.pval < thr, "block"].nunique()))
    print(f"  {cls:11s} {axis:5s} [{level}] {len(d):>9,} rows | {r['blocks']:>7,} blocks | Bonferroni {thr:.2e}: "
          f"{r['hits']:>5,} hits in {r['sig_blocks']:>4,} blocks", flush=True)
    return r


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", nargs="+", default=CLASSES, choices=CLASSES)
    ap.add_argument("--axes", nargs="+", default=AXES)
    ap.add_argument("--outdir", default=OUTDIR)
    a = ap.parse_args()
    lib.check_version_dir(a.outdir, "--outdir")             # never mix panel versions (KMATE_PANEL)
    os.makedirs(a.outdir, exist_ok=True)
    S = pd.DataFrame([one(c, ax, a.outdir) for c in a.classes for ax in a.axes])
    f = f"{a.outdir}/_reblock_summary.csv"
    if os.path.exists(f):                                   # one row per (class, axis), newest wins
        S = pd.concat([pd.read_csv(f), S], ignore_index=True).drop_duplicates(["cls", "axis"], keep="last")
    S.drop(columns=["blocks_raw"], errors="ignore").to_csv(f, index=False)
    print(f"\nwrote {len(a.classes) * len(a.axes)} tables to {a.outdir} (summary: {f}, {len(S)} rows)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Summarise the clq0.9 WZA outputs: BH/Bonferroni per output, cross-class &
cross-model recurrence, and the CAM5 headline block under the finer blocks.

clq0.9 replication naming: wza_{model}_{snp|nonsnp}_gen9_bio1_{deg2|deg7cap2000}.csv
(deg2 = canonical Booker PRIMARY; deg7cap2000 = phase-1-style sensitivity).

CAM5 (AT2G27030, Chr2 ~11.53 Mb) is mapped to its clq0.9 block(s) by intersecting
the TAIR10 gene span with the Chr2 clq0.9 interval table (the same blocks the WZA
windows are defined on), so the ids match.

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY compare_clq90.py --deg deg2            # primary
  $PY compare_clq90.py --deg both
"""
from __future__ import annotations
import argparse, glob, os, re, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PCOL = "Z_pVal"
_NAME_RE = re.compile(
    r"wza_(?P<model>\w+?)_(?P<cls>snp|nonsnp)_gen(?P<gen>\d+)_"
    r"(?P<clim>bio\d+)_(?P<regime>deg2|deg7cap2000)\.csv$")


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float); n = len(p); order = np.argsort(p); q = np.empty(n)
    q[order] = (p[order] * n) / (np.arange(n) + 1)
    q[order] = np.minimum.accumulate(q[order][::-1])[::-1]
    return np.clip(q, 0, 1)


def cam5_clq_blocks(gene_id: str, r2: float) -> tuple[list[str], str, int, int]:
    """clq{r2} block ids whose [start,end] intersect the gene span."""
    g = lib.load_genes(); row = g[g.gene == gene_id]
    if row.empty:
        raise SystemExit(f"gene {gene_id} not in TAIR10 annotation")
    r = row.iloc[0]; ci = int(str(r.chrom).replace("Chr", ""))
    tag = f"clq{r2}"
    bt = pd.read_csv(f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv", sep="\t")
    bt = bt.sort_values("start_pos").reset_index(drop=True)
    hit = bt[(bt.start_pos <= r.end) & (bt.end_pos >= r.start)]
    blks = [f"Chr{ci}_{i}" for i in hit.index]
    return blks, r.chrom, int(r.start), int(r.end)


def load_one(path: str, fdr: float) -> dict | None:
    m = _NAME_RE.search(os.path.basename(path))
    if not m:
        return None
    w = pd.read_csv(path)
    gcol = "gene" if "gene" in w.columns else w.columns[0]
    w = w.rename(columns={gcol: "block"})
    w["block"] = w["block"].astype(str)
    if PCOL not in w.columns:
        return None
    w = w[w[PCOL].notna()].copy()
    w["q"] = bh_fdr(w[PCOL].to_numpy())
    d = m.groupdict()
    d.update(path=path, df=w, n=len(w), n_bh=int((w.q < fdr).sum()),
             n_bonf=int((w[PCOL] < 0.05 / len(w)).sum()))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wza", default=f"{lib.GEA}/phase1_replication/results/clq90/wza")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--gene", default="AT2G27030", help="target gene (default CAM5)")
    ap.add_argument("--r2", type=float, default=0.9)
    ap.add_argument("--deg", default="deg2", choices=["deg2", "deg7cap2000", "both"])
    ap.add_argument("--out", default=f"{lib.GEA}/phase1_replication/results/clq90/compare_clq90.csv")
    args = ap.parse_args()

    paths = sorted(glob.glob(f"{args.wza}/wza_*.csv"))
    outs = [o for o in (load_one(p, args.fdr) for p in paths) if o]
    if args.deg != "both":
        outs = [o for o in outs if o["regime"] == args.deg]
    if not outs:
        sys.exit(f"no usable WZA outputs in {args.wza} (deg={args.deg})")

    tgt, tchrom, tstart, tend = cam5_clq_blocks(args.gene, args.r2)
    print(f"\n=== target {args.gene} -> clq{args.r2} block(s) {tgt} "
          f"({tchrom}:{tstart:,}-{tend:,}) ===\n")

    rows = []
    for o in sorted(outs, key=lambda d: (d["cls"], d["model"], d["regime"])):
        w = o["df"].sort_values(PCOL).reset_index(drop=True)
        hit = w.index[w.block.isin(tgt)]
        if len(hit):
            i = int(hit[0])
            cam = (f"blk {w.loc[i,'block']} rank {i+1}/{len(w)} p={w.loc[i,PCOL]:.2e} "
                   f"q={w.loc[i,'q']:.3f}" + ("  *BH*" if w.loc[i,'q'] < args.fdr else ""))
        else:
            cam = "block absent"
        print(f"{o['model']:9s} {o['cls']:7s} {o['regime']:12s} | {o['n']:6,} blocks | "
              f"BH<{args.fdr}: {o['n_bh']:4d} | Bonf: {o['n_bonf']:3d} | CAM5: {cam}")
        rows.append(dict(model=o["model"], cls=o["cls"], regime=o["regime"],
                         n_blocks=o["n"], n_bh=o["n_bh"], n_bonf=o["n_bonf"], cam5=cam))

    # cross-class + cross-model recurrence of BH-sig blocks (per regime)
    for regime in sorted({o["regime"] for o in outs}):
        print(f"\n=== BH<{args.fdr} block recurrence ({regime}) ===")
        sig = {(o["model"], o["cls"]): set(o["df"].loc[o["df"].q < args.fdr, "block"])
               for o in outs if o["regime"] == regime}
        for model in sorted({k[0] for k in sig}):
            cls_sets = {c: s for (mo, c), s in sig.items() if mo == model}
            if len(cls_sets) >= 2:
                shared = set.intersection(*cls_sets.values())
                print(f"  {model}: snp({len(cls_sets.get('snp',[]))}) ∩ "
                      f"nonsnp({len(cls_sets.get('nonsnp',[]))}) = {len(shared)} shared "
                      f"{sorted(shared)[:10]}")
        for cls in ("snp", "nonsnp"):
            msets = {m: s for (m, c), s in sig.items() if c == cls}
            if len(msets) >= 2:
                allm = set.intersection(*msets.values())
                print(f"  {cls}: blocks BH-sig in ALL models {sorted(msets)}: "
                      f"{len(allm)} {sorted(allm)[:10]}")

    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(f"\n  summary -> {args.out}")


if __name__ == "__main__":
    main()

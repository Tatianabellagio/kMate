#!/usr/bin/env python
"""Intersect WZA block-level hits across class x gen x model, and check CAM5.

Reads every `wza_*.csv` produced by run_wza.py (window col `gene` = LD block id,
`Z_pVal` = empirical block p, retained `chrom`/`pos`), and for each:
  * BH-FDR q-values and Bonferroni threshold over the block p-values,
  * the set of significant blocks (BH q < --fdr, default 0.05),
then reports:
  1. a per-(model,class,gen) summary (n blocks, n BH-sig, n Bonferroni, top block),
  2. which blocks recur across generations (within a class) and across classes,
  3. the CAM5 (AT2G27030, Chr2 ~11.53 Mb) block: its rank and p in every output —
     the phase-1 headline we are trying to reproduce.

CAM5's LD block is found by mapping the TAIR10 gene midpoint through the SAME
lib.assign_ld_blocks used to build the class matrices, so the id matches exactly.

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY compare_blocks.py                       # all wza_*.csv in the wza dir
  $PY compare_blocks.py --fdr 0.1 --gene AT2G27030
"""
from __future__ import annotations
import argparse, glob, os, re, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

PCOL = "Z_pVal"
# wza_{model}_{cls}_gen{g}_{clim}[_{regime}].csv ; regime e.g. deg7nocap / deg7cap2000 /
# deg2nocap (absent = canonical deg-2). clim is always bioNN (no underscore).
_NAME_RE = re.compile(
    r"wza_(?P<model>\w+?)_(?P<cls>snp|sv|smallindel)_gen(?P<gen>\d+)_"
    r"(?P<clim>bio\d+)(?:_(?P<regime>[a-z0-9]+))?\.csv$")


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = (p[order] * n) / (np.arange(n) + 1)
    # enforce monotonicity from the largest p downward
    q[order] = np.minimum.accumulate(q[order][::-1])[::-1]
    return np.clip(q, 0, 1)


def cam5_block(gene_id: str) -> tuple[list[str], str, int]:
    """Return (block_ids, chrom, midpoint) for ALL LD blocks overlapping `gene_id`.

    A gene can span >1 LD block; phase-1 placed CAM5 in block 2_1265 (its 3' end),
    while the gene midpoint maps to 2_1264 — so a midpoint-only lookup misses the
    real signal. We take every block whose mapped LD-SNP falls in [start,end].
    """
    g = lib.load_genes()
    row = g[g.gene == gene_id]
    if row.empty:
        raise SystemExit(f"gene {gene_id} not in TAIR10 annotation")
    r = row.iloc[0]
    mid = int((r.start + r.end) // 2)
    # scan the phase-1 LD-block SNP map for blocks intersecting the gene span
    ld = pd.read_csv(lib.LD_BLOCKS, usecols=["pos", "chrom", "block"])
    ci = r.chrom.replace("Chr", "")
    sub = ld[(ld.chrom.astype(str) == ci) & (ld.pos >= r.start) & (ld.pos <= r.end)]
    blks = [str(b) for b in pd.unique(sub.block)]
    if not blks:  # gene falls between mapped SNPs -> nearest-SNP block of the midpoint
        blks = [str(lib.assign_ld_blocks(np.array([r.chrom]), np.array([mid]))[0])]
    return blks, r.chrom, mid


def load_one(path: str, fdr: float) -> dict | None:
    m = _NAME_RE.search(os.path.basename(path))
    if not m:
        print(f"  (skip unrecognized name) {os.path.basename(path)}")
        return None
    w = pd.read_csv(path)
    gcol = "gene" if "gene" in w.columns else w.columns[0]
    w = w.rename(columns={gcol: "block"})
    w["block"] = w["block"].astype(str).str.replace(r"\.0$", "", regex=True)
    if PCOL not in w.columns:
        print(f"  (no {PCOL}, raw scores only) {os.path.basename(path)}")
        return None
    w = w[w[PCOL].notna()].copy()
    w["q"] = bh_fdr(w[PCOL].to_numpy())
    w["bonf"] = w[PCOL] < 0.05 / len(w)
    d = m.groupdict()
    regime = d.pop("regime") or "deg2"          # absent suffix = canonical deg-2
    deg = "deg7" if regime.startswith("deg7") else "deg2"
    d.update(path=path, regime=regime, deg=deg, df=w,
             n=len(w), n_bh=int((w.q < fdr).sum()), n_bonf=int(w.bonf.sum()))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wza", default=f"{lib.GEA}/phase1_replication/results/wza")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--gene", default="AT2G27030", help="target gene (default CAM5)")
    ap.add_argument("--deg", default="deg2", choices=["deg2", "deg7", "both"],
                    help="which correction to compare on (default deg2 primary)")
    ap.add_argument("--out", default=None, help="optional CSV of the per-output summary")
    args = ap.parse_args()

    paths = sorted(glob.glob(f"{args.wza}/wza_*.csv"))
    if not paths:
        sys.exit(f"no wza_*.csv in {args.wza} — run run_wza.py first")
    outs = [o for o in (load_one(p, args.fdr) for p in paths) if o]
    if args.deg != "both":
        outs = [o for o in outs if o["deg"] == args.deg]
    if not outs:
        sys.exit("no usable WZA outputs after filtering")

    tgt_blocks, tgt_chrom, tgt_mid = cam5_block(args.gene)
    print(f"\n=== target {args.gene} -> LD block(s) {tgt_blocks} ({tgt_chrom}:{tgt_mid:,}) ===")
    print("    (CAM5 line reports the BEST-ranked of these overlapping blocks)\n")

    # 1) per-output summary + CAM5 rank (best overlapping block)
    rows = []
    for o in sorted(outs, key=lambda d: (d["cls"], d["gen"], d["model"], d["regime"])):
        w = o["df"].sort_values(PCOL).reset_index(drop=True)
        hit = w.index[w.block.isin(tgt_blocks)]
        if len(hit):
            i = int(hit[0]); pr = w.loc[i, PCOL]; qr = w.loc[i, "q"]; bb = w.loc[i, "block"]
            cam = f"blk {bb} rank {i+1}/{len(w)}  p={pr:.2e} q={qr:.3f}" + ("  *BH*" if qr < args.fdr else "")
        else:
            cam = "block absent"
        tag = f"{o['model']:9s} {o['cls']:11s} gen{o['gen']} {o['regime']:12s}"
        print(f"{tag} | {o['n']:6,} blocks | BH<{args.fdr}: {o['n_bh']:4d} | "
              f"Bonf: {o['n_bonf']:3d} | CAM5: {cam}")
        rows.append(dict(model=o["model"], cls=o["cls"], gen=o["gen"], regime=o["regime"],
                         deg=o["deg"], n_blocks=o["n"], n_bh=o["n_bh"], n_bonf=o["n_bonf"],
                         cam5_in="block absent" not in cam, cam5=cam))

    # 2) recurrence of BH-significant blocks across gens (within class) and classes
    print(f"\n=== BH-significant block recurrence ({args.deg}) ===")
    sig = {}
    for o in outs:
        if o["deg"] != ("deg2" if args.deg == "both" else args.deg):
            continue
        s = set(o["df"].loc[o["df"].q < args.fdr, "block"])
        sig[(o["model"], o["cls"], o["gen"])] = s
    # within class, across gens (same model)
    for model in sorted({k[0] for k in sig}):
        for cls in sorted({k[1] for k in sig if k[0] == model}):
            gens = sorted(g for (mo, cl, g) in sig if mo == model and cl == cls)
            if len(gens) >= 2:
                inter = set.intersection(*(sig[(model, cls, g)] for g in gens))
                print(f"  {model} {cls}: blocks BH-sig in ALL gens {gens}: "
                      f"{len(inter)}  {sorted(inter)[:12]}{' ...' if len(inter)>12 else ''}")
    # across classes (per model, union over gens)
    for model in sorted({k[0] for k in sig}):
        per_cls = {}
        for (mo, cl, g), s in sig.items():
            if mo == model:
                per_cls.setdefault(cl, set()).update(s)
        if len(per_cls) >= 2:
            shared = set.intersection(*per_cls.values())
            print(f"  {model} shared across classes {sorted(per_cls)}: "
                  f"{len(shared)}  {sorted(shared)[:12]}{' ...' if len(shared)>12 else ''}")

    if args.out:
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print(f"\n  summary -> {args.out}")


if __name__ == "__main__":
    main()

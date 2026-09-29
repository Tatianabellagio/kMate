#!/usr/bin/env python
"""Where does a named gene set sit in the raw GEA -- every class, every axis?

The convergence pool only holds Bonferroni-significant NON-SNP records, so a gene absent
from it may still carry sub-threshold signal, or SNP-only signal. This reads the raw
per-record LFMM tables directly (wza_in_clq09_tile, all classes x 22 axes), takes every
record in gene body +-FLANK, and reports per gene x class: best -log10 p, the axis it
came from, lambda of that axis, and whether it clears that file's own Bonferroni
(-log10(0.05 / n_records), computed per class x axis exactly as the GEA pool did).

Raw p, no inflation correction (memory `raw-lfmm-over-wza-decision`).

Usage: PY scan_gene_set.py AGI[:label] ...     writes results/gene_set_scan.csv
env: kmate. Compute node, several minutes (reads the full SNP tables).
"""
from __future__ import annotations
import os, sys, glob
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                       # noqa: E402

W = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile"
FLANK = 2000


def main(args):
    lab = {a.split(":")[0]: (a.split(":")[1] if ":" in a else a) for a in args}
    G = lib.load_genes()
    G = G[G.gene.isin(lab)].copy()
    G["lo"], G["hi"] = G.start - FLANK, G.end + FLANK
    missing = set(lab) - set(G.gene)
    if missing:
        print("no TAIR10 coordinates:", sorted(missing))
    rows = []
    for f in sorted(glob.glob(f"{W}/lfmm_*_gen9_*.csv")):
        base = os.path.basename(f)[5:-4]              # "<cls>_gen9_<axis>"
        cls, axis = base.split("_gen9_")
        d = pd.read_csv(f, usecols=["chrom", "pos", "MAF", "pval"])
        d = d[d.MAF > 0.05]
        bonf = -np.log10(0.05 / len(d))
        # lambda for context
        from scipy.stats import chi2
        lam = float(np.median(chi2.isf(d.pval.clip(1e-300), 1)) / chi2.ppf(0.5, 1))
        for g in G.itertuples():
            s = d[(d.chrom == g.chrom) & (d.pos >= g.lo) & (d.pos <= g.hi)]
            if not len(s):
                continue
            nlp = -np.log10(s.pval.clip(1e-300))
            rows.append(dict(gene=g.gene, symbol=lab[g.gene], cls=cls, axis=axis,
                             n_records=len(s), best_nlp=float(nlp.max()),
                             best_pos=int(s.pos.values[nlp.values.argmax()]),
                             n_bonf=int((nlp > bonf).sum()), bonf=round(bonf, 2),
                             lam=round(lam, 2)))
        print(f"  {cls:<10} {axis:<6} done", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(f"{HERE}/results/gene_set_scan.csv", index=False)

    print("\n=== per gene x class: best record over all 22 axes ===")
    best = (R.sort_values("best_nlp", ascending=False)
             .groupby(["symbol", "cls"]).head(1)
             .sort_values(["symbol", "cls"]))
    sig = (R[R.n_bonf > 0].groupby(["symbol", "cls"]).axis
             .apply(lambda x: ",".join(sorted(set(x)))))
    best["bonf_axes"] = [sig.get((a, b), "") for a, b in zip(best.symbol, best.cls)]
    print(best[["symbol", "gene", "cls", "best_nlp", "axis", "lam", "bonf",
                "bonf_axes"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1:])
